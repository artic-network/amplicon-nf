#!/usr/bin/env python3

import itertools
from types import SimpleNamespace
import pysam
import sys
import os
from collections import defaultdict
from importlib.metadata import version


# from https://www.geeksforgeeks.org/python-make-a-list-of-intervals-with-sequential-numbers/
# via artic-mask
def intervals_extract(iterable):
    iterable = sorted(set(iterable))
    for _, group in itertools.groupby(enumerate(iterable), lambda t: t[1] - t[0]):
        group = list(group)
        yield [group[0][1], group[-1][1]]

# write the depth mask used with bcftools to turn consensus positions into Ns
def write_depth_mask(out_filename, contig_depths, min_coverage):
    maskfh = open(out_filename, 'w')
    for contig_name, depths in contig_depths.items():
        # from artic-mask, create list of positions that fail the depth check
        mask_vector = []
        for pos, depth in enumerate(depths):
            if depth < min_coverage:
                mask_vector.append(pos)

        # get the intervals from the mask_vector
        intervals = list(intervals_extract(mask_vector))

        for i in intervals:
            maskfh.write("%s\\t%s\\t%s\\n" % (contig_name, i[0]+1, i[1]+1))
    maskfh.close()

# calculate the variant allele fraction for each alt allele using freebayes' read/alt observation tags
def calculate_vafs(record):
    vafs = list()
    total_depth = float(record.info["DP"])
    for i in range(0, len(record.alts)):
        alt_reads = int(record.info["AO"][i])
        vaf = float(alt_reads) / float(record.info["DP"])
        vafs.append(vaf)
    return vafs

# make a simple VCF record with the minimal information needed to make the consensus sequence
def make_simple_record(vcf_header, parent_record, position, ref, alt, vaf):
    r = vcf_header.new_record()
    r.chrom = parent_record.chrom
    r.pos = position
    r.ref = ref
    r.alts = [ alt ]
    r.info["DP"] = parent_record.info["DP"]
    r.info["VAF"] = vaf
    return r

# process indel variants found by freebayes into a variant that should be
# applied to the consensus sequence
def handle_indel(vcf_header, record):
    output = list()
    vafs = calculate_vafs(record)

    # special case, if we have evidence for multiple possible indels (eg CTTT -> C, CTTT -> CT)
    # we decide whether to apply an indel based on the summed VAF across all alt alleles, then
    # apply the most frequent ALT. This is because there is evidence for /an/ indel but it is
    # ambiguous which one. We can't represent ambiguous indels in a consensus fasta so this
    # is the best we can do.
    if sum(vafs) < 0.5:
        return output

    # argmax without bringing in numpy
    max_idx = None
    max_vaf = 0.0
    for idx, value in enumerate(vafs):
        if value > max_vaf:
            max_vaf = value
            max_idx = idx

    r = make_simple_record(vcf_header, record, record.pos, record.ref, record.alts[max_idx], [ max_vaf ])
    output.append(r)
    return output

# return the base with the highest value in vaf_by_base,
# optionally skipping a character (eg. the reference base)
def base_max(vaf_by_base, skip=None):
    max_vaf = 0.0
    max_b = None
    for b in "ACGT":
        if b != skip and vaf_by_base[b] > max_vaf:
            max_vaf = vaf_by_base[b]
            max_b = b
    return max_b

def handle_sub(vcf_header, record):
    output = list()

    # this code is general enough to handle multi-allelic MNPs
    # and the typical case of a biallelic SNP
    sub_length = len(record.ref)

    vafs = calculate_vafs(record)

    # calculate the VAF of each base at each position of the MNP
    base_frequency = list()
    for i in range(0, sub_length):
        base_frequency.append({ "A":0.0, "C":0.0, "G":0.0, "T":0.0 })

    for alt, vaf in zip(record.alts, vafs):
        assert(len(alt) == sub_length)
        for i,b in enumerate(alt):
            base_frequency[i][b] += vaf

    # construct output records
    for i in range(0, sub_length):

        # choose base with highest frequency, skipping the reference
        max_b = base_max(base_frequency[i], record.ref[i])
        if max_b is None:
            continue
        r = make_simple_record(vcf_header, record, record.pos + i, record.ref[i], max_b, base_frequency[i][max_b])
        output.append(r)
    return output

def main():

    args = SimpleNamespace(mask_output="mask.txt", variants_output="variants.vcf", consensus_sites_output="consensus.vcf",
                            min_depth=${params.min_coverage_depth}, lower_ambiguity_frequency=${params.lower_ambiguity_frequency}, upper_ambiguity_frequency=${params.upper_ambiguity_frequency})

    vcf = pysam.VariantFile(open("${gvcf}","r"))

    # Initalize depth mask to all zeros for all contigs
    contig_depth = defaultdict(list)
    for r in vcf.header.records:
        if r.type == "CONTIG":
            contig_depth[r['ID']] = [0] * int(r['length'])

    out_header = vcf.header

    # open the output file with the filtered variant sites
    out_header.info.add("VAF", number="A", type='Float', description="Variant allele fraction, called from observed reference/alt reads")
    variants_out = pysam.VariantFile(args.variants_output, 'w', header=out_header)

    # open the output file with the changes to apply to the consensus fasta
    # this includes an additional tag in the VCF file
    out_header.info.add("ConsensusTag", number=1, type='String', description="The type of base to be included in the consensus sequence (IUPAC or Fixed)")
    consensus_sites_out = pysam.VariantFile(args.consensus_sites_output, 'w', header=out_header)

    for record in vcf:

        is_gvcf_ref = record.alts[0] == "<*>"

        # set depth for this part of the genome
        # this works for both gVCF blocks and regular variants
        # because pos/stop are set appropriately
        v_start = record.pos
        v_end = record.stop
        depth = record.info["DP"]

        # disallow gvcf records that are longer than a single base
        assert(not is_gvcf_ref or v_start == v_end)

        # update depth mask
        for i in range(v_start, v_end + 1):
            assert(i > 0)
            # VCF coordinates are 1-based, we record the depth vector as 0-based
            # to be consistent with artic-mask
            contig_depth[record.chrom][i - 1] = depth

        # do nothing else with ref records, or records that don't meet our minimum depth
        if is_gvcf_ref or depth < args.min_depth:
            continue

        # determine if any allele in the variant is an indel
        has_indel = False
        for i in range(0, len(record.alts)):
            has_indel = has_indel or len(record.ref) != len(record.alts[i])

        # process the input variant record to handle multi-allelic variants and MNPs
        out_records = list()
        if has_indel:
            # indels need to be handle specially as we can't apply ambiguity codes
            out_records = handle_indel(out_header, record)
        else:
            out_records = handle_sub(out_header, record)

        # classify variants using VAF cutoffs for IUPAC ambiguity codes, etc
        accept_variant = False
        for out_r in out_records:

            # at this point we should have resolved multi-allelic variants
            assert(len(out_r.alts) == 1)

            vaf = out_r.info["VAF"][0]
            is_indel = len(out_r.ref) != len(out_r.alts[0])

            # discard low frequency variants
            if vaf < args.lower_ambiguity_frequency:
                continue

            # Write a tag describing what to do with the variant
            consensus_tag = "None"

            # high-frequency subs and indels are always applied without ambiguity
            # we don't have to do an indel VAF check here as it is dealt with in handle_indel
            if vaf > args.upper_ambiguity_frequency or is_indel:
                # always apply these to the consensus
                consensus_tag = "fixed"
            else:
                # record ambiguous SNPs in the consensus sequence with IUPAC codes
                consensus_tag = "ambiguous"
            out_r.info["ConsensusTag"] = consensus_tag
            consensus_sites_out.write(out_r)
            accept_variant = True

        if accept_variant:
            record.info["VAF"] = calculate_vafs(record)
            variants_out.write(record)

    write_depth_mask(args.mask_output, contig_depth, args.min_depth)

    with open("versions.yml", "wt") as versions_fh:
        versions_fh.write("${task.process}:\\n  pysam: " + str(version("pysam")))


main()
