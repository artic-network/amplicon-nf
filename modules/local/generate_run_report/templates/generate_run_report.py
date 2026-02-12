#!/usr/bin/env python3

import plotly.graph_objects as go
import plotly.io as pio
import plotly.express as px

import numpy as np
import pathlib
from importlib.metadata import version

from primalbedtools.scheme import Scheme
from primalbedtools.amplicons import create_amplicons
from primalbedtools.bedfiles import BedLine

from collections import Counter
import pandas as pd
from Bio import SeqIO
import csv
from enum import Enum
from glob import glob
import math

from pathlib import Path
from datetime import datetime
from jinja2 import Environment, FileSystemLoader, select_autoescape

ALL_DNA = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
}

AMBIGUOUS_DNA_COMPLEMENT = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "M": "K",
    "R": "Y",
    "W": "W",
    "S": "S",
    "Y": "R",
    "K": "M",
    "V": "B",
    "H": "D",
    "D": "H",
    "B": "V",
    "X": "X",
    "N": "N",
    "-": "-",
}

IUPAC_ALL_ALLOWED_DNA = {
    "A",
    "G",
    "K",
    "Y",
    "B",
    "S",
    "N",
    "H",
    "C",
    "W",
    "D",
    "R",
    "M",
    "T",
    "V",
    "-",
}

plotly_colour_scale = [(0, "#345c67"), (1, "#f3edca")]


class MappingType(Enum):
    """
    Enum for the mapping type
    """

    FIRST = "first"
    CONSENSUS = "consensus"


class PlotlyText:
    """
    A class to hold the text for a plotly heatmap.
    """

    primer_name: str
    primer_seq: str
    genome_seq: str

    def __init__(
        self,
        primer_name: str,
        primer_seq: str,
        genome_seq: str,
    ):
        self.primer_name = primer_name
        self.primer_seq = primer_seq
        self.genome_seq = genome_seq

    def format_str(self) -> str:
        # parsedseqs
        cigar = []
        for p, g in zip(self.primer_seq[::-1], self.genome_seq[::-1]):
            if p == g:
                cigar.append("|")
            else:
                cigar.append(".")
        cigar = "".join(cigar)[::-1]
        return f"5'{self.primer_seq}: {self.primer_name}<br>5'{cigar}<br>5'{self.genome_seq[-len(self.primer_seq) :]}"


def remove_end_insertion(msa_array: np.ndarray) -> np.ndarray:
    """
    Removes leading and trailing "-" from an msa
    """
    tmp_array = msa_array
    ncols = tmp_array.shape[1]
    for row_index in range(0, tmp_array.shape[0]):
        # Solves the 5' end
        for col_index in range(0, ncols):
            if tmp_array[row_index, col_index] == "-":
                tmp_array[row_index, col_index] = ""
            else:
                break
        for rev_col_index in range(ncols - 1, 0, -1):
            if tmp_array[row_index, rev_col_index] == "-":
                tmp_array[row_index, rev_col_index] = ""
            else:
                break
    return tmp_array


def parse_msa(msa_path: pathlib.Path) -> tuple[np.ndarray, dict]:
    """
    Parses a multiple sequence alignment (MSA) file in FASTA format.

    This function reads an MSA file, validates its format and content, and returns a numpy array of the sequences
    and a dictionary with additional information. It checks for sequences of different lengths, empty columns,
    and non-DNA characters. It also removes end insertions from the sequences.

    Args:
        msa_path (pathlib.Path): The path to the MSA file to be parsed.

    Returns:
        tuple: A tuple containing two elements:
            - np.ndarray: A 2D numpy array where each row represents a sequence in the MSA and each column represents a position in the alignment.
            - dict: A dictionary with additional information about the MSA (currently not implemented, returns an empty dict).

    Raises:
        MSAFileInvalidLength: If the MSA contains sequences of different lengths.
        MSAFileInvalid: If the MSA file is empty or not in FASTA format.
        ValueError: If the MSA contains empty columns.
        MSAFileInvalidBase: If the MSA contains non-DNA characters.
    """
    try:
        records_index = SeqIO.index(
            str(msa_path),
            "fasta",
        )
    except ValueError as e:
        raise Exception(f"{msa_path.name}: {e}") from e

    try:
        array = np.array(
            [record.seq.upper() for record in records_index.values()],
            dtype="U1",
            ndmin=2,  # Enforce 2D array even if one genome
        )
    except ValueError as e:
        raise Exception(
            f"MSA ({msa_path.name}): contains sequences of different lengths"
        ) from e

    # Check for empty MSA, caused by no records being parsed
    if array.size == 0:
        raise Exception(
            f"No sequences in MSA ({msa_path.name}). Please ensure the MSA uses .fasta format."
        )

    empty_set = {"", "-"}

    empty_col_indexes = []
    # Check for empty columns and non DNA characters
    for col_index in range(0, array.shape[1]):
        slice: set[str] = set(array[:, col_index])
        # Check for empty columns
        if slice.issubset(empty_set):
            empty_col_indexes.append(col_index)
        # Check for non DNA characters
        if slice.difference(IUPAC_ALL_ALLOWED_DNA):
            base_str = ", ".join(slice.difference(IUPAC_ALL_ALLOWED_DNA))
            raise Exception(
                f"MSA ({msa_path.name}) contains non DNA characters ({base_str}) at column: {col_index}"
            )
    # Remove empty columns
    array = np.delete(array, empty_col_indexes, axis=1)

    # Remove end insertions
    array = remove_end_insertion(array)

    return array, dict(records_index)


def extend_ambiguous_base(base: str) -> list[str]:
    """Return list of all possible sequences given an ambiguous DNA input"""
    return [*ALL_DNA.get(base, "N")]


def reverse_complement(kmer_seq: str) -> str:
    rev_seq = kmer_seq[::-1]
    return "".join(AMBIGUOUS_DNA_COMPLEMENT[base.upper()] for base in rev_seq)


def create_mapping(
    msa: np.ndarray, mapping_index: int = 0
) -> tuple[np.ndarray, np.ndarray]:
    """
    This returns a tuple of two items:
        mapping_array: list[int | None]
        truncated_msa: np.ndarray
    mapping_array: Each position in the list corresponds to the same index in the MSA, The value in the list is the position in the reference genome
    """
    # As NP is modified in place, returning is not necessary but is done for clarity
    # Create the empty mapping array
    mapping_list = [None] * msa.shape[1]
    mapping_array = np.array(mapping_list)
    # Select the reference genome
    reference_genome = msa[mapping_index]
    # Iterate over the msa genome
    current_ref_index = 0
    for col_index in range(msa.shape[1]):
        # If the base is not a gap, assign the mapping
        if reference_genome[col_index] not in {"", "-"}:
            mapping_array[col_index] = current_ref_index
            # increase reference index
            current_ref_index += 1
    return (mapping_array, msa)


def generate_consensus(msa: np.ndarray) -> str:
    """
    Generates a consensus sequence from an msa
    """
    consensus = []
    # For each column in the msa
    for col in range(msa.shape[1]):
        # Create the counter
        col_counter = Counter()
        # For each row in the msa
        for row in range(msa.shape[0]):
            # Update the counter with the de-ambiguous bases
            col_counter.update(extend_ambiguous_base(msa[row, col]))

        # Remove invalid bases if other bases are available
        col_counter.pop("N", None)

        if len(col_counter) == 0:
            consensus.append("N")
        else:
            consensus.append(col_counter.most_common(1)[0][0])
    return "".join(consensus)


def generate_reference(msa: np.ndarray) -> str:
    """
    Generates a reference string from the first row of a multiple sequence alignment (MSA) array.

    Args:
    - msa (np.ndarray): A numpy array representing a multiple sequence alignment.

    Returns:
    - str: A string representing the reference sequence, obtained by joining the characters in the first row of the MSA and removing any gaps represented by hyphens ("-").
    """

    return "".join(msa[0]).replace("-", "")


def ref_index_to_msa(mapping_array: np.ndarray) -> dict[int, int]:
    """
    Convert a reference index to an MSA index
    """
    ref_dict = {x: i for i, x in enumerate(list(mapping_array)) if x is not None}
    ref_dict[max(ref_dict.keys()) + 1] = (
        max(ref_dict.values()) + 1
    )  # This ensures that an fprimer with non-inclusive end will not cause key error.

    return ref_dict


def check_for_end_on_gap(ref_index_to_msa: dict[int, int], ref_index) -> bool:
    """
    Check if a slice of a mapping array ends on a gap
    Returns True if the slice ends on a gap, False otherwise

    # Example

           5' AGAGTGTGGGGGTAGTGTTACG          > MPXV_142_LEFT_1 170931:170953
    TTTTTTTTATAGAGTGTGGGGGTAGTGTTACG-------GAT >MT903345
    TTTTTTTTATAGAGTGT-GGGGTAGTGTTACGGATATCTGAT >KJ642613.1
                                           ^ 173598: MSA
                                           ^ 170953: ref

    In this example, the slice of the primer ends on a gap. So slicing the array with
    - array[:, ref_to_msa[170931]:ref_to_msa[170953]] will return "GGGGTAGTGTTACG-------" as the non exclusive end captures the gap
    - fix_end_on_gap() will return in indexes to slice the array without the gap ie "TGTGGGGGTAGTGTTACG"

    """
    exclusive_msa_end = ref_index_to_msa[ref_index]
    inclusive_msa_end = ref_index_to_msa[ref_index - 1]
    return exclusive_msa_end - inclusive_msa_end != 1


def fix_end_on_gap(ref_index_to_msa: dict[int, int], ref_index) -> int:
    """
    Returns the MSA index of the non-inclusive end of a slice with the gap removed
    """
    return ref_index_to_msa[ref_index - 1] + 1


def get_primers_from_msa(
    array: np.ndarray, index: int, forward: bool = True, length: int = 20, row=0
) -> dict[int, str | None]:
    """
    Get a primer from an MSA array.
    """
    row_data = {}
    if forward:
        for row in range(array.shape[0]):
            gaps = 0
            row_data[row] = None
            while index - length - gaps >= 0:
                # Get slice
                initial_slice = array[row, index - length - gaps : index]
                # Check for gaps on set base
                if initial_slice[-1] == "-":
                    break
                sequence = "".join(initial_slice).replace("-", "")
                # Covered removed gaps
                if "" in initial_slice:
                    break
                # Check for gaps in the slice
                if len(sequence) == length:
                    row_data[row] = sequence
                    break
                # Walk left
                gaps += 1
    else:
        for row in range(array.shape[0]):
            gaps = 0
            row_data[row] = None
            while index + length + gaps <= array.shape[1]:
                # Get slice
                initial_slice = array[row, index : index + length + gaps]
                # Check for gaps on set base
                if initial_slice[0] == "-":
                    break
                sequence = "".join(initial_slice).replace("-", "")
                # Covered removed gaps
                if "" in initial_slice:
                    break
                # Check for gaps in the slice
                if len(sequence) == length:
                    row_data[row] = reverse_complement(sequence)
                    break
                # Walk right
                gaps += 1
    return row_data


def calc_primer_hamming(seq1, seq2) -> int:
    """
    Calculate the hamming distance between two sequences of equal length. Ignores N.
    :param seq1: The primer sequence in 5' to 3' orientation.
    :param seq2: The primer sequence in 5' to 3' orientation.
    :return: The number of mismatches between the two sequences.
    """
    dif = 0
    for seq1b, seq2b in zip(seq1[::-1], seq2[::-1]):
        seq1b_exp = set(extend_ambiguous_base(seq1b))
        seq2b_exp = set(extend_ambiguous_base(seq2b))

        if not seq1b_exp & seq2b_exp and (seq1b != "N" and seq2b != "N"):
            dif += 1

    return dif


def render_qc_report(
    payload: dict,
    template_path: Path,
    output_path: Path,
    svg_logo_path: Path,
    bootstrap_css_path: Path,
    bootstrap_bundle_js_path: Path,
    plotly_js_path: Path,
):
    """
    Render a standalone HTML QC report using a Jinja2 template and a structured payload.

    Args:
        payload (dict): Input data with all required metadata and contig info.
        template_path (Path): Path to Jinja2 HTML template.
        output_path (Path): Path to write the output HTML.
        svg_logo_path (Path): Path to SVG logo file.
        bootstrap_css_path (Path): Path to Bootstrap 5 CSS.
        bootstrap_bundle_js_path (Path): Path to Bootstrap 5 JS (bundle version).
        plotly_js_path (Path): Path to Plotly JS.
    """

    payload = payload.copy()
    payload["embedded_logo_svg"] = Path(svg_logo_path).read_text(encoding="utf-8")
    payload["bootstrap_css"] = Path(bootstrap_css_path).read_text(encoding="utf-8")
    payload["bootstrap_bundle_js"] = Path(bootstrap_bundle_js_path).read_text(
        encoding="utf-8"
    )
    payload["plotly_js"] = Path(plotly_js_path).read_text(encoding="utf-8")

    env = Environment(
        loader=FileSystemLoader(template_path.parent),
        autoescape=select_autoescape(["html", "xml"]),
    )
    template = env.get_template(template_path.name)

    rendered_html = template.render(**payload)
    output_path.write_text(rendered_html, encoding="utf-8")
    print(f"[✓] QC report written to: {output_path}")


def amplicon_depth_heatmap(
    amplicon_depths: pd.DataFrame, scheme_str: str, chrom_name: str
) -> str:

    # Hovertemplate string
    # if include_seqs:
    #     hovertemplatestr = "%{text}<br>" + "<b>Mismatches: %{z}</b><br>"
    # else:
    #     hovertemplatestr = ""

    # amplicon_depths.sort_values(
    #     by="amplicon",
    #     inplace=True,
    #     key=lambda x: int(x.str.replace("Amplicon ", "")),
    # )
    # columns = [int(x) for x in amplicon_depths.amplicon.to_list()]

    amplicon_depths = amplicon_depths.pivot(index="sample", columns="amplicon")[
        "mean_depth"
    ].fillna(0)

    # amplicon_depths.reindex(columns=columns)

    fig = px.imshow(
        amplicon_depths,
        x=amplicon_depths.columns,
        y=amplicon_depths.index,
        color_continuous_scale=plotly_colour_scale,
        labels=dict(x="Amplicon", y="Sample", color="Mean Depth"),
        aspect="auto",
    )

    # # Create the heatmap
    # fig = go.Figure(
    #     data=go.Heatmap(
    #         z=amplicon_depths["mean_depth"],
    #         x=amplicon_depths["amplicon"],
    #         y=amplicon_depths["sample"],
    #         colorscale="Viridis",
    #         xgap=0.1,
    #         ygap=0.1,
    #         name="Amplicon Depths",
    #         hoverinfo="text",

    #         # labels=dict(x="Amplicon", y="Sample", z="Mean Depth"),
    #         legend=dict(x="Amplicon", y="Sample", z="Mean Depth"),
    #     )
    # )
    fig.update_layout(
        font=dict(family="Roboto, monospace", size=16),
        hoverlabel=dict(font_family="Roboto, monospace"),
        title_text=f"Amplicon Depths: {scheme_str} - {chrom_name}",
        coloraxis=dict(cmax=int("${params.min_coverage_depth}"), cmin=0),
    )
    # fig.update_yaxes(autorange="reversed")
    fig.update_xaxes(type="category")

    # Remove unnecessary plot elements
    fig.update_layout(
        modebar_remove=[
            "select2d",
            "lasso2d",
            "select",
        ]
    )

    return pio.to_html(
        fig,
        include_plotlyjs=False,
        full_html=False,
        default_height="3000px",
        default_width="100%",
    )


def primer_mismatch_heatmap(
    array: np.ndarray,
    seqdict: dict,
    bedfile: pathlib.Path,
    include_seqs: bool = True,
    mapping: MappingType = MappingType.FIRST,
) -> str:
    """
    Create a heatmap of primer mismatches in an MSA.
    :param array: The MSA array.
    :param seqdict: The sequence dictionary.
    :param bedfile: The bedfile of primers.
    :param include_seqs: Reduces plot size by removing hovertext.
    :raises: click.UsageError
    """
    # Read in the bedfile

    primer_scheme = Scheme.from_file(bedfile)

    scheme_headers = primer_scheme.header_dict

    # Find the mapping genome
    bed_chrom_names = {bedline.chrom for bedline in primer_scheme.bedlines}

    # Reference genome
    primary_ref = bed_chrom_names.intersection(seqdict.keys())

    if len(primary_ref) == 0:
        # Try to fix a common issue with Jalview
        parsed_seqdict = {"_".join(k.split("/")): v for k, v in seqdict.items()}
        primary_ref = bed_chrom_names.intersection(parsed_seqdict.keys())
        seqdict = parsed_seqdict

    # handle errors if mapping is set to first
    if len(primary_ref) == 0 and mapping == MappingType.FIRST:
        raise Exception(
            f"Primer chrom names ({', '.join(bed_chrom_names)}) not found in MSA ({', '.join(seqdict.keys())})"
        )
    # If consensus mapping ensure only one chrom in bedfile
    elif mapping == MappingType.CONSENSUS:
        primary_ref = ["Consensus"]
        if len(bed_chrom_names) > 1:
            raise Exception(
                f"Primer chrom names ({', '.join(bed_chrom_names)}) not found in MSA ({', '.join(seqdict.keys())})"
            )
    else:  # mapping == MappingType.FIRST & len(primary_ref) > 0
        # Filter the bedlines for only the reference genome
        bedlines = [
            bedline
            for bedline in primer_scheme.bedlines
            if bedline.chrom in primary_ref
        ]

    kmers_names = [bedline.primername for bedline in bedlines]

    # Create mapping array
    # Find index of primary ref
    if mapping == MappingType.FIRST:
        mapping_index = [
            i for i, (k, v) in enumerate(seqdict.items()) if k in primary_ref
        ][0]
        mapping_array, array = create_mapping(array, mapping_index)
    else:
        mapping_array = np.array([x for x in range(array.shape[1])])

    ref_index_to_msa_dict = ref_index_to_msa(mapping_array)

    # Group Primers by basename
    basename_to_line: dict[str, set[BedLine]] = {
        "_".join(name.split("_")[:-1]): set() for name in kmers_names
    }
    for bedline in bedlines:
        basename = "_".join(bedline.primername.split("_")[:-1])
        basename_to_line[basename].add(bedline)

    basename_to_index = {bn: i for i, bn in enumerate(basename_to_line.keys())}

    seq_to_primername = {line.sequence: line.primername for line in bedlines}

    # Create the scoremap
    scoremap = np.empty((array.shape[0], len(basename_to_line)))
    scoremap.fill(None)
    textmap = np.empty((array.shape[0], len(basename_to_line)), dtype="str")
    textmap.fill("None")
    textmap = textmap.tolist()

    # get FPrimer sequences for each basename
    for bn, lines in basename_to_line.items():
        # Get primer size
        primer_len_max = max(len(line.sequence) for line in lines)

        # Set the direction
        if "LEFT" in bn:
            forward = True
            primer_end = list(lines)[0].end
            # Check for the end on a gap edge case and fix it
            if check_for_end_on_gap(ref_index_to_msa_dict, primer_end):
                msa_index = fix_end_on_gap(ref_index_to_msa_dict, primer_end)
            else:
                msa_index = ref_index_to_msa_dict[list(lines)[0].end]
        else:
            forward = False
            msa_index = ref_index_to_msa_dict[list(lines)[0].start]

        # Get the primer sequences
        msa_data = get_primers_from_msa(array, msa_index, forward, primer_len_max)

        # Get the score for each genome
        primer_seqs = {line.sequence for line in lines}

        for genome_index, genome_seq in msa_data.items():
            # Caused by gaps in the msa
            if genome_seq is None:
                if forward:
                    slice = array[genome_index, msa_index - primer_len_max : msa_index]
                    slice[slice == ""] = "-"
                    genome_seq = "".join(slice)
                else:
                    slice = array[genome_index, msa_index : msa_index + primer_len_max]
                    slice[slice == ""] = "-"
                    genome_seq = reverse_complement("".join(slice))

                textmap[genome_index][basename_to_index[bn]] = PlotlyText(
                    primer_seq=[x for x in primer_seqs][0],
                    genome_seq="".join(slice),
                    primer_name=bn,
                ).format_str()
                continue
            # Quick check for exact match
            if genome_seq in primer_seqs:
                scoremap[genome_index, basename_to_index[bn]] = 0
                primer_seq = "".join(primer_seqs.intersection({genome_seq}))
                textmap[genome_index][basename_to_index[bn]] = PlotlyText(
                    primer_seq=primer_seq,
                    genome_seq=genome_seq,
                    primer_name=seq_to_primername.get(primer_seq, "Unknown"),
                ).format_str()
                continue
            # Calculate the hamming distance between all
            seq_to_scores: dict[str, int] = {}
            for primer_seq in primer_seqs:
                seq_to_scores[primer_seq] = calc_primer_hamming(primer_seq, genome_seq)
            scoremap[genome_index, basename_to_index[bn]] = min(
                seq_to_scores.values(),  # type: ignore
            )
            primer_seq = "".join(
                [
                    k
                    for k, v in seq_to_scores.items()
                    if v == scoremap[genome_index, basename_to_index[bn]]
                ][0]
            )
            textmap[genome_index][basename_to_index[bn]] = PlotlyText(
                genome_seq=genome_seq,
                primer_seq=primer_seq,
                primer_name=seq_to_primername.get(primer_seq, "Unknown"),
            ).format_str()

    # Hovertemplate string
    if include_seqs:
        hovertemplatestr = "%{text}<br>" + "<b>Mismatches: %{z}</b><br>"
    else:
        hovertemplatestr = ""

    # Create the heatmap
    fig = go.Figure(
        data=go.Heatmap(
            z=scoremap,
            x=list(basename_to_line.keys()),
            y=[x for x in seqdict.keys()],
            colorscale=plotly_colour_scale,
            text=textmap if include_seqs else None,  # only show text if not minimal
            hovertemplate=hovertemplatestr,
            xgap=0.1,
            ygap=0.1,
            name="Primer Mismatches",
            zmin=0,
            zmax=10,
        )
    )
    chrom_label = (
        f"{scheme_headers.get(list(primary_ref)[0])} ({list(primary_ref)[0]})"
        if scheme_headers.get(list(primary_ref)[0])
        else list(primary_ref)[0]
    )
    fig.update_layout(
        font=dict(family="Roboto, monospace", size=16),
        hoverlabel=dict(font_family="Roboto, monospace"),
        title_text=f"Primer Mismatches: {chrom_label}",
    )
    fig.update_yaxes(autorange="reversed")

    # Remove unnecessary plot elements
    fig.update_layout(
        modebar_remove=[
            "select2d",
            "lasso2d",
            "select",
        ]
    )

    return pio.to_html(
        fig,
        include_plotlyjs=False,
        full_html=False,
        default_height="3000px",
        default_width="100%",
    )


if "${meta.scheme}" != "[]":
    scheme_version_str = "${meta.scheme}"
elif "${meta.custom_scheme_name}" != "[]":
    scheme_version_str = "${meta.custom_scheme_name}"
else:
    scheme_version_str = "Unknown Scheme Name"

payload = {
    "sample_id": "${meta.id}",
    "pipeline_version": "${workflow.manifest.version}",
    "primer_scheme_version": scheme_version_str,
    "minimum_coverage_depth": "${params.min_coverage_depth}",
    "tool_name": "amplicon-nf",
    "tool_version": "${workflow.manifest.version}",
    "citation_link": "https://doi.org/10.5281/zenodo.17522200",
    "contact_email": "s.a.j.wilkinson@bham.ac.uk",
    "funder_statement": "This pipeline has been created as part of the ARTIC network project funded by the Wellcome Trust (collaborator award – 313694/Z/24/Z and discretionary award – 206298/Z/17/Z) and is distributed as open source and open access. All non-code files are made available under a Creative Commons CC-BY licence unless otherwise specified. Please acknowledge or cite this repository or associated publications if used in derived work so we can provide our funders with evidence of impact in the field.",
    "timestamp": datetime.now().strftime("%Y-%m-%d_%H:%M:%S"),
    "qc_table_info": {},
    "single_plots": [],
    "nested_plots": [],
}

samples = set()

primer_scheme = Scheme.from_file("${bed}")

primer_pairs = create_amplicons(primer_scheme.bedlines)

samplesheet_df = pd.read_csv("${samplesheet_csv}", index_col=False)
if "${meta.scheme}" != "[]":
    scheme_samplesheet_df = samplesheet_df[
        samplesheet_df["scheme_name"].str.contains("${meta.scheme}", na=False)
    ]
else:
    scheme_samplesheet_df = samplesheet_df[
        samplesheet_df["custom_scheme_path"].str.contains(
            "${meta.custom_scheme}", na=False
        )
    ]

depth_tsvs = glob("depth_tsvs/*.tsv")
for tsv_path in depth_tsvs:
    sample_name = tsv_path.split("/")[-1].split(".")[0]
    payload["qc_table_info"].setdefault(sample_name, {})
    df = pd.read_csv(
        tsv_path, sep="\t", index_col=None, names=["chrom", "pos", "depth"]
    )
    payload["qc_table_info"][sample_name]["mean_depth"] = (
        round(df["depth"].mean(), 2) if len(df) > 0 else 0.0
    )
    payload["qc_table_info"][sample_name]["primer_scheme"] = scheme_version_str
    payload["qc_table_info"][sample_name]["coverage"] = (
        round(
            len(df[df["depth"] >= int("${params.min_coverage_depth}")]) / len(df) * 100,
            2,
        )
        if len(df) > 0
        else 0.0
    )
    if payload["qc_table_info"][sample_name]["coverage"] >= int(
        "${params.qc_pass_high_coverage}"
    ):
        payload["qc_table_info"][sample_name]["qc_result"] = "pass"
    elif payload["qc_table_info"][sample_name]["coverage"] >= int(
        "${params.qc_pass_min_coverage}"
    ):
        payload["qc_table_info"][sample_name]["qc_result"] = "warning"
    else:
        payload["qc_table_info"][sample_name]["qc_result"] = "fail"

coverage_tsvs = glob("coverage_tsvs/*.txt")
for tsv_path in coverage_tsvs:
    sample_name = tsv_path.split("/")[-1].split(".")[0]
    samples.add(sample_name)
    with open(tsv_path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\\t")
        reads = {x["#rname"]: x["numreads"] for x in reader}
        payload["qc_table_info"][sample_name]["total_reads"] = sum(
            int(x) for x in reads.values()
        )

amplicon_depth_rows = []
amp_depth_tsvs = glob("amplicon_depth_tsvs/*.tsv")
for tsv_path in amp_depth_tsvs:
    sample_name = tsv_path.split("/")[-1].split(".")[0]
    samples.add(sample_name)
    with open(tsv_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = []
        for row in reader:
            row_dict = {"sample": sample_name, **row}
            row_dict["amplicon"] = int(row_dict["amplicon"])
            # row_dict["amplicon"] = f"Amplicon {row_dict['amplicon']}"
            rows.append(row_dict)

        amplicon_depth_rows.extend(rows)

        for x in primer_pairs:
            if (
                len(
                    [
                        y
                        for y in amplicon_depth_rows
                        if y["chrom"] == str(x.chrom)
                        and y["amplicon"] == x.amplicon_number
                    ]
                )
                == 0
            ):
                amplicon_depth_rows.append(
                    {
                        "sample": sample_name,
                        "chrom": str(x.chrom),
                        "amplicon": x.amplicon_number,
                        "mean_depth": 0.0,
                    }
                )

    payload["qc_table_info"][sample_name]["total_amp_dropouts"] = (
        len(
            [
                x
                for x in rows
                if float(x["mean_depth"]) < int("${params.min_coverage_depth}")
            ]
        )
        if len(rows) > 0
        else len(primer_pairs)
    )

for row in scheme_samplesheet_df.itertuples():
    if not payload["qc_table_info"].get(row.sample):
        samples.add(row.sample)
        payload["qc_table_info"].setdefault(row.sample, {})
        payload["qc_table_info"][row.sample]["primer_scheme"] = scheme_version_str
        payload["qc_table_info"][row.sample]["coverage"] = 0.0
        payload["qc_table_info"][row.sample]["mean_depth"] = 0.0
        payload["qc_table_info"][row.sample]["total_reads"] = 0
        payload["qc_table_info"][row.sample]["total_amp_dropouts"] = len(primer_pairs)
        payload["qc_table_info"][row.sample]["qc_result"] = "fail"

    if len([x for x in amplicon_depth_rows if x["sample"] == row.sample]) == 0:
        for x in primer_pairs:
            amplicon_depth_rows.append(
                {
                    "sample": row.sample,
                    "chrom": str(x.chrom),
                    "amplicon": x.amplicon_number,
                    "mean_depth": 0.0,
                }
            )

amplicon_depth_rows = sorted(
    amplicon_depth_rows,
    key=lambda x: (x["sample"], x["chrom"], x["amplicon"]),
)

# Sort the qc table info by sample name
payload["qc_table_info"] = dict(
    sorted(
        payload["qc_table_info"].items(),
        key=lambda item: item[0],
    )
)
with open(f"{scheme_version_str.replace("/", "_")}_qc_results.tsv", "w") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "sample",
            "primer_scheme",
            "coverage",
            "mean_depth",
            "total_reads",
            "total_amp_dropouts",
            "qc_result",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    for sample, row in payload["qc_table_info"].items():
        row["sample"] = sample
        writer.writerow(row)

# amplicon_depth_rows.sort(key=lambda x: int(x["amplicon"].replace("Amplicon ", "")))
amplicon_depth_df = pd.DataFrame(amplicon_depth_rows)

chroms = amplicon_depth_df.chrom.unique()
amp_depth_heatmaps = {"name": "Amplicon Depths", "plots": []}
for chrom in chroms:
    chrom_df = amplicon_depth_df.loc[amplicon_depth_df["chrom"] == chrom]
    chrom_label = (
        f"{primer_scheme.header_dict.get(chrom)} ({chrom})"
        if primer_scheme.header_dict.get(chrom)
        else chrom
    )
    amp_depth_heatmaps["plots"].append(
        {
            "name": chrom,
            "plot_html": amplicon_depth_heatmap(
                amplicon_depths=chrom_df,
                scheme_str=scheme_version_str,
                chrom_name=chrom_label,
            ),
        }
    )
payload["nested_plots"].append(amp_depth_heatmaps)

msa_list = glob("msas/*.fa*")
if len(msa_list) > 0:
    primer_mismatch_heatmaps = {"name": "Primer Mismatches", "plots": []}
    for msa_path in msa_list:
        msa, seqdict = parse_msa(msa_path)
        contig_name = msa_path.split("/")[-1].split("_")[0]

        primer_mismatch_heatmaps["plots"].append(
            {
                "name": contig_name,
                "plot_html": primer_mismatch_heatmap(
                    array=msa, seqdict=seqdict, bedfile="${bed}"
                ),
            }
        )
    payload["nested_plots"].append(primer_mismatch_heatmaps)

render_qc_report(
    payload=payload,
    template_path=Path("${report_template}"),
    output_path=Path(
        f"{scheme_version_str.replace("/", "_")}_amplicon-nf_run-report.html"
    ),
    bootstrap_css_path=Path("${bootstrap_bundle_min_css}"),
    bootstrap_bundle_js_path=Path("${bootstrap_bundle_min_js}"),
    plotly_js_path=Path("${plotly_js}"),
    svg_logo_path=Path("${artic_logo_svg}"),
)

with open("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    for package in ("plotly", "primalbedtools", "pandas", "jinja2"):
        f.write(f"  {package}: {version(package)}\\n")
