#!/usr/bin/env python3
#
# Per-sample half of the Pathoplexus submission template (on unless --skip_pathoplexus_template).
# Writes <id>.pathoplexus.json (one metadata row + bookkeeping) and <id>.pathoplexus.fasta
# (consensus records renamed to their Pathoplexus FASTA IDs) and writes raw_reads/: copies of
# the dehosted read files referenced in files.rawReads, with the FASTQ headers stripped to the
# read ID (run, flow cell, barcode, start time and basecaller details removed).
#
# NOTE: this file is a Nextflow template. Every literal dollar sign and backslash would be
# interpreted by Nextflow, so the script deliberately uses neither outside the block of
# Nextflow values below (tab/newline come from chr(), regexes avoid backslash classes).

import csv
import gzip
import json
import os
import platform
import re
import statistics
import sys
from datetime import datetime
from pathlib import Path

# ---- values injected by Nextflow ----------------------------------------------------------
SAMPLE_ID = "${meta.id}"
PLATFORM = "${meta.platform}"
SCHEME = "${meta.scheme ?: ''}"
CUSTOM_SCHEME_NAME = "${meta.custom_scheme_name ?: ''}"
PIPELINE_NAME = "${workflow.manifest.name}"
PIPELINE_VERSION = "${workflow.manifest.version}"
MIN_DEPTH = int("${params.min_coverage_depth}")
QC_PASS_HIGH = float("${params.qc_pass_high_coverage}")
QC_PASS_MIN = float("${params.qc_pass_min_coverage}")
SKIP_DEHOSTING = "${params.skip_dehosting}" == "true"
DEHOSTING_INDEX = "${params.dehosting_index}"
MIN_ONT_READ_QUALITY = "${params.min_ont_read_quality ?: ''}"
MIN_ONT_READ_LENGTH = "${params.min_ont_read_length ?: ''}"
MAX_ONT_READ_LENGTH = "${params.max_ont_read_length ?: ''}"
SKIP_ONT_QUALITY_CHECK = "${params.skip_ont_quality_check}" == "true"
MULTILINE_FASTA = "${params.multiline_fasta}" == "true"
CONSENSUS = "${consensus}"
BED = "${bed}"
QC_TSV = "${qc_tsv}"
TASK_PROCESS = "${task.process}"
# -------------------------------------------------------------------------------------------

TAB = chr(9)
NL = chr(10)
NL_BYTES = NL.encode()
READS_DIR = Path("reads")  # staged reads (read from, never published)
UPLOAD_READS_DIR = Path("raw_reads")  # header-stripped copies for upload
MAX_RAW_READ_BYTES = 2 * 1024**3  # Pathoplexus per-file upload limit
FILENAME_OK = "[A-Za-z0-9_.-]+"  # Pathoplexus raw read file name rule
INSDC_ACCESSION = "(?:[A-Z]{1,2}[0-9]{5,8}|[A-Z]{2}_[0-9]{6,9})[.][0-9]+"
ILLUMINA_HEADERS_TO_SCAN = 1000

# ONT flow cell ID -> flow cell type. Only the Flongle format is documented by ONT (ABC123);
# the F... (MinION/GridION) and P... (PromethION) prefixes are a widely observed convention.
ONT_FLOW_CELL_TYPES = [
    ("[A-Z]{3}[0-9]{3}", "Flongle"),
    ("F[A-Z]{2}[0-9]{5}", "MinION/GridION"),
    ("P[A-Z]{2}[0-9]{5}", "PromethION"),
]

# Illumina instrument ID (first field of the read name) -> Pathoplexus sequencingInstrument.
# Community-maintained prefix table (fastp PR #508, Biostars), refined with Illumina's
# DRAGEN TSO 500 input-check table where it is more specific (NextSeq 500 vs 550).
# Prefixes the community table leaves ambiguous stay at platform level (ILLUMINA):
# VH/VL (NextSeq 1000 or 2000), LH (NovaSeq X or X Plus), SH (MiSeq i100, no vocabulary option).
ILLUMINA_INSTRUMENT_PREFIXES = [
    ("MN[0-9]", "Illumina MiniSeq"),
    ("M[0-9]", "Illumina MiSeq"),
    ("FS[0-9]", "Illumina iSeq 100"),
    ("N[SBL]50", "NextSeq 500"),
    ("N[SB]55", "NextSeq 550"),
    ("NA[0-9]", "Illumina NovaSeq 6000"),
    ("A[0-9]", "Illumina NovaSeq 6000"),
    ("LL[0-9]", "Illumina NovaSeq X"),
]

QC_DETERMINATION = {
    "pass": "sequence passed quality control",
    "warning": "sequence flagged for potential quality control issues",
    "fail": "sequence failed quality control",
}
QC_SEVERITY = {"pass": 0, "warning": 1, "fail": 2}

warnings = []


def warn(message):
    warnings.append(message)
    print(f"WARNING [pathoplexus] {SAMPLE_ID}: {message}", file=sys.stderr)


def is_set(value):
    return value not in ("", "[]", "null")


def open_text(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == bytes([0x1F, 0x8B]):
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")


def read_files():
    """Staged read files: a single file, an Illumina pair, or (dehosting skipped) an ONT directory."""
    files = []
    if not READS_DIR.is_dir():
        return files
    for entry in sorted(READS_DIR.iterdir()):
        if entry.is_dir():
            # Mirror guppyplex/deacon: top-level *.fastq* only, skipping *.temp
            files += [
                f
                for f in sorted(entry.iterdir())
                if f.is_file() and ".fastq" in f.name and not f.name.endswith(".temp")
            ]
        elif entry.is_file():
            files.append(entry)
    return files


def scan_ont_headers(files):
    """Scan every read header: read order is not chronological (multiple runs, deacon, guppyplex)."""
    earliest, earliest_raw = None, None
    flow_cells, models = set(), set()
    run_flow_cells, run_starts = {}, {}  # per MinKNOW run: its flow cell(s) and earliest read start
    n_reads = 0
    for path in files:
        with open_text(path) as fh:
            for i, line in enumerate(fh):
                if i % 4 != 0 or not line.startswith("@"):
                    continue
                n_reads += 1
                fields = {}
                for token in line.split()[1:]:
                    key, sep, value = token.partition("=")
                    if sep:
                        fields[key] = value
                run = fields.get("runid", "")
                if run:
                    run_flow_cells.setdefault(run, set())
                if fields.get("flow_cell_id"):
                    flow_cells.add(fields["flow_cell_id"])
                    if run:
                        run_flow_cells[run].add(fields["flow_cell_id"])
                if fields.get("basecall_model_version_id"):
                    models.add(fields["basecall_model_version_id"])
                if fields.get("start_time"):
                    try:
                        stamp = datetime.fromisoformat(fields["start_time"])
                    except ValueError:
                        continue
                    try:
                        if earliest is None or stamp < earliest:
                            earliest, earliest_raw = stamp, fields["start_time"]
                        if run and (run not in run_starts or stamp < run_starts[run]):
                            run_starts[run] = stamp
                    except TypeError:  # naive vs aware timestamps, keep the first seen
                        pass
    ordered_runs = sorted(
        run_flow_cells, key=lambda r: (r not in run_starts, run_starts[r].timestamp() if r in run_starts else 0, r)
    )
    return {
        "n_reads": n_reads,
        # The date as recorded by the sequencing device (its own UTC offset)
        "sequencing_date": earliest_raw[:10] if earliest_raw else "",
        "flow_cell_ids": sorted(flow_cells),
        "basecall_models": sorted(models),
        # In run start order; run IDs shortened to 8 characters as in MinKNOW file names
        "runs": [{"run": r[:8], "flow_cells": sorted(run_flow_cells[r])} for r in ordered_runs],
    }


def scan_illumina_headers(files):
    instruments = set()
    if files:
        with open_text(files[0]) as fh:  # R1 sorts first
            for i, line in enumerate(fh):
                if i >= 4 * ILLUMINA_HEADERS_TO_SCAN:
                    break
                if i % 4 == 0 and line.startswith("@"):
                    name = line[1:].split()[0] if line[1:].split() else ""
                    if ":" in name:
                        instruments.add(name.split(":")[0])
    return {"instrument_ids": sorted(instruments)}


def flow_cell_type(flow_cell_id):
    for pattern, fc_type in ONT_FLOW_CELL_TYPES:
        if re.fullmatch(pattern, flow_cell_id):
            return fc_type
    return None


def flow_cell_label(flow_cell_ids):
    # "MinION/GridION flow cell FBH86191"; the type is only stated when all IDs share one known type
    types = sorted({flow_cell_type(x) or "unknown" for x in flow_cell_ids})
    type_text = f"{types[0]} " if len(types) == 1 and types[0] != "unknown" else ""
    return f"{type_text}flow cell{'s' if len(flow_cell_ids) > 1 else ''} {', '.join(flow_cell_ids)}"


def illumina_model(instrument_id):
    for pattern, model in ILLUMINA_INSTRUMENT_PREFIXES:
        if re.match(pattern, instrument_id):
            return model
    return None


def read_fasta(path):
    records, header, chunks = [], None, []
    with open(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header, chunks = line[1:], []
            elif line:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def chrom_from_header(header):
    # Consensus headers are ">{id} {chrom}_artic-network/amplicon-nf_{version}"
    parts = header.split(None, 1)
    if len(parts) < 2:
        return parts[0] if parts else ""
    description = parts[1]
    if "_artic-network/" in description:
        return description.rsplit("_artic-network/", 1)[0]
    return description.split()[0]


def safe_id(text):
    return re.sub("[^A-Za-z0-9_.-]", "_", text)


def wrap(seq, width=60):
    if not MULTILINE_FASTA:
        return seq
    return NL.join(seq[i : i + width] for i in range(0, len(seq), width))


def amplicon_size():
    parts = SCHEME.split("/") if is_set(SCHEME) else []
    if len(parts) == 3 and parts[1].isdigit():
        return f"{parts[1]}bp"
    # Custom scheme: median amplicon length (LEFT start -> RIGHT end) rounded to the nearest 100
    starts, ends = {}, {}
    with open(BED, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) < 4:
                continue
            match = re.search("_([0-9]+)_(LEFT|RIGHT)", cols[3])
            if not match:
                continue
            key = (cols[0], int(match.group(1)))
            if match.group(2) == "LEFT":
                starts[key] = min(starts.get(key, int(cols[1])), int(cols[1]))
            else:
                ends[key] = max(ends.get(key, int(cols[2])), int(cols[2]))
    lengths = [ends[k] - starts[k] for k in starts if k in ends and ends[k] > starts[k]]
    if not lengths:
        warn("could not determine amplicon size from the primer BED")
        return ""
    return f"{max(100, int(statistics.median(lengths) / 100 + 0.5) * 100)}bp"


def scheme_name():
    # The samplesheet scheme name. Reference-selecting schemes (e.g. artic-pan-ebola) use the same primers
    # for every sample, artic only remaps reads to the closest variant reference, so no variant suffix.
    if is_set(SCHEME):
        return SCHEME
    return CUSTOM_SCHEME_NAME if is_set(CUSTOM_SCHEME_NAME) else ""


def strip_fastq_headers(src, dest):
    """Copy a gzipped FASTQ, reducing each header to "@<read ID>" and each separator line to "+"."""
    # Works on blocks of whole records (much faster than line by line). Compression dominates the
    # run time: level 3 is ~3x faster than gzip's default 6 for ~8% larger files. mtime=0 keeps
    # the output byte-identical between runs.
    carry = []
    with gzip.open(src, "rb") as fin, gzip.GzipFile(dest, "wb", compresslevel=3, mtime=0) as fout:
        while True:
            block = fin.readlines(16 * 1024**2)
            if not block:
                if carry:
                    raise ValueError(f"{src.name}: truncated FASTQ record at the end of the file")
                break
            lines = carry + block
            n = len(lines) - len(lines) % 4
            records, carry = lines[:n], lines[n:]
            headers, separators = records[0::4], records[2::4]
            if not all(h.startswith(b"@") for h in headers) or not all(s.startswith(b"+") for s in separators):
                raise ValueError(f"{src.name}: not a 4-line FASTQ file")
            records[0::4] = [(h.split(None, 1) or [b"@"])[0] + NL_BYTES for h in headers]
            records[2::4] = [b"+" + NL_BYTES] * len(separators)
            fout.write(b"".join(records))


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


# ---- consensus -----------------------------------------------------------------------------
records = []
for header, seq in read_fasta(CONSENSUS):
    seq = seq.upper()
    non_n = len(seq) - seq.count("N") - seq.count("-")
    records.append({"chrom": chrom_from_header(header), "seq": seq, "non_n": non_n})

kept = [r for r in records if r["non_n"] > 0]
for r in records:
    if r["non_n"] == 0:
        warn(f"consensus record for {r['chrom']} is entirely N, it is left out of the Pathoplexus files")

excluded = not kept
exclusion_reason = ""
if excluded:
    exclusion_reason = "no consensus sequence with called bases"
    warn("excluded from the Pathoplexus template (consensus is entirely N), its reads are not prepared either")

if len(kept) == 1:
    fasta_ids = [SAMPLE_ID]
else:
    fasta_ids = [f"{SAMPLE_ID}_{safe_id(r['chrom'])}" for r in kept]
    if len(set(fasta_ids)) != len(fasta_ids):
        fasta_ids = [f"{SAMPLE_ID}_{i + 1}" for i in range(len(kept))]

# ---- run information -----------------------------------------------------------------------
reads = read_files()
if PLATFORM == "nanopore":
    header_info = scan_ont_headers(reads)
else:
    header_info = scan_illumina_headers(reads)

qc_rows = {}
with open(QC_TSV, "rt") as fh:
    for row in csv.DictReader(fh, delimiter=TAB):
        qc_rows[row["contig"]] = row

# ---- metadata row --------------------------------------------------------------------------
pipeline_label = f"{PIPELINE_NAME} v{PIPELINE_VERSION}"
scheme_used = scheme_name()
size = amplicon_size()

if PLATFORM == "nanopore":
    fc_ids = header_info["flow_cell_ids"]
    fc_types = sorted({flow_cell_type(x) or "unknown" for x in fc_ids})
    instrument = "PromethION" if fc_ids and fc_types == ["PromethION"] else "OXFORD_NANOPORE"
    chemistries = sorted(
        {m.group(1).upper() for m in (re.search("_(r[0-9.]+)_", x) for x in header_info["basecall_models"]) if m}
    )
    platform_text = "Oxford Nanopore sequencing"
    if chemistries:
        platform_text += f" ({', '.join(chemistries)} chemistry)"
    runs = header_info["runs"]
    if len(runs) > 1:
        # Each run with its own flow cell(s), e.g. "run 993a8f43 on MinION/GridION flow cell FBH86191"
        per_run = [f"run {r['run']}" + (f" on {flow_cell_label(r['flow_cells'])}" if r["flow_cells"] else "") for r in runs]
        platform_text += f", {len(runs)} sequencing runs ({', '.join(per_run)})"
    elif fc_ids:
        platform_text += f" on {flow_cell_label(fc_ids)}"

    if SKIP_ONT_QUALITY_CHECK:
        filters = ["quality check skipped"]
    else:
        filters = [f"min Q{MIN_ONT_READ_QUALITY}"] if is_set(MIN_ONT_READ_QUALITY) else []
    if is_set(MIN_ONT_READ_LENGTH):
        filters.append(f"min length {MIN_ONT_READ_LENGTH}")
    if is_set(MAX_ONT_READ_LENGTH):
        filters.append(f"max length {MAX_ONT_READ_LENGTH}")
    processing = []
    if header_info["basecall_models"]:
        processing.append(f"basecalled with {', '.join(header_info['basecall_models'])}")
    processing.append(f"reads filtered with artic guppyplex ({', '.join(filters) or 'default settings'}) in {pipeline_label}")
    sequencing_date = header_info["sequencing_date"]
else:
    ids = header_info["instrument_ids"]
    models = {illumina_model(x) for x in ids}
    instrument = models.pop() if len(models) == 1 and None not in models else "ILLUMINA"
    platform_text = "Illumina paired-end sequencing" + (f" on {instrument}" if instrument != "ILLUMINA" else "")
    processing = [f"adapter and quality trimming with fastp in {pipeline_label}"]
    sequencing_date = ""

# Wet lab / sequencing only: data processing is described by the processing, dehosting and consensus fields
protocol = ["Tiled amplicon sequencing" + (f" with primer scheme {scheme_used}" if scheme_used else "") + (f" ({size} amplicons)" if size else "")]
protocol.append(platform_text)

# QC: depth/breadth are pre-normalisation (normalisation caps depth); the pass/warning/fail call
# uses the same measure as the QC report (post-normalisation when normalisation ran).
kept_qc = [(r, qc_rows.get(r["chrom"])) for r in kept]
missing_qc = [r["chrom"] for r, q in kept_qc if q is None]
if missing_qc and not excluded:
    warn(f"no QC row for {', '.join(missing_qc)}")
kept_qc = [(r, q) for r, q in kept_qc if q is not None]


def weighted(column):
    pairs = [(to_float(q.get(column)), len(r["seq"])) for r, q in kept_qc]
    pairs = [(v, w) for v, w in pairs if v is not None]
    total = sum(w for _, w in pairs)
    return round(sum(v * w for v, w in pairs) / total) if total else ""


qc_result = ""
issues, details = [], []
if kept_qc:
    qc_result = max((q.get("qc_result", "") for _, q in kept_qc), key=lambda x: QC_SEVERITY.get(x, -1))
    low_coverage, dropout = False, False
    for r, q in kept_qc:
        post = to_float(q.get("post_norm_coverage"))
        measure = post if post is not None else to_float(q.get("pre_norm_coverage")) or 0.0
        label = "post-normalisation" if post is not None else "pre-normalisation"
        dropouts = [x for x in (q.get("amplicon_dropouts") or "").split(",") if x]
        low_coverage = low_coverage or measure < QC_PASS_HIGH
        dropout = dropout or bool(dropouts) or (to_float(q.get("total_amp_dropouts")) or 0) > 0
        prefix = f"{r['chrom']}: " if len(kept_qc) > 1 else ""
        details.append(
            f"{prefix}{measure:.2f}% of reference at >={MIN_DEPTH}x ({label} depth; pass >={QC_PASS_HIGH:g}%, "
            f"warning >={QC_PASS_MIN:g}%); mean depth {q.get('pre_norm_mean_depth', '')}x; "
            f"amplicon dropouts (<{MIN_DEPTH}x): {', '.join(dropouts) if dropouts else 'none'}"
        )
    if qc_result in ("warning", "fail"):
        if low_coverage:
            issues.append("low genome coverage")
        if dropout:
            issues.append("amplicon dropout")

accessions = [r["chrom"] for r in kept if re.fullmatch(INSDC_ACCESSION, r["chrom"])]

dehosting = ""
if not SKIP_DEHOSTING:
    index_name = os.path.basename(DEHOSTING_INDEX.rstrip("/"))
    if index_name.endswith(".idx"):
        index_name = index_name[: -len(".idx")]
    dehosting = f"Deacon (index {index_name}) in {pipeline_label}"

# AUTO: every field this module fills; anything else in the template config is left for the user
AUTO = {
    "id": SAMPLE_ID,
    "sequencingInstrument": instrument,
    "sequencingAssayType": "AMPLICON",
    "sequencingDate": sequencing_date,
    "sequencingProtocol": "; ".join(protocol) + ".",
    "ampliconPcrPrimerScheme": scheme_used,
    "ampliconSize": size,
    "assemblyReferenceGenomeAccession": "; ".join(accessions),
    "rawSequenceDataProcessingMethod": "; ".join(processing),
    "dehostingMethod": dehosting,
    "consensusSequenceSoftwareName": PIPELINE_NAME,
    "consensusSequenceSoftwareVersion": PIPELINE_VERSION,
    "depthOfCoverage": weighted("pre_norm_mean_depth"),
    "breadthOfCoverage": weighted("pre_norm_coverage"),
    "qualityControlMethodName": PIPELINE_NAME,
    "qualityControlMethodVersion": PIPELINE_VERSION,
    "qualityControlDetermination": QC_DETERMINATION.get(qc_result, ""),
    "qualityControlIssues": "; ".join(issues),
    "qualityControlDetails": " | ".join(details),
    "files.rawReads": "",
}

# ---- raw reads -----------------------------------------------------------------------------
# Only the dehosted read files of included samples are written, so raw_reads/ holds exactly the
# files referenced in files.rawReads (Pathoplexus blocks folder uploads with unreferenced files).
# Headers are stripped to the read ID: the rest of a header can carry run, flow cell, barcode,
# sample and timing details that are not needed for submission.
raw_read_files = []
if READS_DIR.is_dir() and not SKIP_DEHOSTING and not excluded:
    candidates = [e for e in sorted(READS_DIR.iterdir()) if e.is_file() and e.name.endswith(".dehosted.fastq.gz")]
    bad = [e.name for e in candidates if not re.fullmatch(FILENAME_OK, e.name)]
    if bad:
        warn(f"read file names {', '.join(bad)} contain characters Pathoplexus does not accept, reads not prepared")
    elif candidates:
        UPLOAD_READS_DIR.mkdir()
        for e in candidates:
            dest = UPLOAD_READS_DIR / e.name
            strip_fastq_headers(e, dest)
            if dest.stat().st_size > MAX_RAW_READ_BYTES:
                warn(f"{e.name} is larger than the 2 GiB Pathoplexus upload limit")
            raw_read_files.append(e.name)
AUTO["files.rawReads"] = " ".join(raw_read_files)

# ---- outputs -------------------------------------------------------------------------------
with open(f"{SAMPLE_ID}.pathoplexus.fasta", "w") as fh:
    if not excluded:
        for fasta_id, r in zip(fasta_ids, kept):
            fh.write(f">{fasta_id}{NL}{wrap(r['seq'])}{NL}")

with open(f"{SAMPLE_ID}.pathoplexus.json", "w") as fh:
    json.dump(
        {
            "id": SAMPLE_ID,
            "excluded": excluded,
            "exclusion_reason": exclusion_reason,
            "fasta_ids": [] if excluded else fasta_ids,
            "auto_fields": sorted(AUTO),
            "row": {} if excluded else AUTO,
            "raw_read_files": raw_read_files,
            "warnings": warnings,
        },
        fh,
        indent=2,
    )

with open("versions.yml", "w") as fh:
    fh.write(f'"{TASK_PROCESS}":{NL}  python: {platform.python_version()}{NL}')
