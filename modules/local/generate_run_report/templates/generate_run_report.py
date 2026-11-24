#!/usr/bin/env python3

import csv
import pathlib
from collections import Counter
from datetime import datetime
from enum import Enum
from glob import glob
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from Bio import SeqIO
from jinja2 import Environment, FileSystemLoader, select_autoescape
from primalbedtools.amplicons import create_amplicons
from primalbedtools.bedfiles import BedLine
from primalbedtools.scheme import Scheme

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
        "mean_fragment_depth"
    ].fillna(0) # type: ignore

    # amplicon_depths.reindex(columns=columns)

    fig = px.imshow(
        amplicon_depths,
        x=amplicon_depths.columns,
        y=amplicon_depths.index,
        color_continuous_scale=plotly_colour_scale,
        labels=dict(x="Amplicon", y="Sample", color="Mean Fragment Depth"),
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

    primer_scheme = Scheme.from_file(bedfile) # type: ignore

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


# align_trim (>= 2.0) per-sample QC reports, written with --qc-prefix
ALIGNTRIM_REPORT_KINDS = (
    "summary",
    "amplicons",
    "primer_pairs",
    "split_reads",
    "overlaps",
)
PRIMER_PAIR_TYPES = [
    "same_amplicon",
    "at_overlap_primers",
    "inside_overlap",
    "spans_amplicons_same_pool",
    "spans_amplicons_other_pool",
    "different_references",
    "other_combination",
]
# Recorded in addition to a fragment's own pair type, so never summed with the above
PRIMER_PAIR_EXTRA_TYPES = ["left_primer_at_both_ends", "right_primer_at_both_ends"]
SPLIT_READ_TYPES = ["amplicon_chimera", "joined_inside_amplicon", "low_mapq_part"]
def param_text(name: str, value: str, unit: str = "") -> str:
    """Parameter name as HTML code, with this run's value when it is a positive integer.
    Definitions are trusted constants rendered unescaped in the report keys."""
    code = f"<code>{name}</code>"
    return f"{code} ({value}{unit} in this run)" if value.isdigit() and int(value) > 0 else code


PRIMER_THRESHOLD = param_text("--primer_match_threshold", "${params.primer_match_threshold}", " bp")
MIN_MAPQ = param_text("--min_mapping_quality", "${params.min_mapping_quality}")
NORMALISE_DEPTH = param_text("--normalise_depth", "${params.normalise_depth}")

# What each definition is judged on, shown as a tag in the report keys
BASIS_LABELS = {
    "aligner": "aligner flags / MAPQ",
    "position": "primer position",
    "sequence": "primer sequence",
    "depth": "depth",
    "bases": "base frequencies",
}

# Read fates in align_trim pipeline order: fixed colour (stable across runs and schemes,
# not encoding good / bad), what the fate is judged on, and the definition shown in the report key
READ_FATES = {
    "kept": (
        "#345c67",
        None,
        "Passed every check below: assigned to an amplicon, primer-trimmed and written to the BAM used for variant calling.",
    ),
    "unmapped": (
        "#999999",
        "aligner",
        "The read (or either mate of a pair) has no alignment to the reference.",
    ),
    "not_proper_pair": (
        "#bab0ac",
        "aligner",
        "Illumina only: both mates aligned, but the aligner did not set the proper-pair flag.",
    ),
    "missing_mate": (
        "#d4a6c8",
        "aligner",
        "Illumina only: only one mate of the pair is present in the alignment.",
    ),
    "low_mapq": (
        "#9d7660",
        "aligner",
        f"Mapping quality of the alignment (of either mate) is below {MIN_MAPQ}.",
    ),
    "no_alignment": (
        "#79706e",
        "aligner",
        "Flagged as mapped, but with no aligned reference span.",
    ),
    "no_primer_match": (
        "#e69f00",
        "position",
        f"An aligned end of the fragment lies outside every amplicon (more than {PRIMER_THRESHOLD} beyond the outer primers), "
        "so no primer could be assigned to it. The read sequence is not checked.",
    ),
    "ambiguous_amplicon": (
        "#56b4e9",
        "position",
        "The aligned ends were assigned LEFT n+1 and RIGHT n: the fragment lies entirely within the overlap of two neighbouring "
        "amplicons, so it cannot be placed in either.",
    ),
    "mismatched_primers": (
        "#cc79a7",
        "position",
        "The primers assigned to the two aligned ends belong to different amplicons (other than the overlap case above).",
    ),
    "trim_failed": (
        "#0072b2",
        "position",
        "Soft-masking the assigned primers' reference coordinates from the alignment failed.",
    ),
    "primer_only": (
        "#f0e442",
        "position",
        "No aligned bases remained once the assigned primers' reference coordinates were soft-masked: the alignment lies "
        "entirely within primer sites.",
    ),
    "not_full_length": (
        "#009e73",
        "position",
        "The alignment does not reach both primer sites of its amplicon (only applied when full-length fragments are required).",
    ),
    "by_normalise": (
        "#d3b88c",
        "depth",
        f"Passed every other check, but not kept because its amplicon had already reached the fragment depth target, {NORMALISE_DEPTH}.",
    ),
}
# Fates added by a newer align_trim get these, assigned in sorted name order
UNKNOWN_FATE_COLOURS = ["#4d4d4d", "#737373", "#a6a6a6", "#cccccc"]
UNKNOWN_DEFINITION = "Reported by align_trim; no description available in this pipeline version."

# Definitions for the keys under the other Primer Trimming QC plots / tables, as (basis, text)
# Amplicons key: the Read Fates definition plus which amplicon the fragment is counted against
AMPLICON_METRIC_DEFINITIONS = {
    "mismatched_primers": (
        "position",
        READ_FATES["mismatched_primers"][2]
        + " Counted against the amplicon of the primer assigned to the left end.",
    ),
    "no_primer_match": (
        "position",
        READ_FATES["no_primer_match"][2]
        + " Counted against the amplicon of the primer assigned to the other end; fragments with neither end assigned "
        "are only counted in Read Fates.",
    ),
}
# Metric shown first (the default) in the Amplicons heatmap dropdown
AMPLICON_HEATMAP_DEFAULT = "removed_mismatched_primers"
PRIMER_PAIR_DEFINITIONS = {
    "same_amplicon": (
        "position",
        "Both aligned ends were assigned primers of the same amplicon: the intended product.",
    ),
    "at_overlap_primers": (
        "position",
        f"Ends assigned LEFT n+1 and RIGHT n, and both aligned ends are within {PRIMER_THRESHOLD} of those primers' 5' ends: "
        "a short product spanning just the overlap between neighbouring amplicons.",
    ),
    "inside_overlap": (
        "position",
        f"Ends assigned LEFT n+1 and RIGHT n, but the aligned ends are not within {PRIMER_THRESHOLD} of those primers' 5' ends: "
        "part of a normal amplicon that lies within the overlap.",
    ),
    "spans_amplicons_same_pool": (
        "position",
        "Ends assigned LEFT n and RIGHT n+1 or further, both from the same pool: the fragment spans the amplicon(s) in between.",
    ),
    "spans_amplicons_other_pool": (
        "position",
        "Ends assigned LEFT n and RIGHT n+1 or further, from different pools: the fragment spans the amplicon(s) in between.",
    ),
    "different_references": (
        "position",
        "The two ends were assigned primers on different reference sequences.",
    ),
    "other_combination": (
        "position",
        "Any other combination of assigned primers, e.g. LEFT n with RIGHT n-2 or further back.",
    ),
    "left_primer_at_both_ends": (
        "sequence",
        "The right end of the read carries a left primer's sequence, and the left end is at a left primer (by position, or by the "
        "same sequence search). Found by searching both read ends, including soft-clipped bases, for the last 12 bases of every "
        "primer. Recorded in addition to the fragment's own type; it does not change what happens to the fragment.",
    ),
    "right_primer_at_both_ends": (
        "sequence",
        "The left end of the read carries a right primer's sequence, and the right end is at a right primer (by position, or by "
        "the same sequence search). Found by searching both read ends, including soft-clipped bases, for the last 12 bases of "
        "every primer. Recorded in addition to the fragment's own type; it does not change what happens to the fragment.",
    ),
}
SPLIT_READ_DEFINITIONS = {
    "amplicon_chimera": (
        "position",
        f"Both breakpoints of the join are within {PRIMER_THRESHOLD} of a primer 5' end (an amplicon end): "
        "two amplicon products joined end to end.",
    ),
    "joined_inside_amplicon": (
        "position",
        f"At least one breakpoint of the join is more than {PRIMER_THRESHOLD} from every primer 5' end.",
    ),
    "low_mapq_part": (
        "aligner",
        f"One of the joined alignment parts has mapping quality below {MIN_MAPQ}; its breakpoint positions are not classified.",
    ),
}
OVERLAP_METRIC_DEFINITIONS = {
    "max_freq_difference": (
        "bases",
        "Base frequencies (A, C, G, T, deletion and insertion) are counted separately for each amplicon from its kept, "
        "correctly paired fragments, after primer trimming and normalisation. At every position in the overlap (where both "
        "amplicon inserts cover the reference, primers excluded) with a depth of at least 10 fragments in both amplicons, "
        "the largest difference in any base's frequency is taken; this is the maximum over those positions. 0 means "
        "identical base frequencies, 1 completely different. Blank (NA) where no position reached that depth in both amplicons.",
    ),
}


# amplicons.tsv columns not offered in the Amplicons heatmap dropdown (still in the TSV)
AMPLICON_HEATMAP_EXCLUDED = {
    "removed_trim_failed",
    "removed_not_full_length",
    "removed_by_normalise",
}


def key_items(names: list, definitions: dict) -> list:
    """Key entries (name, basis label, definition) in display order, for the report template."""
    items = []
    for name in names:
        basis, definition = definitions.get(name, (None, UNKNOWN_DEFINITION))
        items.append(
            {"name": name, "basis": BASIS_LABELS.get(basis), "definition": definition}
        )
    return items
LOG_TICKS = [0, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000]
PLOTLY_MODEBAR_REMOVE = ["select2d", "lasso2d", "select"]


def read_aligntrim_reports(directory: str = "aligntrim_qc") -> dict:
    """
    Read the align_trim QC TSVs staged in directory into {kind: {sample: DataFrame}}.
    Zero-byte files (stub runs) are skipped, header-only files give an empty DataFrame.
    """
    reports = {kind: {} for kind in ALIGNTRIM_REPORT_KINDS}
    for path in sorted(Path(directory).glob("*.tsv")):
        for kind in ALIGNTRIM_REPORT_KINDS:
            suffix = f".{kind}.tsv"
            if not path.name.endswith(suffix):
                continue
            if path.stat().st_size > 0:
                reports[kind][path.name[: -len(suffix)]] = pd.read_csv(
                    path, sep="\\t", keep_default_na=False, na_values=["NA"]
                )
            break
    return reports


def ordered_types(known: list, observed) -> list:
    """Known types in their canonical order, then any types a newer align_trim adds."""
    return list(known) + sorted(set(observed) - set(known))


def read_fate_colours(columns: list) -> dict:
    """Fixed colour per read fate, so the same fate is always the same colour."""
    unknown = sorted(c for c in columns if c not in READ_FATES)
    colours = {c: READ_FATES[c][0] for c in columns if c in READ_FATES}
    for i, c in enumerate(unknown):
        colours[c] = UNKNOWN_FATE_COLOURS[i % len(UNKNOWN_FATE_COLOURS)]
    return colours


def aligntrim_fates_bar(rows: list, columns: list) -> tuple[str, int]:
    """Horizontal stacked bar of the % of fragments in each fate, one bar per sample."""
    colours = read_fate_colours(columns)
    samples = [row["sample"] for row in rows]
    fig = go.Figure()
    for i, column in enumerate(columns):
        counts = [row["cells"][i]["n"] for row in rows]
        pcts = [row["cells"][i]["pct"] or 0.0 for row in rows]
        fig.add_trace(
            go.Bar(
                y=samples,
                x=pcts,
                name=column,
                orientation="h",
                marker_color=colours[column],
                customdata=counts,
                hovertemplate="%{y}<br>"
                + column
                + ": %{customdata:,} (%{x:.1f}%)<extra></extra>",
            )
        )
    # tall enough for the full legend to the right of the plot without scrolling
    height = max(160 + 28 * len(samples), 180 + 24 * len(columns))
    fig.update_layout(
        barmode="stack",
        height=height,
        margin=dict(t=60, b=60),
        font=dict(family="Roboto, monospace", size=14),
        hoverlabel=dict(font_family="Roboto, monospace"),
        title_text="Read Fates",
        legend=dict(
            orientation="v",
            x=1.01,
            xanchor="left",
            y=1,
            yanchor="top",
            traceorder="normal",
            font=dict(size=12),
        ),
        xaxis=dict(title="% of fragments", range=[0, 100]),
        yaxis=dict(type="category", autorange="reversed"),
        modebar_remove=PLOTLY_MODEBAR_REMOVE,
    )
    html = pio.to_html(
        fig,
        include_plotlyjs=False,
        full_html=False,
        default_height=f"{height}px",
        default_width="100%",
    )
    return html, height


def metric_heatmap(metrics: list, title: str, x_title: str) -> tuple[str, int]:
    """
    Heatmap (rows = samples) with a dropdown to switch between metrics.
    metrics: list of (label, DataFrame indexed by sample, scale) where scale is
    "log" (colour by log10(1 + n), ticks labelled with raw counts) or a (zmin, zmax) tuple.
    NaN cells are left blank.
    """
    fig = go.Figure()
    for i, (label, matrix, scale) in enumerate(metrics):
        raw = matrix.to_numpy(dtype=float)
        if scale == "log":
            z = np.log10(1 + raw)
            finite = z[np.isfinite(z)]
            zmax = max(float(finite.max()) if finite.size else 1.0, 1.0)
            ticks = [t for t in LOG_TICKS if np.log10(1 + t) <= zmax]
            zmin = 0.0
            colorbar = dict(
                title=label,
                tickvals=[float(np.log10(1 + t)) for t in ticks],
                ticktext=[f"{t:,}" for t in ticks],
            )
        else:
            z = raw
            zmin, zmax = scale
            colorbar = dict(title=label)
        # Pre-formatted hover values: 2D numpy customdata is not resolved by
        # heatmap hovertemplates in the bundled plotly.js, plain string text is
        hover_text = [
            [
                "NA"
                if not np.isfinite(v)
                else (f"{int(v):,}" if scale == "log" else f"{v:.4g}")
                for v in row
            ]
            for row in raw
        ]
        fig.add_trace(
            go.Heatmap(
                z=z,
                x=[str(c) for c in matrix.columns],
                y=[str(s) for s in matrix.index],
                text=hover_text,
                zmin=zmin,
                zmax=zmax,
                colorscale=plotly_colour_scale,
                colorbar=colorbar,
                xgap=1,
                ygap=1,
                hoverongaps=False,
                visible=(i == 0),
                hovertemplate="Sample: %{y}<br>"
                + x_title
                + ": %{x}<br>"
                + label
                + ": %{text}<extra></extra>",
            )
        )
    buttons = [
        dict(
            label=label,
            method="update",
            args=[{"visible": [j == i for j in range(len(metrics))]}],
        )
        for i, (label, _, _) in enumerate(metrics)
    ]
    n_rows = len(metrics[0][1].index) if metrics else 0
    height = max(350, 180 + 22 * n_rows)
    fig.update_layout(
        height=height,
        font=dict(family="Roboto, monospace", size=14),
        hoverlabel=dict(font_family="Roboto, monospace"),
        title_text=title,
        # metric dropdown only when there is more than one metric to switch between
        updatemenus=[
            dict(
                buttons=buttons,
                direction="down",
                showactive=True,
                x=1,
                xanchor="right",
                y=1.02,
                yanchor="bottom",
            )
        ]
        if len(metrics) > 1
        else [],
        modebar_remove=PLOTLY_MODEBAR_REMOVE,
    )
    fig.update_xaxes(type="category", title=x_title)
    fig.update_yaxes(type="category", autorange="reversed")
    html = pio.to_html(
        fig,
        include_plotlyjs=False,
        full_html=False,
        default_height=f"{height}px",
        default_width="100%",
    )
    return html, height


def chrom_display_name(chrom: str) -> str:
    header = primer_scheme.header_dict.get(chrom)
    return f"{header} ({chrom})" if header else chrom


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
    "nextclade_table": {},
    "run_warnings": [],
}

COVERAGE_DISAGREEMENT_TOLERANCE = 1.0  # percentage points


def amplicon_region_coverage(df, amplicons, min_depth) -> float:
    """
    Percentage of positions inside amplicon inserts (union of primer-trimmed regions,
    [coverage_start, coverage_end) 0-based half-open) with depth >= min_depth.
    df is a per-base coverage table with 1-based pos (chrom, pos, depth).
    """
    regions = {}
    for amp in amplicons:
        if amp.coverage_end > amp.coverage_start:
            regions.setdefault(str(amp.chrom), []).append(
                (amp.coverage_start, amp.coverage_end)
            )

    covered = 0
    total = 0
    for chrom, spans in regions.items():
        length = max(end for _, end in spans)
        mask = np.zeros(length, dtype=bool)
        for start, end in spans:
            mask[start:end] = True

        depth = np.zeros(length, dtype=int)
        chrom_df = df[df["chrom"].astype(str) == chrom]
        idx = chrom_df["pos"].to_numpy(dtype=int) - 1
        vals = chrom_df["depth"].to_numpy(dtype=int)
        in_range = (idx >= 0) & (idx < length)
        depth[idx[in_range]] = vals[in_range]

        covered += int(np.sum(depth[mask] >= min_depth))
        total += int(mask.sum())

    return round(covered / total * 100, 2) if total > 0 else 0.0

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

pre_norm_tsvs = glob("pre_norm_tsvs/*.tsv")
for tsv_path in pre_norm_tsvs:
    sample_name = tsv_path.split("/")[-1].split(".")[0]
    payload["qc_table_info"].setdefault(sample_name, {})
    df = pd.read_csv(tsv_path, sep="\t", index_col=None)
    payload["qc_table_info"][sample_name]["pre_norm_mean_depth"] = (
        round(df["depth"].mean(), 2) if len(df) > 0 else 0.0
    )
    payload["qc_table_info"][sample_name]["primer_scheme"] = scheme_version_str
    payload["qc_table_info"][sample_name]["pre_norm_coverage"] = (
        round(
            len(df[df["depth"] >= int("${params.min_coverage_depth}")]) / len(df) * 100,
            2,
        )
        if len(df) > 0
        else 0.0
    )
    payload["qc_table_info"][sample_name]["pre_norm_amplicon_coverage"] = (
        amplicon_region_coverage(
            df, primer_pairs, int("${params.min_coverage_depth}")
        )
    )
    coverage_for_qc =payload["qc_table_info"][sample_name]["pre_norm_coverage"]
    if coverage_for_qc >= int("${params.qc_pass_high_coverage}"):
        payload["qc_table_info"][sample_name]["qc_result"] = "pass"
    elif coverage_for_qc >= int("${params.qc_pass_min_coverage}"):
        payload["qc_table_info"][sample_name]["qc_result"] = "warning"
    else:
        payload["qc_table_info"][sample_name]["qc_result"] = "fail"

post_norm_tsvs = glob("post_norm_tsvs/*.tsv")
for tsv_path in post_norm_tsvs:
    sample_name = tsv_path.split("/")[-1].split(".")[0]
    payload["qc_table_info"].setdefault(sample_name, {})
    df = pd.read_csv(tsv_path, sep="\t", index_col=None)
    payload["qc_table_info"][sample_name]["post_norm_mean_depth"] = (
        round(df["depth"].mean(), 2) if len(df) > 0 else 0.0
    )
    payload["qc_table_info"][sample_name]["post_norm_coverage"] = (
        round(
            len(df[df["depth"] >= int("${params.min_coverage_depth}")]) / len(df) * 100,
            2,
        )
        if len(df) > 0
        else 0.0
    )
    payload["qc_table_info"][sample_name]["post_norm_amplicon_coverage"] = (
        amplicon_region_coverage(
            df, primer_pairs, int("${params.min_coverage_depth}")
        )
    )

has_post_norm = len(post_norm_tsvs) > 0

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
            # align_trim >= 2.0 reports mean_fragment_depth (a read pair counts once)
            row_dict = {
                "sample": sample_name,
                "chrom": row["chrom"],
                "amplicon": int(row["amplicon"]),
                "mean_fragment_depth": float(
                    row.get("mean_fragment_depth", row.get("mean_depth"))
                ),
            }
            # row_dict["amplicon"] = f"Amplicon {row_dict['amplicon']}"
            rows.append(row_dict)

        amplicon_depth_rows.extend(rows)

        covered = {(y["chrom"], y["amplicon"]) for y in rows}
        for x in primer_pairs:
            if (str(x.chrom), x.amplicon_number) not in covered:
                amplicon_depth_rows.append(
                    {
                        "sample": sample_name,
                        "chrom": str(x.chrom),
                        "amplicon": x.amplicon_number,
                        "mean_fragment_depth": 0.0,
                    }
                )

    payload["qc_table_info"][sample_name]["total_amp_dropouts"] = (
        len(
            [
                x
                for x in rows
                if float(x["mean_fragment_depth"]) < int("${params.min_coverage_depth}")
            ]
        )
        if len(rows) > 0
        else len(primer_pairs)
    )

# nextclade parsing tsv
for tsv in glob("nextclade_tsv/*.tsv"):
    with open(tsv, "r") as file:
        for row in csv.DictReader(file, delimiter="	"):
            sample_name = row.get("seqName", "").split(" ")[0]
            raw_coverage = row.get("coverage", "")
            payload["nextclade_table"][sample_name] = {
                "qc_status": row.get("qc.overallStatus", ""),
                "qc_score": row.get("qc.overallScore", ""),
                "clade": row.get("clade_display", ""),
                "lineage": row.get("Nextclade_pango", ""),
                "coverage": f"{round(float(raw_coverage) * 100, 2)}%" if raw_coverage else "",
                "qc_missing_status": row.get("qc.missingData.status", ""),
                "qc_missing_count": row.get("qc.missingData.totalMissing", ""),
                "qc_mixedsites_status": row.get("qc.mixedSites.status", ""),
                "qc_mixedsites_count": row.get("qc.mixedSites.totalMixedSites", ""),
                "qc_privatemut_status": row.get("qc.privateMutations.status", ""),
                "qc_privatemut_count": row.get("qc.privateMutations.total", ""),
                "qc_snpclust_status": row.get("qc.snpClusters.status", ""),
                "qc_snpclust_count": row.get("qc.snpClusters.totalSNPs", ""),
                "qc_framshift_status": row.get("qc.frameShifts.status", ""),
                "qc_framshift_count": row.get("qc.frameShifts.totalFrameShifts", ""),
                "qc_stopcodon_status": row.get("qc.stopCodons.status", ""),
                "qc_stopcodon_count": row.get("qc.stopCodons.totalStopCodons", ""),
            }

payload["nextclade_table"] = dict(
    sorted(
        payload["nextclade_table"].items(),
        key=lambda item: item[0],
    )
)

samples_with_amp_depths = {row["sample"] for row in amplicon_depth_rows}

for row in scheme_samplesheet_df.itertuples():
    if not payload["qc_table_info"].get(row.sample):
        samples.add(row.sample)
        payload["qc_table_info"].setdefault(row.sample, {})
        payload["qc_table_info"][row.sample]["primer_scheme"] = scheme_version_str
        payload["qc_table_info"][row.sample]["pre_norm_coverage"] = 0.0
        payload["qc_table_info"][row.sample]["pre_norm_mean_depth"] = 0.0
        payload["qc_table_info"][row.sample]["pre_norm_amplicon_coverage"] = 0.0
        if has_post_norm:
            payload["qc_table_info"][row.sample]["post_norm_coverage"] = 0.0
            payload["qc_table_info"][row.sample]["post_norm_mean_depth"] = 0.0
            payload["qc_table_info"][row.sample]["post_norm_amplicon_coverage"] = 0.0
        payload["qc_table_info"][row.sample]["total_reads"] = 0
        payload["qc_table_info"][row.sample]["total_amp_dropouts"] = len(primer_pairs)
        payload["qc_table_info"][row.sample]["qc_result"] = "fail"

    if row.sample not in samples_with_amp_depths:
        for x in primer_pairs:
            amplicon_depth_rows.append(
                {
                    "sample": row.sample,
                    "chrom": str(x.chrom),
                    "amplicon": x.amplicon_number,
                    "mean_fragment_depth": 0.0,
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

# Unified HTML coverage metric: amplicon-region coverage from the BAM used for
# consensus (post-normalisation when normalising, otherwise pre-normalisation).
# Pre and post values stay separate in the TSV.
coverage_disagreements = []
for sample, info in payload["qc_table_info"].items():
    pre_amp_cov = info.get("pre_norm_amplicon_coverage", 0.0)
    post_amp_cov = info.get("post_norm_amplicon_coverage")
    info["amplicon_coverage"] = (
        post_amp_cov if has_post_norm and post_amp_cov is not None else pre_amp_cov
    )
    if (
        post_amp_cov is not None
        and abs(pre_amp_cov - post_amp_cov) > COVERAGE_DISAGREEMENT_TOLERANCE
    ):
        coverage_disagreements.append(
            f"{sample} (pre {pre_amp_cov}% / post {post_amp_cov}%)"
        )

if coverage_disagreements:
    warning = (
        f"Pre- and post-normalisation amplicon coverage (>= ${params.min_coverage_depth}x) "
        f"differ by more than {COVERAGE_DISAGREEMENT_TOLERANCE}% for "
        f"{len(coverage_disagreements)} sample(s): {', '.join(coverage_disagreements)}."
    )
    try:
        normalise_depth = int("${params.normalise_depth}")
    except ValueError:
        normalise_depth = None
    if normalise_depth and normalise_depth < int("${params.min_coverage_depth}"):
        warning += (
            f" --normalise_depth ({normalise_depth}) is lower than --min_coverage_depth "
            "(${params.min_coverage_depth}), so normalisation is removing coverage needed "
            "to call consensus at these positions."
        )
    payload["run_warnings"].append(warning)

qc_fieldnames = [
    "sample",
    "primer_scheme",
    "pre_norm_coverage",
    "pre_norm_mean_depth",
]
if has_post_norm:
    qc_fieldnames += ["post_norm_coverage", "post_norm_mean_depth"]
qc_fieldnames += ["total_reads", "total_amp_dropouts", "qc_result"]
# New columns are appended so existing column positions are unchanged
qc_fieldnames += ["pre_norm_amplicon_coverage"]
if has_post_norm:
    qc_fieldnames += ["post_norm_amplicon_coverage"]

with open(f"{scheme_version_str.replace('/', '_')}_qc_results.tsv", "w") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=qc_fieldnames,
        delimiter="\t",
        extrasaction="ignore",
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

# Primer Trimming QC tab: summaries of the align_trim QC reports, values only
aligntrim_reports = read_aligntrim_reports()
aligntrim_samples = sorted(aligntrim_reports["summary"])
sample_platforms = (
    dict(zip(scheme_samplesheet_df["sample"], scheme_samplesheet_df["platform"]))
    if "platform" in scheme_samplesheet_df.columns
    else {}
)
payload["aligntrim"] = {
    "primer_threshold": PRIMER_THRESHOLD,
    "samples": aligntrim_samples,
    "no_report": sorted(
        set(scheme_samplesheet_df["sample"].astype(str)) - set(aligntrim_samples)
    ),
}

# Read fates (summary.tsv): kept and each removal reason as n and % of fragments
fate_columns = []
for df in aligntrim_reports["summary"].values():
    for column in df.columns:
        if (column == "kept" or column.startswith("removed_")) and column not in fate_columns:
            fate_columns.append(column)
fate_rows = []
for sample in aligntrim_samples:
    df = aligntrim_reports["summary"][sample]
    if len(df) == 0:
        continue
    record = df.iloc[0]
    fragments = int(record["fragments"])
    cells = []
    for column in fate_columns:
        n = int(record[column]) if column in df.columns else 0
        cells.append({"n": n, "pct": n / fragments * 100 if fragments else None})
    fate_rows.append(
        {
            "sample": sample,
            "platform": sample_platforms.get(sample, ""),
            "fragments": fragments,
            "cells": cells,
        }
    )
fate_labels = [column.replace("removed_", "") for column in fate_columns]
fate_colours = read_fate_colours(fate_labels)
payload["aligntrim"]["read_fates"] = {
    "columns": fate_labels,
    "rows": fate_rows,
    "plot": None,
    "key": [
        {
            "name": label,
            "colour": fate_colours[label],
            "basis": BASIS_LABELS.get(READ_FATES[label][1]) if label in READ_FATES else None,
            "definition": READ_FATES[label][2] if label in READ_FATES else UNKNOWN_DEFINITION,
        }
        for label in fate_labels
    ],
}
if fate_rows:
    fates_html, fates_height = aligntrim_fates_bar(fate_rows, fate_labels)
    payload["aligntrim"]["read_fates"]["plot"] = {
        "html": fates_html,
        "height": fates_height,
    }

# Amplicons (amplicons.tsv): per-amplicon fragment fates, one heatmap per chrom
amplicon_frames = [
    df.assign(sample=sample)
    for sample, df in aligntrim_reports["amplicons"].items()
    if len(df) > 0
]
payload["aligntrim"]["amplicon_plots"] = []
payload["aligntrim"]["amplicon_key"] = []
if amplicon_frames:
    amplicons_df = pd.concat(amplicon_frames, ignore_index=True)
    amplicon_metrics = [
        c
        for c in amplicons_df.columns
        if c.startswith("removed_") and c not in AMPLICON_HEATMAP_EXCLUDED
    ]
    # default metric first, the rest in align_trim column order
    amplicon_metrics.sort(key=lambda c: c != AMPLICON_HEATMAP_DEFAULT)
    payload["aligntrim"]["amplicon_key"] = key_items(
        [m.replace("removed_", "") for m in amplicon_metrics],
        {
            **{fate: (basis, text) for fate, (_, basis, text) in READ_FATES.items()},
            **AMPLICON_METRIC_DEFINITIONS,
        },
    )
    for chrom in sorted(amplicons_df["chrom"].astype(str).unique()):
        chrom_df = amplicons_df[amplicons_df["chrom"].astype(str) == chrom]
        metrics = []
        for metric in amplicon_metrics:
            matrix = chrom_df.pivot_table(
                index="sample", columns="amplicon", values=metric, aggfunc="sum"
            ).sort_index()
            matrix = matrix[sorted(matrix.columns, key=int)]
            metrics.append((metric.replace("removed_", ""), matrix, "log"))
        html, height = metric_heatmap(
            metrics,
            title=f"Removed Fragments per Amplicon: {chrom_display_name(chrom)}",
            x_title="Amplicon",
        )
        payload["aligntrim"]["amplicon_plots"].append({"html": html, "height": height})

# Primer pairs (primer_pairs.tsv): fragments per pair type
pair_frames = [df for df in aligntrim_reports["primer_pairs"].values() if len(df) > 0]
observed_pair_types = (
    set(pd.concat(pair_frames)["type"]) if pair_frames else set()
) - set(PRIMER_PAIR_EXTRA_TYPES)
pair_columns = ordered_types(PRIMER_PAIR_TYPES, observed_pair_types)
pair_rows = []
for sample in sorted(aligntrim_reports["primer_pairs"]):
    df = aligntrim_reports["primer_pairs"][sample]
    by_type = df.groupby("type")["fragments"].sum() if len(df) > 0 else pd.Series(dtype=int)
    total = int(sum(by_type.get(t, 0) for t in pair_columns))
    pair_rows.append(
        {
            "sample": sample,
            "total": total,
            "cells": [
                {
                    "n": int(by_type.get(t, 0)),
                    "pct": by_type.get(t, 0) / total * 100 if total else None,
                }
                for t in pair_columns
            ],
            "extra": [int(by_type.get(t, 0)) for t in PRIMER_PAIR_EXTRA_TYPES],
        }
    )
payload["aligntrim"]["primer_pairs"] = {
    "columns": pair_columns,
    "extra_columns": PRIMER_PAIR_EXTRA_TYPES,
    "rows": pair_rows,
    "key": key_items(pair_columns + PRIMER_PAIR_EXTRA_TYPES, PRIMER_PAIR_DEFINITIONS),
}

# Split reads (split_reads.tsv): fragment-joins and distinct joins per type
split_frames = [df for df in aligntrim_reports["split_reads"].values() if len(df) > 0]
split_columns = ordered_types(
    SPLIT_READ_TYPES, set(pd.concat(split_frames)["type"]) if split_frames else set()
)
split_rows = []
for sample in sorted(aligntrim_reports["split_reads"]):
    df = aligntrim_reports["split_reads"][sample]
    cells = []
    for split_type in split_columns:
        type_df = df[df["type"] == split_type] if len(df) > 0 else df
        cells.append(
            {
                "fragments": int(type_df["fragments"].sum()) if len(type_df) > 0 else 0,
                "joins": len(type_df),
            }
        )
    split_rows.append({"sample": sample, "cells": cells})
payload["aligntrim"]["split_reads"] = {
    "columns": split_columns,
    "rows": split_rows,
    "key": key_items(split_columns, SPLIT_READ_DEFINITIONS),
}

# Overlaps (overlaps.tsv): one heatmap per chrom, columns are amplicon n/n+1 overlaps
overlap_frames = [
    df.assign(sample=sample)
    for sample, df in aligntrim_reports["overlaps"].items()
    if len(df) > 0
]
overlap_metrics = [
    ("max_freq_difference", "max", (0, 1)),
]
payload["aligntrim"]["overlap_plots"] = []
payload["aligntrim"]["overlap_key"] = key_items(
    [metric for metric, _, _ in overlap_metrics], OVERLAP_METRIC_DEFINITIONS
)
if overlap_frames:
    overlaps_df = pd.concat(overlap_frames, ignore_index=True)
    overlaps_df["overlap"] = (
        overlaps_df["amplicon"].astype(str) + "/" + overlaps_df["next_amplicon"].astype(str)
    )
    for chrom in sorted(overlaps_df["chrom"].astype(str).unique()):
        chrom_df = overlaps_df[overlaps_df["chrom"].astype(str) == chrom]
        overlap_order = [
            label
            for _, label in sorted(
                {
                    (int(a), int(n)): f"{a}/{n}"
                    for a, n in zip(chrom_df["amplicon"], chrom_df["next_amplicon"])
                }.items()
            )
        ]
        metrics = []
        for metric, aggfunc, scale in overlap_metrics:
            matrix = chrom_df.pivot_table(
                index="sample",
                columns="overlap",
                values=metric,
                aggfunc=aggfunc,
                dropna=False,
            ).sort_index()
            metrics.append((metric, matrix.reindex(columns=overlap_order), scale))
        html, height = metric_heatmap(
            metrics,
            title=f"Amplicon Overlaps: {chrom_display_name(chrom)}",
            x_title="Overlap",
        )
        payload["aligntrim"]["overlap_plots"].append({"html": html, "height": height})

msa_list = glob("msas/*.fa*")
if len(msa_list) > 0:
    primer_mismatch_heatmaps = {"name": "Primer Mismatches", "plots": []}
    for msa_path in msa_list:
        msa, seqdict = parse_msa(msa_path) # type: ignore
        contig_name = msa_path.split("/")[-1].split("_")[0]

        primer_mismatch_heatmaps["plots"].append(
            {
                "name": contig_name,
                "plot_html": primer_mismatch_heatmap(
                    array=msa, seqdict=seqdict, bedfile="${bed}" # type: ignore
                ),
            }
        )
    payload["nested_plots"].append(primer_mismatch_heatmaps)

payload["has_post_norm"] = has_post_norm

render_qc_report(
    payload=payload,
    template_path=Path("${report_template}"),
    output_path=Path(
        f"{scheme_version_str.replace('/', '_')}_amplicon-nf_run-report.html"
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
