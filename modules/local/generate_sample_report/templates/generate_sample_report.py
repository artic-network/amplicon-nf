#!/usr/bin/env python3

import plotly.express as px
import plotly.subplots as subplots
import plotly.graph_objects as go
import plotly.io as pio

import pandas as pd

import csv
from importlib.metadata import version

from primalbedtools.scheme import Scheme
from primalbedtools.amplicons import create_amplicons


from pathlib import Path
from datetime import datetime
from jinja2 import Environment, FileSystemLoader, select_autoescape


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


def amplicon_depth_plot(amplicon_depth_tsv_path: str, min_depth: int = 20):

    df = pd.read_csv(amplicon_depth_tsv_path, sep="\\t")

    chroms = df["chrom"].unique()

    figs = {}

    for chrom in chroms:
        df_chrom = df[df["chrom"] == chrom]
        fig = px.bar(
            df_chrom,
            x="amplicon",
            y="mean_depth",
            title=f"Amplicon Depths for reference: {chrom}",
            labels={"amplicon": "Amplicon", "mean_depth": "Mean Depth"},
        )
        fig.add_hline(y=min_depth, line_color="red", line_dash="dash")

        fig.update_xaxes(
            ticks="outside",
            tickcolor="black",
            showticklabels=True,
            linecolor="black",
            gridcolor="lightgrey",
        )
        fig.update_yaxes(
            ticks="outside",
            tickcolor="black",
            showticklabels=True,
            linecolor="black",
            gridcolor="lightgrey",
        )
        fig.update_layout(plot_bgcolor="whitesmoke")
        fig.update_traces(marker_color="skyblue")

        figs.setdefault(chrom, [])
        figs[chrom].append(fig)

    return figs


def read_depth_plot(
    depth_df: pd.DataFrame,
    scheme_df: pd.DataFrame,
    primer_pairs: list,
    min_depth: int = 20,
):

    chroms = scheme_df["chrom"].unique()
    figs = {}

    for chrom in chroms:
        depth_df_chrom = depth_df[depth_df["chrom"] == chrom]
        chrom_primer_pairs = [x for x in primer_pairs if x.chrom == chrom]
        chrom_alias = scheme_df[scheme_df["chrom"] == chrom]["chrom_alias"].values[0]

        chrom_label = f"{chrom_alias} ({chrom})" if chrom_alias else chrom

        fig = subplots.make_subplots(
            cols=1,
            rows=2,
            shared_xaxes=True,
            row_heights=[4, 0.5],
            specs=[
                [{"secondary_y": True}],
                [{"secondary_y": True}],
            ],
            vertical_spacing=0.02,
        )

        fig.add_trace(
            px.line(
                depth_df_chrom,
                x="pos",
                y="depth",
                title=f"Read Depth for reference: {chrom_label}",
                labels={"pos": "Position", "depth": "Read Depth"},
            ).data[0],
            row=1,
            col=1,
        )
        fig.add_hline(y=min_depth, line_color="red", line_dash="dash", row=1, col=1)
        fig.update_xaxes(
            ticks="",
            tickcolor="black",
            showticklabels=False,
            linecolor="black",
            gridcolor="lightgrey",
            row=1,
            col=1,
        )
        fig.update_yaxes(
            ticks="outside",
            tickcolor="black",
            showticklabels=True,
            linecolor="black",
            gridcolor="lightgrey",
            row=1,
            col=1,
            title="Read Depth",
        )

        for pp in chrom_primer_pairs:
            fig.add_shape(
                type="line",
                y0=pp.pool,
                y1=pp.pool,
                x0=pp.amplicon_start,
                x1=(
                    pp.amplicon_end
                    if not pp.is_circular
                    else depth_df_chrom["pos"].max()
                ),
                line=dict(color="LightSeaGreen", width=5),
                name=f"amplicon {pp.amplicon_number}",
                row=2,
                col=1,
            )
            # Handle circular genomes
            if pp.is_circular:
                fig.add_shape(
                    type="line",
                    y0=pp.pool,
                    y1=pp.pool,
                    x0=0,
                    x1=pp.amplicon_end,
                    line=dict(color="LightSeaGreen", width=5),
                    name=f"amplicon {pp.amplicon_number}",
                    row=2,
                    col=1,
                )

            fig.add_shape(
                type="rect",
                y0=pp.pool - 0.05,
                y1=pp.pool + 0.05,
                x0=pp.amplicon_start,
                x1=pp.coverage_start,
                fillcolor="LightSalmon",
                line=dict(color="darksalmon", width=3),
                name=pp.left[0].primername,
                row=2,
                col=1,
            )
            fig.add_shape(
                type="rect",
                y0=pp.pool - 0.05,
                y1=pp.pool + 0.05,
                x0=pp.coverage_end,
                x1=pp.amplicon_end,
                fillcolor="LightSalmon",
                line=dict(color="darksalmon", width=3),
                name=pp.right[0].primername,
                row=2,
                col=1,
            )

        fig.add_trace(
            go.Scattergl(
                x=[x.coverage_start for x in chrom_primer_pairs],
                y=[x.pool for x in chrom_primer_pairs],
                opacity=0,
                name="Forward Primer",
                hovertext=[f"{x.left[0].primername}" for x in chrom_primer_pairs],
                showlegend=False,
            ),
            row=2,
            col=1,
        )
        fig.add_trace(
            go.Scattergl(
                x=[x.coverage_end for x in chrom_primer_pairs],
                y=[x.pool for x in chrom_primer_pairs],
                opacity=0,
                name="Reverse Primer",
                hovertext=[f"{x.right[0].primername}" for x in chrom_primer_pairs],
                showlegend=False,
            ),
            row=2,
            col=1,
        )

        fig.update_layout(plot_bgcolor="whitesmoke")

        fig.update_xaxes(
            showline=True,
            mirror=True,
            ticks="outside",
            linewidth=2,
            linecolor="black",
            tickformat=",d",
            range=[0, depth_df_chrom["pos"].max()],
            title="",
        )
        fig.update_xaxes(
            title="Position",
            row=2,
            col=1,
        )
        fig.update_yaxes(
            showline=True,
            mirror=True,
            ticks="outside",
            linewidth=2,
            linecolor="black",
            fixedrange=True,
        )
        pools = sorted({x.pool for x in chrom_primer_pairs})
        fig.update_yaxes(
            range=[pools[0] - 0.5, pools[-1] + 0.5],
            # title="Amplicon Pools",
            tickmode="array",
            tickvals=pools,
            row=2,
            col=1,
        )

        figs[chrom] = fig

    return figs


with open("${coverage_report}", "rt") as f:
    reader = csv.DictReader(f, delimiter="\\t")
    reads = {x["#rname"]: x["numreads"] for x in reader}

scheme = Scheme.from_file("${bed}")

scheme_headers = scheme.header_dict

primer_pairs = create_amplicons(scheme.bedlines)

scheme_df = pd.DataFrame(
    data=[
        [
            x.chrom,
            scheme_headers.get(x.chrom),
            x.amplicon_number,
            x.pool,
            x.coverage_start,
            x.coverage_end,
        ]
        for x in primer_pairs
    ],
    columns=[
        "chrom",
        "chrom_alias",
        "amplicon",
        "pool",
        "amplicon_start",
        "amplicon_end",
    ],
)

depth_df = pd.read_csv("${depth_tsv}", sep="\\t", names=["chrom", "pos", "depth"])

with open("${amp_depth_tsv}", "rt") as f:
    reader = csv.DictReader(f, delimiter="\\t")
    amplicon_depths = [
        {
            "chrom": row["chrom"],
            "amplicon": row["amplicon"],
            "mean_depth": float(row["mean_depth"]),
        }
        for row in reader
    ]

    for x in primer_pairs:
        if (
            len(
                [
                    y
                    for y in amplicon_depths
                    if y["chrom"] == str(x.chrom)
                    and y["amplicon"] == str(x.amplicon_number)
                ]
            )
            == 0
        ):
            amplicon_depths.append(
                {
                    "chrom": str(x.chrom),
                    "amplicon": str(x.amplicon_number),
                    "mean_depth": 0.0,
                }
            )

amplicon_depths.sort(key=lambda x: (x["chrom"], int(x["amplicon"])))

plot = read_depth_plot(
    depth_df=depth_df,
    scheme_df=scheme_df,
    primer_pairs=primer_pairs,
    min_depth=int("${params.min_coverage_depth}"),
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
    "sequencing_platform": "${meta.platform}",
    "tool_name": "amplicon-nf",
    "tool_version": "${workflow.manifest.version}",
    "citation_link": "https://doi.org/10.5281/zenodo.17522200",
    "contact_email": "s.a.j.wilkinson@bham.ac.uk",
    "funder_statement": "This pipeline has been created as part of the ARTIC network project funded by the Wellcome Trust (collaborator award – 313694/Z/24/Z and discretionary award – 206298/Z/17/Z) and is distributed as open source and open access. All non-code files are made available under a Creative Commons CC-BY licence unless otherwise specified. Please acknowledge or cite this repository or associated publications if used in derived work so we can provide our funders with evidence of impact in the field.",
    "timestamp": datetime.now().strftime("%Y-%m-%d_%H:%M:%S"),
    "contigs": [],
}

for chrom, fig in plot.items():

    contig_plot_html = pio.to_html(
        fig,
        include_plotlyjs=False,
        full_html=False,
        default_height="1500px",
        default_width="100%",
    )

    bases_above_min_depth = depth_df[
        (depth_df["chrom"] == chrom)
        & (depth_df["depth"] >= int("${params.min_coverage_depth}"))
    ].shape[0]
    total_bases = depth_df[depth_df["chrom"] == chrom].shape[0]

    percent_coverage = (
        (bases_above_min_depth / total_bases) * 100 if total_bases > 0 else 0.0
    )

    amplicon_dropouts = [
        str(x["amplicon"])
        for x in amplicon_depths
        if x["chrom"] == chrom and x["mean_depth"] < int("${params.min_coverage_depth}")
    ]

    mean_depth = (
        round(depth_df[depth_df["chrom"] == chrom]["depth"].mean(), 2)
        if not depth_df[depth_df["chrom"] == chrom].empty
        else 0.0
    )

    contig_alias = scheme_df[scheme_df["chrom"] == chrom]["chrom_alias"].values[0]
    contig_label = f"{contig_alias} ({chrom})" if contig_alias else chrom

    # payload["contigs"].append(
    contig_payload = {
        "name": chrom,
        "contig_alias": contig_alias,
        "contig_label": contig_label,
        "percent_coverage": round(percent_coverage, 2),
        "amplicon_dropouts": amplicon_dropouts,
        "total_amplicon_dropouts": len(amplicon_dropouts),
        "amplicon_mean_depths": {
            x["amplicon"]: round(x["mean_depth"], 2)
            for x in amplicon_depths
            if x["chrom"] == chrom
        },
        "amplicon_coords": {
            str(x.amplicon_number): {
                "start": str(x.coverage_start),
                "end": str(x.coverage_end),
                "pool": str(x.pool),
            }
            for x in primer_pairs
            if x.chrom == chrom
        },
        "qc_status": "",
        "total_reads": reads[chrom],
        "average_depth": mean_depth,
        "plotly_html": contig_plot_html,
    }

    if contig_payload["percent_coverage"] >= int("${params.qc_pass_high_coverage}"):
        contig_payload["qc_status"] = "pass"
    elif contig_payload["percent_coverage"] >= int("${params.qc_pass_min_coverage}"):
        contig_payload["qc_status"] = "warning"
    else:
        contig_payload["qc_status"] = "fail"

    payload["contigs"].append(contig_payload)


render_qc_report(
    payload=payload,
    template_path=Path("${report_template}"),
    output_path=Path(f"{"${meta.id}"}_amplicon-nf_sample-report.html"),
    bootstrap_css_path=Path("${bootstrap_bundle_min_css}"),
    bootstrap_bundle_js_path=Path("${bootstrap_bundle_min_js}"),
    plotly_js_path=Path("${plotly_js}"),
    svg_logo_path=Path("${artic_logo_svg}"),
)

with open("versions.yml", "w") as f:
    f.write("${task.process}:\\n")
    for package in ("plotly", "primalbedtools", "pandas", "jinja2"):
        f.write(f"  {package}: {version(package)}\\n")
