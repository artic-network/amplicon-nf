import plotly.express as px
import plotly.subplots as subplots
import plotly.graph_objects as go

import pandas as pd

import primalbedtools
import primalbedtools.bedfiles
import primalbedtools.primerpairs


def amplicon_depth_plot(amplicon_depth_tsv_path: str, min_depth: int = 20):

    df = pd.read_csv(amplicon_depth_tsv_path, sep="\t")

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
        figs[chrom].append(fig.to_html(full_html=False))

    return figs


def read_depth_plot(pileup_path: str, scheme_bed_path: str, min_depth: int = 20):

    header_lines, bed_lines = primalbedtools.bedfiles.BedLineParser.from_file(
        scheme_bed_path
    )

    primer_pairs = primalbedtools.primerpairs.create_primerpairs(bed_lines)

    scheme_df = pd.DataFrame(
        data=[
            [x.chrom, x.amplicon_number, x.pool, x.coverage_start, x.coverage_end]
            for x in primer_pairs
        ],
        columns=["chrom", "amplicon", "pool", "amplicon_start", "amplicon_end"],
    )

    pileup_df = pd.read_csv(pileup_path, sep="\t", names=["chrom", "pos", "depth"])

    chroms = scheme_df["chrom"].unique()
    figs = {}

    for chrom in chroms:
        pileup_df_chrom = pileup_df[pileup_df["chrom"] == chrom]
        chrom_primer_pairs = [x for x in primer_pairs if x.chrom == chrom]

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
                pileup_df_chrom,
                x="pos",
                y="depth",
                title=f"Read Depth for reference: {chrom}",
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
                    else pileup_df_chrom["pos"].max()
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
                name=pp.fbedlines[0].primername,
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
                name=pp.rbedlines[0].primername,
                row=2,
                col=1,
            )

        fig.add_trace(
            go.Scattergl(
                x=[x.coverage_start for x in chrom_primer_pairs],
                y=[x.pool for x in chrom_primer_pairs],
                opacity=0,
                name="Forward Primer",
                hovertext=[f"{x.fbedlines[0].primername}" for x in chrom_primer_pairs],
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
                hovertext=[f"{x.rbedlines[0].primername}" for x in chrom_primer_pairs],
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
            range=[0, pileup_df_chrom["pos"].max()],
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
            title="Amplicon Pools",
            tickmode="array",
            tickvals=pools,
            row=2,
            col=1,
        )

        figs.setdefault(chrom, [])
        # figs[chrom].append(fig.to_html(full_html=False))
        figs[chrom].append(fig)

    return figs


figs = read_depth_plot(
    pileup_path="/Users/sam/Library/CloudStorage/OneDrive-SharedLibraries-onedrive/bioinformatics/artic-repos/artic-network-fieldbioinformatics-nf/modules/local/generate_report/nanopore.pileup",
    scheme_bed_path="/Users/sam/bio_bulk/store_dir/fieldbioinformatics-nf/primer-schemes/artic-sars-cov-2/400/v3.0.0/primer.bed",
)

depth_figs = amplicon_depth_plot(
    amplicon_depth_tsv_path="/Users/sam/Library/CloudStorage/OneDrive-SharedLibraries-onedrive/bioinformatics/artic-repos/artic-network-fieldbioinformatics-nf/modules/local/generate_report/illumina_amplicon_test.amplicon_depths.tsv",
)

# figs["depth"] = depth_figs["MN908947.3"]

for chrom, fig_list in figs.items():
    for fig in fig_list:
        fig.show()

# import jinja2
# from jinja2 import Environment, FileSystemLoader
# import os

# templateLoader = jinja2.FileSystemLoader(
#     searchpath=os.path.dirname(os.path.realpath(__file__))
# )
# templateEnv = jinja2.Environment(loader=templateLoader)

# template = templateEnv.get_template("report_template.html")

# body_html = template.render(page_title="test-sample QC Report", figs=figs)

# with open("report.html", "w") as f:
#     f.write(body_html)
