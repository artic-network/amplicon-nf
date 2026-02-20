import argparse
import os
import re
import sys
from importlib.resources import files

import pandas as pd
import plotly.graph_objects as go

# from bokeh.sampledata.autompg2 import autompg2
# from bokeh.layouts import gridplot
# from bokeh.models import ColumnDataSource, NumeralTickFormatter, Whisker
# from bokeh.plotting import figure
# from bokeh.resources import CDN
from jinja2 import Template
from plotly.offline import plot
from plotly.subplots import make_subplots

# from bokeh.embed import components


class A4:
    width = 794
    height = 1123 * 1.3


class Colors:
    """Some colours that someone thought were nice."""

    cerulean = "#0084A9"
    not_black = "#001A21"
    feldgrau = "#455556"
    dim_gray = "#666666"
    light_cornflower_blue = "#90C5E7"
    dark_gray = "#B5AEA7"
    isabelline = "#F0EFED"
    medium_spring_bud = "#B8E986"
    cinnabar = "#EF4134"
    sandstorm = "#F5CC49"
    fandango = "#A53F96"
    green = "#17BB75"
    verdigris = "#54B8B1"


def get_barcode(x):
    pattern = r"/(barcode\d+)/"
    match = re.search(pattern, x)
    if match:
        barcode = match.group(1)
        return barcode
    else:
        return "Unknown"


def plot_reads_count_plotly(
    df, sample_col, top=1200, bottom=50, title="Test reads Count"
):
    # Filter and count reads
    df_count = (
        df[(df.read_length >= bottom) & (df.read_length <= top)]
        .groupby(sample_col)
        .size()
        .reset_index(name="reads_count")
        .sort_values(by=sample_col)
    )

    fig = go.Figure(
        go.Bar(
            x=df_count[sample_col], y=df_count["reads_count"], marker_color="steelblue"
        )
    )

    fig.update_layout(
        title=title,
        xaxis_title="Sample Name",
        yaxis_title="Number of Reads",
        # height=350,
        # width=1200,
        xaxis_tickangle=-90,
        yaxis=dict(rangemode="tozero", tickformat=",.0f"),
        bargap=0.01,
    )

    return fig


# # for bokeh
# def plot_reads_count(df, sample_col, top=1200, bottom=50, title="Test reads Count"):
#     df_count = df[((df.read_length >= bottom) & (df.read_length <= top))][
#         [sample_col, "read_length"]
#     ]
#     df_count = df_count.groupby(by=sample_col).count()
#     df_count = df_count.reset_index()
#     df_count.columns = [sample_col, "reads_count"]

#     df_count = df_count.sort_values(by=sample_col)

#     # print(df_count.head())
#     p = figure(
#         x_range=df_count[sample_col],
#         title=title,
#         toolbar_location=None,
#         tools="",
#         y_axis_label="Number of Reads",
#         x_axis_label="Sample Name",
#     )

#     p.vbar(x=df_count[sample_col], top=df_count["reads_count"], width=0.9)

#     p.xgrid.grid_line_color = None
#     p.y_range.start = 0
#     p.xaxis.major_label_orientation = 1.5708
#     # p.xgrid.grid_line_color = None
#     p.axis.major_label_text_font_size = "14px"
#     p.axis.axis_label_text_font_size = "12px"
#     p.yaxis.formatter = NumeralTickFormatter(format="0.0a")

#     return p


def plot_box_plotly(
    df,
    sample_col,
    value_col,
    yaxis="value",
    title="Unknown",
    color="blue",
    outlier=True,
):
    """
    Convert Bokeh boxplot to Plotly boxplot.
    Plotly’s built-in box trace already computes quartiles, whiskers and outliers,
    so we only need to group the data and draw one trace per sample.
    """
    # Ensure categorical order
    df[sample_col] = df[sample_col].astype(str)
    df = df.sort_values(by=sample_col)
    samples = df[sample_col].unique()

    fig = go.Figure()

    for sample in samples:
        y = df.loc[df[sample_col] == sample, value_col].values
        fig.add_trace(
            go.Box(
                y=y,
                name=sample,
                boxmean=False,  # no mean line
                marker_color=color,
                line_color="grey",
                showlegend=False,
                boxpoints="outliers" if outlier else False,
                hovertemplate="Median: %{median}<extra></extra>",
                pointpos=0,
                jitter=0,
            )
        )

    fig.update_layout(
        title=title,
        xaxis_title="Sample Name",
        yaxis_title=yaxis,
        # height=300,
        # width=1200,
        xaxis=dict(tickangle=-90),
        yaxis=dict(rangemode="tozero"),
        bargap=0.01,
    )

    return fig


# def plot_box_bokeh(
#     df,
#     sample_col,
#     value_col,
#     yaxis="value",
#     title="Unknown",
#     color="blue",
#     outlier=True,
# ):
#     df[sample_col] = df[sample_col].astype(str)
#     df = df.sort_values(by=sample_col)
#     samples = df[sample_col].unique()
#     grouper = df.groupby(by=sample_col)
#     qs = grouper[value_col].quantile([0.25, 0.5, 0.75]).unstack().reset_index()
#     qs.columns = [sample_col, "q1", "q2", "q3"]

#     # compute IQR outlier bounds
#     iqr = qs.q3 - qs.q1
#     qs["upper"] = qs.q3 + 1.5 * iqr
#     qs["lower"] = qs.q1 - 1.5 * iqr
#     qs["lower"] = qs["lower"].apply(lambda x: x if x > 1 else 1)
#     for sample, group in grouper:
#         qs_idx = qs[qs[sample_col] == sample].index[0]
#         data = group[value_col]

#         # the upper whisker is the maximum between p3 and upper
#         q3 = qs.loc[qs_idx, "q3"]
#         upper = qs.loc[qs_idx, "upper"]
#         wiskhi = group[(q3 <= data) & (data <= upper)][value_col]
#         qs.loc[qs_idx, "upper"] = q3 if len(wiskhi) == 0 else wiskhi.max()

#         # the lower whisker is the minimum between q1 and lower
#         q1 = qs.loc[qs_idx, "q1"]
#         lower = qs.loc[qs_idx, "lower"]
#         wisklo = group[(lower <= data) & (data <= q1)][value_col]
#         if value_col == "read_length":
#             qs.loc[qs_idx, "lower"] = (
#                 q1 + 10 if len(wisklo) == 0 else (wisklo.min() + 10)
#             )
#         else:
#             qs.loc[qs_idx, "lower"] = q1 if len(wisklo) == 0 else wisklo.min()

#     df = pd.merge(df, qs, on=sample_col, how="left")

#     qs = qs.sort_values(by=sample_col)
#     source = ColumnDataSource(qs)

#     p = figure(
#         x_range=samples,
#         tools="",
#         toolbar_location=None,
#         title=title,
#         # height=300,
#         # width=1200,
#         y_axis_label=yaxis,
#         x_axis_label="Sample Name",
#     )

#     # outlier range
#     whisker = Whisker(base=sample_col, upper="upper", lower="lower", source=source)
#     whisker.upper_head.size = whisker.lower_head.size = 15
#     p.add_layout(whisker)

#     # quantile boxes
#     # cmap = factor_cmap("kind", "TolRainbow7", kinds)
#     p.vbar(sample_col, 0.7, "q2", "q3", source=source, color=color, line_color="black")
#     p.vbar(sample_col, 0.7, "q1", "q2", source=source, color=color, line_color="black")

#     # outliers
#     if outlier:
#         outliers = df[~df[value_col].between(df.lower, df.upper)]
#         p.scatter(
#             sample_col, value_col, source=outliers, size=6, color=color, alpha=0.3
#         )

#     p.xaxis.major_label_orientation = 1.5708
#     # p.xgrid.grid_line_color = None
#     p.axis.major_label_text_font_size = "14px"
#     p.axis.axis_label_text_font_size = "12px"
#     p.y_range.start = 0
#     p.y_range.end = df.upper.max() * 1.1

#     return p


# def plot_summary(reads_stats_path):
#     colors = Colors()
#     # load the data frame
#     df = pd.DataFrame()
#     try:
#         df = pd.read_csv(reads_stats_path, sep="\t")
#     except Exception:
#         print(
#             "Error to read fastcat per reads stats file, please check your input or fastcat running stats"
#         )
#         return None, None, None

#     # make the plots
#     df_read = df[["sample_name", "read_length"]].copy()
#     df_quality = df[["sample_name", "mean_quality"]].copy()
#     p1 = plot_box_bokeh(
#         df_read,
#         "sample_name",
#         "read_length",
#         yaxis="Reads Length",
#         title="Boxplot for reads length",
#         color=colors.fandango,
#         outlier=False,
#     )
#     p2 = plot_box_bokeh(
#         df_quality,
#         "sample_name",
#         "mean_quality",
#         yaxis="Qscore",
#         title="Boxplot for reads Quality",
#         color=colors.cerulean,
#         outlier=False,
#     )
#     p3 = plot_reads_count(df, "sample_name", title="Bar Plot for Reads Count")
#     return (p1, p2, p3)


def plot_summary_plotly(reads_stats_path, colors=Colors):
    """
    Plotly-only version of plot_summary.
    Returns three Plotly Figure objects:
      fig_len  – box plot of read lengths
      fig_qual – box plot of mean quality
      fig_cnt  – bar plot of read counts
    """

    df = pd.read_csv(reads_stats_path, sep="\t")
    if "sample_name" not in list(df.columns):
        df["sample_name"] = df["filename"].apply(get_barcode)

    def _filter_outlier(df, sample_col, value_col):
        grouper = df.groupby(by=sample_col)
        qs = grouper[value_col].quantile([0.25, 0.5, 0.75]).unstack().reset_index()
        qs.columns = [sample_col, "q1", "q2", "q3"]
        iqr = qs.q3 - qs.q1
        qs["upper"] = qs.q3 + 1.5 * iqr
        df = pd.merge(df, qs, on=sample_col, how="left")
        df = df[df[value_col] <= df["upper"]]
        return df

    def _box_plotly(df_in, sample_col, value_col, yaxis, title, color):

        a4 = A4()
        # print(df_in.head())
        df_in = df_in[df_in[sample_col] != "Unknown"]
        df_in = _filter_outlier(df_in, sample_col, value_col)
        df_in[sample_col] = df_in[sample_col].astype(str)
        df_in = df_in.sort_values(by=sample_col)
        samples = df_in[sample_col].unique()
        fig = go.Figure()
        for s in samples:
            fig.add_trace(
                go.Box(
                    y=df_in.loc[df_in[sample_col] == s, value_col].tolist(),
                    name=s,
                    boxmean=False,
                    marker_color=color,
                    line_color="black",
                    line_width=1,
                    fillcolor=color,
                    showlegend=False,
                    hoverinfo="skip",
                )
            )
        fig.update_traces(boxpoints=False)
        fig.update_layout(
            title=title,
            xaxis_title="Sample Name",
            yaxis_title=yaxis,
            xaxis_tickangle=-90,
            plot_bgcolor="white",
            paper_bgcolor="white",
            yaxis=dict(
                rangemode="tozero",
                gridcolor="#d3d3d3",
                gridwidth=0.3,
            ),
            bargap=0.01,
            bargroupgap=0.01,
            # autosize=True,
            width=a4.width,
            height=a4.height / 4,
        )
        return fig

    # ---- barplot helper ----
    def _bar_plotly(df_in, sample_col, title):
        a4 = A4()
        df_in = df_in[df_in[sample_col] != "Unknown"]
        counts = (
            df_in.groupby(sample_col)
            .size()
            .reset_index(name="reads_count")
            .sort_values(by=sample_col)
        )
        # print(counts.head())
        fig = go.Figure(
            go.Bar(
                x=counts[sample_col].tolist(),
                y=counts["reads_count"].tolist(),
                marker_color="steelblue",
                hovertemplate="%{y}<extra></extra>",
                showlegend=False,
            )
        )
        fig.update_layout(
            title=title,
            xaxis_title="Sample Name",
            yaxis_title="Number of Reads",
            xaxis_tickangle=-90,
            plot_bgcolor="white",
            paper_bgcolor="white",
            yaxis=dict(
                rangemode="tozero",
                gridcolor="#d3d3d3",
                gridwidth=0.3,
            ),
            bargap=0.01,
            width=a4.width,
            height=a4.height / 4,
        )
        return fig

    # build figures
    fig_len = _box_plotly(
        df,
        "sample_name",
        "read_length",
        yaxis="Read Length",
        title="Boxplot for Read Length",
        color=colors.fandango,
    )
    # fig_len.write_image("fig_len.png")
    fig_qual = _box_plotly(
        df,
        "sample_name",
        "mean_quality",
        yaxis="Qscore",
        title="Boxplot for Read Quality",
        color=colors.cerulean,
    )
    # fig_qual.write_image("fig_quanl.png")
    fig_cnt = _bar_plotly(df, "sample_name", title="Bar Plot for Read Count")
    # fig_cnt.write_image("fig_cnt.png")

    return fig_len, fig_qual, fig_cnt


def plot_coverage_plotly(
    bed_path: str, threshold: float, xlim: float, ylim: float, ncols: int, colors
):
    """
    Returns a single Plotly Figure containing the grid of coverage plots.
    """
    if not os.path.exists(bed_path):
        print(f"{bed_path} does not exist, please double-check your input")
        sys.exit(1)

    cols = ["chrome", "start", "end", "depth", "pool", "sample"]
    try:
        df = pd.read_csv(bed_path, sep="\t", header=None, names=cols)
        if df.iloc[0, 0] == "chrome":  # header row present
            df = pd.read_csv(bed_path, sep="\t")
    except Exception:
        print(f"Error loading bed file {bed_path}; check format")
        sys.exit(1)

    df["pos"] = (df["start"] + df["end"]) / 2

    samples = sorted(df["sample"].unique())
    nplots = len(samples)
    nrows = (nplots + ncols - 1) // ncols

    fig = make_subplots(
        rows=nrows,
        cols=ncols,
        subplot_titles=[f"sample{i}" for i in range(nplots)],
        shared_xaxes=False,
        shared_yaxes=False,
        horizontal_spacing=0.1,
        vertical_spacing=0.09,
    )
    # a4 = A4()

    for idx, sample in enumerate(samples):
        row = (idx // ncols) + 1
        col = (idx % ncols) + 1

        df_sub = df[df["sample"] == sample].copy()
        df_sum = (
            df_sub[["start", "end", "depth", "pos"]]
            .groupby(["start", "end", "pos"])
            .sum()
            .reset_index()
        )
        mean_depth = df_sum["depth"].mean()
        pass_ratio = 100 * (df_sum["depth"] >= threshold).sum() / len(df_sum)
        title = f"{sample}: {mean_depth:.0f}X, {pass_ratio:.1f}% > {threshold}X"

        # pool-1 area + line
        p1 = df_sub[df_sub["pool"] == 1]
        fig.add_trace(
            go.Scatter(
                x=p1["pos"].tolist(),
                y=p1["depth"].tolist(),
                mode="lines",
                line=dict(color=colors.dark_gray),
                showlegend=False,
                hovertemplate="Position: %{x}<br>Pool-1: %{y}<extra></extra>",
            ),
            row=row,
            col=col,
        )
        fig.add_trace(
            go.Scatter(
                x=p1["pos"].tolist(),
                y=p1["depth"].tolist(),
                mode="none",
                fill="tozeroy",
                fillcolor=colors.dark_gray,
                showlegend=False,
                hoverinfo="skip",
            ),
            row=row,
            col=col,
        )

        # pool-2 area + line
        p2 = df_sub[df_sub["pool"] == 2]
        fig.add_trace(
            go.Scatter(
                x=p2["pos"].tolist(),
                y=p2["depth"].tolist(),
                mode="lines",
                line=dict(color=colors.verdigris),
                showlegend=False,
                hovertemplate="Position: %{x}<br>Pool-2: %{y}<extra></extra>",
            ),
            row=row,
            col=col,
        )
        fig.add_trace(
            go.Scatter(
                x=p2["pos"].tolist(),
                y=p2["depth"].tolist(),
                mode="none",
                fill="tozeroy",
                fillcolor=colors.verdigris,
                showlegend=False,
                hoverinfo="skip",
            ),
            row=row,
            col=col,
        )

        fig.update_xaxes(range=[0, xlim], title="position", row=row, col=col)
        fig.update_yaxes(range=[0, ylim], title="depth", row=row, col=col)
        fig.update_annotations(font_size=10)
        fig.layout.annotations[idx].text = title  # set subplot title

    fig.update_layout(height=300 * nrows, width=250 * ncols)
    return fig


# def plot_coverage_bokeh(
#     bed_path: str, threshold: float, xlim: float, ylim: float, ncols: int, colors
# ):
#     ## PART C plot the coverage
#     if not os.path.exists(bed_path):
#         print(f"{bed_path} does not exist, please double check your input")
#         sys.exit()
#     # read the bed file with pandas
#     df = pd.DataFrame()
#     cols = ["chrome", "start", "end", "depth", "pool", "sample"]
#     try:
#         df = pd.read_csv(bed_path, sep="\t", header=None, names=cols)
#         # print(df.iloc[0,0])
#         if df.iloc[0, 0] == "chrome":
#             df = pd.read_csv(bed_path, sep="\t")
#     except Exception:
#         print(
#             f"Error when loading the bed file {bed_path}, please check the input format"
#         )
#         sys.exit()
#     df["pos"] = (df["start"] + df["end"]) / 2

#     plot_list = []
#     for sample in list(df["sample"].unique()):
#         df_sub = df[df["sample"] == sample]
#         df_sum = (
#             df_sub[["start", "end", "depth", "pos"]]
#             .groupby(by=["start", "end", "pos"])
#             .sum()
#         )
#         mean_depth = df_sum.depth.mean()
#         pass_ratio = 100 * (df_sum.depth >= threshold).sum() / len(df_sum.depth)
#         title = "{}: {:.0f}X, {:.1f}% > {}X".format(
#             "test", mean_depth, pass_ratio, threshold
#         )
#         source = ColumnDataSource(
#             data=dict(
#                 depth_pool_1=df_sub[df_sub.pool == 1]["depth"],
#                 depth_pool_2=df_sub[df_sub.pool == 2]["depth"],
#                 position=df_sub.pos.unique(),
#             )
#         )
#         TOOLTIPS = [
#             ("Position", "@position"),
#             ("Pool-1", "@depth_pool_1"),
#             ("Pool-2", "@depth_pool_2"),
#         ]
#         p = figure(
#             width=400,
#             height=250,
#             tooltips=TOOLTIPS,
#             title=title,
#             x_axis_label="position",
#             y_axis_label="depth",
#             x_range=(0, xlim),
#             y_range=(0, ylim),
#         )
#         p.line(x="position", y="depth_pool_1", color=Colors.dark_gray, source=source)
#         p.varea(
#             x="position",
#             y1=0,
#             y2="depth_pool_1",
#             source=source,
#             fill_color=Colors.dark_gray,
#             alpha=0.7,
#             muted_color=Colors.dark_gray,
#             muted_alpha=0.2,
#         )
#         p.line(x="position", y="depth_pool_2", source=source, color=Colors.verdigris)
#         p.varea(
#             x="position",
#             y1=0,
#             y2="depth_pool_2",
#             source=source,
#             fill_color=Colors.verdigris,
#             alpha=0.7,
#             muted_color=Colors.verdigris,
#             muted_alpha=0.2,
#         )
#         plot_list.append(p)

#     coverage_plot_bokeh = gridplot(plot_list, ncols=ncols)
#     print("-- coverage plots are created successfully -- ")
#     return coverage_plot_bokeh


def combine_and_save_plotly(
    figures: list,
    outfile: str,
    ncols: int = 2,
    width: int = 1200,
    height: int = 600,
    **kwargs,
):
    """
    Combine a list of Plotly Figure objects into a single HTML file.
    figures   : iterable of go.Figure (or compatible)
    outfile   : path for the saved HTML
    ncols     : number of columns in the grid
    width/height: dimensions of the final layout
    kwargs    : passed to make_subplots (e.g. subplot_titles, shared_xaxes, …)
    """
    nplots = len(figures)
    nrows = (nplots + ncols - 1) // ncols

    # Build the container grid
    big = make_subplots(rows=nrows, cols=ncols, **kwargs)

    # Copy each figure’s traces into the appropriate subplot cell
    for idx, fig in enumerate(figures):
        row = (idx // ncols) + 1
        col = (idx % ncols) + 1
        for tr in fig.data:
            big.add_trace(tr, row=row, col=col)

        # Apply original axes ranges / labels
        big.update_xaxes(fig.layout.xaxis, row=row, col=col)
        big.update_yaxes(fig.layout.yaxis, row=row, col=col)

    big.update_layout(width=width, height=height)
    big.write_html(outfile)
    print(f"-- Saved combined plot to {os.path.abspath(outfile)} --")


def load_template():
    package_root = files("amp_depth_viz")
    template_path = package_root / "template" / "template_plotly.html"
    mytemplate = os.path.abspath(template_path)
    with open(mytemplate, "r") as f:
        template_content = f.read()

    template = Template(template_content)
    return template


# use plotly to make plots and use jinja2 to render for html template
def plotly_plot_combine(plot_list: list, template, output_path: str):

    plot_divs = []
    for p in plot_list:
        plot_divs.append(plot(p, output_type="div", include_plotlyjs=False))

    html_output = template.render({"divs": plot_divs})
    with open(output_path, "w") as f:
        f.write(html_output)
    print(f"-- Saved combined plotly plots to {output_path} --")


def main():
    parser = argparse.ArgumentParser(
        description="create wf-artic like amplicon coverage plots using bokeh"
    )
    parser.add_argument(
        "coveragebed",
        help="a full sample bed files contains all genome depth information from mosdepth",
    )
    # input either the per-read-stats.tsv or the fastq_pass folder
    group = parser.add_mutually_exclusive_group()

    # Add arguments to the exclusive group
    group.add_argument(
        "--fastcat_perreads", help="per-read-stats.tsv output from fastcat"
    )
    group.add_argument("--fastq_pass", help="the fastq_folder from a sars cov run")

    # Add template parameter
    parser.add_argument(
        "--template",
        help="jinja2 template to render -- default in template/template.html",
    )
    parser.add_argument(
        "--output", default="amplicon_coverage.html", help="html output file"
    )
    parser.add_argument(
        "--xlim", default=30000, type=int, help="max of the genome position"
    )
    parser.add_argument("--ylim", default=800, type=int, help="the expected max depth")
    parser.add_argument(
        "--threshold",
        default=20,
        type=int,
        help="depth threshold for passing QC , default as 20",
    )
    parser.add_argument(
        "--threads",
        default=10,
        type=int,
        help="max number of cpus to use for fastcat analysis",
    )
    parser.add_argument(
        "--ncols",
        default=3,
        type=int,
        help="number of columns to grid the plots in the html pages, default as 3",
    )

    parser.add_argument(
        "--plot_package",
        "-p",
        type=str,
        choices=["bokeh", "plotly"],
        default="bokeh",
        help="Choose a python plot package to make the plots,from the allowed options: %(choices)s",
    )
    args = parser.parse_args()

    ## PART A Read the template
    if args.template:
        mytemplate = os.path.abspath(args.template)
    else:
        package_root = files("amp_depth_viz")
        template_path = package_root / "template" / f"template_{args.plot_package}.html"
        mytemplate = os.path.abspath(template_path)

    ## load template
    with open(mytemplate, "r") as f:
        template_content = f.read()

    template = Template(template_content)
    output_path = os.path.abspath(args.output)

    print("-- template loaded successfully --")

    ## PART B.1 Run fastcat
    if args.fastq_pass:
        # read_stats_path = run_fastcat_on_folder(args.fastq_pass, args.threads)
        pass
    else:
        read_stats_path = os.path.abspath(args.fastcat_perreads)

    if args.plot_package == "bokeh":
        pass
    #     ## PART B.2 Plot summary
    #     p1, p2, p3 = plot_summary(read_stats_path)
    #     # fig_len, fig_qual, fig_cnt = plot_summary_plotly(read_stats_path)

    #     if p1 is None:
    #         print(
    #             "Something wrong with the summary plot, please check errors and input"
    #         )
    #         sys.exit()
    #     # if fig_len is None:
    #     # print("Something wrong with the summary plot, please check errors and input")
    #     # sys.exit()
    #     print("-- reads summary plots are created successfully --")

    #     ### HERE PART C FUNC GOES
    #     coverage_plot_bokeh = plot_coverage_bokeh(
    #         args.coveragebed,
    #         threshold=args.threshold,
    #         xlim=args.xlim,
    #         ylim=args.ylim,
    #         ncols=args.ncols,
    #         colors=Colors(),
    #     )

    #     ## PART D templete render to create html
    #     script, divs = components((p1, p2, p3, coverage_plot_bokeh))

    #     html = template.render(
    #         script=script,
    #         divs=divs,  # Pass the tuple or dict
    #         resources=CDN.render(),  # Optional: loads Bokeh JS/CSS from CDN
    #     )

    #     # Save to file or return in a web route
    #     print(f"-- Saving the plot to {output_path} --")
    #     with open(output_path, "w") as f:
    #         f.write(html)
    #     print("-- amp-depth-viz pipeline ran successfully without error -- ")

    elif args.plot_package == "plotly":
        print("-- making plotly plots --")
        fig_len, fig_qual, fig_cnt = plot_summary_plotly(read_stats_path)
        if fig_len is None:
            print(
                "Something wrong with the summary plot, please check errors and input"
            )
            sys.exit()
        coverage_plot_plotly = plot_coverage_plotly(
            args.coveragebed,
            threshold=args.threshold,
            xlim=args.xlim,
            ylim=args.ylim,
            ncols=args.ncols,
            colors=Colors(),
        )

        plotly_plot_combine(
            [fig_len, fig_qual, fig_cnt, coverage_plot_plotly],
            template=template,
            output_path=output_path,
        )
        print("-- amp-depth-viz pipeline ran successfully without error -- ")

        # print(" saving plotly plot ")
        # combine_and_save_plotly([fig_len, fig_qual, fig_cnt], "test-plotly.html", ncols=1)
        # combine_and_save_plotly([fig_len, fig_qual, fig_cnt], "test-plotly.html", ncols=1)
