import os
import re
import subprocess
import sys
from importlib.resources import files

import pandas as pd
from bokeh.embed import components
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource, NumeralTickFormatter, Whisker
from bokeh.plotting import figure
from bokeh.resources import CDN
from jinja2 import Template


class _colors:
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


def bash_command(cmd):
    p = subprocess.Popen(cmd, shell=True)
    while True:
        return_code = p.poll()
        if return_code is not None:
            break
    return


# run fastcat on the fastq_pass folder
def run_fastcat_on_folder(folder, threads):
    # pass
    print("-- running fastcat on fastq_pass folder, may take a while --")
    folder_pass = os.path.abspath(folder)
    cmd = "rm -rf histograms"
    bash_command(cmd)
    cwd = os.getcwd()
    cmd = f"fastcat -x -t {threads} -f per-file-stats.tsv   -i per-file-runids.tsv -l \
    per-file-basecallers.tsv -r per-read-stats.tsv --histograms histograms\
      {folder_pass} > log.txt"
    bash_command(cmd)
    cmd = "rm -rf log.txt"
    bash_command(cmd)
    print("-- fastcat finish running --")
    return os.path.join(cwd, "per-read-stats.tsv")


def get_barcode(x):
    pattern = r"/(barcode\d+)/"
    match = re.search(pattern, x)
    if match:
        barcode = match.group(1)
        return barcode
    else:
        return "Unknown"


def plot_reads_count(df, sample_col, top=1200, bottom=50, title="Test reads Count"):
    df_count = df[((df.read_length >= bottom) & (df.read_length <= top))][
        [sample_col, "read_length"]
    ]
    df_count = df_count.groupby(by=sample_col).count()
    df_count = df_count.reset_index()
    df_count.columns = [sample_col, "reads_count"]

    df_count = df_count.sort_values(by=sample_col)

    # print(df_count.head())
    p = figure(
        x_range=df_count[sample_col],
        height=350,
        width=1200,
        title=title,
        toolbar_location=None,
        tools="",
        y_axis_label="Number of Reads",
        x_axis_label="SAMPLE",
    )

    p.vbar(x=df_count[sample_col], top=df_count["reads_count"], width=0.9)

    p.xgrid.grid_line_color = None
    p.y_range.start = 0
    p.xaxis.major_label_orientation = 1.5708
    # p.xgrid.grid_line_color = None
    p.axis.major_label_text_font_size = "14px"
    p.axis.axis_label_text_font_size = "12px"
    p.yaxis.formatter = NumeralTickFormatter(format="0.0a")

    return p


def plot_box_bokeh(
    df,
    sample_col,
    value_col,
    yaxis="value",
    title="Unknown",
    color="blue",
    outlier=True,
):
    df[sample_col] = df[sample_col].astype(str)
    df = df.sort_values(by=sample_col)
    samples = df[sample_col].unique()
    grouper = df.groupby(by=sample_col)
    qs = grouper[value_col].quantile([0.25, 0.5, 0.75]).unstack().reset_index()
    qs.columns = [sample_col, "q1", "q2", "q3"]

    # compute IQR outlier bounds
    iqr = qs.q3 - qs.q1
    qs["upper"] = qs.q3 + 1.5 * iqr
    qs["lower"] = qs.q1 - 1.5 * iqr
    qs["lower"] = qs["lower"].apply(lambda x: x if x > 1 else 1)
    for sample, group in grouper:
        qs_idx = qs[qs[sample_col] == sample].index[0]
        data = group[value_col]

        # the upper whisker is the maximum between p3 and upper
        q3 = qs.loc[qs_idx, "q3"]
        upper = qs.loc[qs_idx, "upper"]
        wiskhi = group[(q3 <= data) & (data <= upper)][value_col]
        qs.loc[qs_idx, "upper"] = q3 if len(wiskhi) == 0 else wiskhi.max()

        # the lower whisker is the minimum between q1 and lower
        q1 = qs.loc[qs_idx, "q1"]
        lower = qs.loc[qs_idx, "lower"]
        wisklo = group[(lower <= data) & (data <= q1)][value_col]
        if value_col == "read_length":
            qs.loc[qs_idx, "lower"] = (
                q1 + 10 if len(wisklo) == 0 else (wisklo.min() + 10)
            )
        else:
            qs.loc[qs_idx, "lower"] = q1 if len(wisklo) == 0 else wisklo.min()

    df = pd.merge(df, qs, on=sample_col, how="left")

    qs = qs.sort_values(by=sample_col)
    source = ColumnDataSource(qs)

    p = figure(
        x_range=samples,
        tools="",
        toolbar_location=None,
        title=title,
        height=300,
        width=1200,
        y_axis_label=yaxis,
        x_axis_label="SAMPLE",
    )

    # outlier range
    whisker = Whisker(base=sample_col, upper="upper", lower="lower", source=source)
    whisker.upper_head.size = whisker.lower_head.size = 15
    p.add_layout(whisker)

    # quantile boxes
    # cmap = factor_cmap("kind", "TolRainbow7", kinds)
    p.vbar(sample_col, 0.7, "q2", "q3", source=source, color=color, line_color="black")
    p.vbar(sample_col, 0.7, "q1", "q2", source=source, color=color, line_color="black")

    # outliers
    if outlier:
        outliers = df[~df[value_col].between(df.lower, df.upper)]
        p.scatter(
            sample_col, value_col, source=outliers, size=6, color=color, alpha=0.3
        )

    p.xaxis.major_label_orientation = 1.5708
    # p.xgrid.grid_line_color = None
    p.axis.major_label_text_font_size = "14px"
    p.axis.axis_label_text_font_size = "12px"
    p.y_range.start = 0
    p.y_range.end = df.upper.max() * 1.1

    return p


def plot_summary(reads_stats_path, samplesheet=None):
    Colors = _colors()
    # load the data frame
    df = pd.DataFrame()
    try:
        df = pd.read_csv(reads_stats_path, sep="\t")
    except:
        print(
            "Error to read fastcat per reads stats file, please check your input or fastcat running stats"
        )
        return None, None, None
    df["barcode"] = df["filename"].apply(get_barcode)
    with_samplesheet = False
    df_sheet_short = pd.DataFrame()
    if samplesheet:
        df_sheet = pd.read_csv(samplesheet)
        try:
            df_sheet_short = df_sheet[["sampleID", "barcode"]]
            with_samplesheet = True
        except:
            print(
                "Did not find sampleID or barcode columns in the samplesheet csv file, please check your input"
            )
            print("Will ignore the input samplesheet and use barcode as sampleID")
    if with_samplesheet:
        df = df.merge(df_sheet_short, how="left", on="barcode")
        df = df.dropna(subset=["sampleID"])
    else:
        df["sampleID"] = df["barcode"]
        df = df[df["sampleID"] != "Unknown"]

    # make the plots
    df_read = df[["sampleID", "read_length"]].copy()
    df_quality = df[["sampleID", "mean_quality"]].copy()
    p1 = plot_box_bokeh(
        df_read,
        "sampleID",
        "read_length",
        yaxis="Reads Length",
        title="Boxplot for reads length",
        color=Colors.fandango,
        outlier=False,
    )
    p2 = plot_box_bokeh(
        df_quality,
        "sampleID",
        "mean_quality",
        yaxis="Qscore",
        title="Boxplot for reads Quality",
        color=Colors.cerulean,
        outlier=False,
    )
    p3 = plot_reads_count(df, "sampleID", title="Bar Plot for Reads Count")
    return (p1, p2, p3)


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
        "--samplesheet",
        help="An optional samplesheet to rename the barcode, at least two columns required: sampleID, barcode",
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
    args = parser.parse_args()

    Colors = _colors()
    ## PART A Read the template
    mytemplate = ""
    if args.template:
        mytemplate = os.path.abspath(args.template)
    else:
        package_root = files("amp_depth_viz")
        template_path = package_root / "template" / "template.html"
        mytemplate = os.path.abspath(template_path)

    # print(mytemplate)

    ## load template
    with open(mytemplate, "r") as f:
        template_content = f.read()

    template = Template(template_content)

    print("-- template loaded successfully --")
    # sys.exit()

    # html = template.render(
    # script=script,
    # divs=divs,  # Pass the tuple or dict
    # resources=CDN.render()  # Optional: loads Bokeh JS/CSS from CDN
    # )
    # Save to file or return in a web route
    # with open("multiple_bokeh_2.html", "w") as f:
    # f.write(html)

    ## PART B.1 Run fastcat
    read_stats_path = ""
    if args.fastq_pass:
        read_stats_path = run_fastcat_on_folder(args.fastq_pass, args.threads)
    else:
        read_stats_path = os.path.abspath(args.fastcat_perreads)

    # print(read_stats_path)
    df = pd.read_csv(read_stats_path, sep="\t")
    ## PART B.2 Plot summary
    p1, p2, p3 = None, None, None
    if args.samplesheet:
        sheetpath = os.path.abspath(args.samplesheet)
        p1, p2, p3 = plot_summary(read_stats_path, sheetpath)
    else:
        p1, p2, p3 = plot_summary(read_stats_path)

    if p1 == None:
        print("Something wrong with the summary plot, please check errors and input")
        sys.exit()

    print("-- reads summary plots are created successfully --")

    ## PART C plot the coverage
    ## check the inputs
    bed_path = args.coveragebed

    if not os.path.exists(bed_path):
        print(f"{bed_path} does not exist, please double check your input")
        sys.exit()
    # read the bed file with pandas
    df = pd.DataFrame()
    cols = ["chrome", "start", "end", "depth", "pool", "sample"]
    try:
        df = pd.read_csv(bed_path, sep="\t", header=None, names=cols)
        # print(df.iloc[0,0])
        if df.iloc[0, 0] == "chrome":
            df = pd.read_csv(bed_path, sep="\t")
    except:
        print(
            f"Error when loading the bed file {bed_path}, please check the input format"
        )
        sys.exit()
    df["pos"] = (df["start"] + df["end"]) / 2

    # Colors = _colors()
    plot_list = []
    for sample in list(df["sample"].unique()):
        df_sub = df[df["sample"] == sample]
        df_sum = (
            df_sub[["start", "end", "depth", "pos"]]
            .groupby(by=["start", "end", "pos"])
            .sum()
        )
        mean_depth = df_sum.depth.mean()
        pass_ratio = 100 * (df_sum.depth >= args.threshold).sum() / len(df_sum.depth)
        title = "{}: {:.0f}X, {:.1f}% > {}X".format(
            "test", mean_depth, pass_ratio, args.threshold
        )
        source = ColumnDataSource(
            data=dict(
                depth_pool_1=df_sub[df_sub.pool == 1]["depth"],
                depth_pool_2=df_sub[df_sub.pool == 2]["depth"],
                position=df_sub.pos.unique(),
            )
        )
        TOOLTIPS = [
            ("Position", "@position"),
            ("Pool-1", "@depth_pool_1"),
            ("Pool-2", "@depth_pool_2"),
        ]
        p = figure(
            width=400,
            height=250,
            tooltips=TOOLTIPS,
            title=title,
            x_axis_label="position",
            y_axis_label="depth",
            x_range=(0, args.xlim),
            y_range=(0, args.ylim),
        )
        p.line(x="position", y="depth_pool_1", color=Colors.dark_gray, source=source)
        p.varea(
            x="position",
            y1=0,
            y2="depth_pool_1",
            source=source,
            fill_color=Colors.dark_gray,
            alpha=0.7,
            muted_color=Colors.dark_gray,
            muted_alpha=0.2,
        )
        p.line(x="position", y="depth_pool_2", source=source, color=Colors.verdigris)
        p.varea(
            x="position",
            y1=0,
            y2="depth_pool_2",
            source=source,
            fill_color=Colors.verdigris,
            alpha=0.7,
            muted_color=Colors.verdigris,
            muted_alpha=0.2,
        )
        plot_list.append(p)

    coverage_plot_bokeh = gridplot(plot_list, ncols=args.ncols)
    # output_file(filename=args.output, title="Amplicon Coverage")
    # save(grid_plot)
    print("-- coverage plots are created successfully -- ")
    ## PART D templete render to create html
    script, divs = components((p1, p2, p3, coverage_plot_bokeh))

    html = template.render(
        script=script,
        divs=divs,  # Pass the tuple or dict
        resources=CDN.render(),  # Optional: loads Bokeh JS/CSS from CDN
    )
    # Save to file or return in a web route
    output_path = os.path.abspath(args.output)
    print(f"-- Saving the plot to {output_path} --")
    with open(output_path, "w") as f:
        f.write(html)
    print("-- amp-depth-viz pipeline ran successfully without error -- ")
