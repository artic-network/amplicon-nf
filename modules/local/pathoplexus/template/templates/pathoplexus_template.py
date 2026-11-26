#!/usr/bin/env python3
#
# Per-scheme half of the Pathoplexus submission template (on unless --skip_pathoplexus_template).
# Combines the PATHOPLEXUS_SAMPLE rows into <prefix>_pathoplexus_metadata.{xlsx,tsv} and
# <prefix>_pathoplexus_sequences.fasta, laid out like the Pathoplexus website template
# (sheets Data / _lists (hidden) / Guidance).
#
# NOTE: this file is a Nextflow template. Every literal dollar sign and backslash would be
# interpreted by Nextflow, so the script deliberately uses neither outside the block of
# Nextflow values below (chr() is used for tab, newline and the Excel absolute-reference anchor).

import json
import platform
import sys
import tomllib
from datetime import datetime
from pathlib import Path

import openpyxl
from openpyxl import Workbook
from openpyxl.comments import Comment
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.datavalidation import DataValidation

# ---- values injected by Nextflow ----------------------------------------------------------
PREFIX = "${prefix}"
SCHEME_LABEL = "${meta.scheme ?: meta.custom_scheme_name}"
PIPELINE_NAME = "${workflow.manifest.name}"
PIPELINE_VERSION = "${pipeline_version}"
TEMPLATE_CONFIG = "${template_config}"
TASK_PROCESS = "${task.process}"
# -------------------------------------------------------------------------------------------

TAB = chr(9)
NL = chr(10)
CR = chr(13)
DOLLAR = chr(36)
REQUIRED_FIELDS = ["id", "sampleCollectionDate", "geoLocCountry", "authors"]
VALIDATION_ROWS = 100000  # dropdowns cover rows users add as well, as on the Pathoplexus template

TIER_FILL = {"required": "FCE4D6", "desired": "DDEBF7", "optional": "F2F2F2"}
USER_HEADER_FILL = "FFC000"
USER_EMPTY_FILL = {"required": "F8CBAD", "desired": "FFF2CC", "optional": "FFF2CC"}
AUTO_FONT = Font(color="595959", italic=True)
HEADER_FONT = Font(bold=True)
THIN = Side(style="thin", color="BFBFBF")
BOX = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)


def fill(colour):
    return PatternFill(fill_type="solid", start_color=colour, end_color=colour)


def clean(value):
    if value is None:
        return ""
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return value
    return str(value).replace(TAB, " ").replace(NL, " ").replace(CR, " ").strip()


def read_fasta_ids(path):
    with open(path, "rt") as fh:
        return [line[1:].split()[0] for line in fh if line.startswith(">") and line[1:].split()]


# ---- config --------------------------------------------------------------------------------
with open(TEMPLATE_CONFIG, "rb") as fh:
    config = tomllib.load(fh)
fields = config["fields"]
names = [f["name"] for f in fields]
duplicates = sorted({n for n in names if names.count(n) > 1})
missing = [n for n in REQUIRED_FIELDS if n not in names]
if duplicates or missing:
    sys.exit(f"Invalid Pathoplexus template config: duplicate fields {duplicates}, missing required fields {missing}")
provenance = config.get("provenance", {})

# ---- samples -------------------------------------------------------------------------------
samples = []
for path in sorted(Path("samples").glob("*.json")):
    with open(path) as fh:
        samples.append(json.load(fh))
samples.sort(key=lambda s: s["id"])
included = [s for s in samples if not s["excluded"]]
excluded = [s for s in samples if s["excluded"]]
auto_fields = set().union(*(s.get("auto_fields", []) for s in samples)) | {"id"}

multi = any(len(s["fasta_ids"]) > 1 for s in included)
has_reads = any(s["row"].get("files.rawReads") for s in included)
rows = []
for s in included:
    row = dict(s["row"])
    if multi:
        row["fastaIds"] = " ".join(s["fasta_ids"])
    rows.append(row)
if multi:
    auto_fields.add("fastaIds")

columns = [
    f for f in fields if not (f["name"] == "fastaIds" and not multi) and not (f["name"] == "files.rawReads" and not has_reads)
]

problems = []
for f in columns:
    options = set(f.get("options", []))
    if not options:
        continue
    for row in rows:
        value = row.get(f["name"], "")
        if value != "" and value not in options:
            problems.append(f"{row['id']}: {f['name']} value '{value}' is not a Pathoplexus vocabulary option")
for message in problems:
    print(f"WARNING [pathoplexus] {message}", file=sys.stderr)

# ---- FASTA (checked against the metadata so the upload links up exactly) ----------------
fasta_by_sample = {}
for path in Path("fastas").glob("*.pathoplexus.fasta"):
    fasta_by_sample[path.name[: -len(".pathoplexus.fasta")]] = path

written_ids = []
with open(f"{PREFIX}_pathoplexus_sequences.fasta", "w") as out:
    for s in included:
        path = fasta_by_sample.get(s["id"])
        if path is None:
            sys.exit(f"Missing Pathoplexus FASTA for sample {s['id']}")
        text = path.read_text()
        written_ids += read_fasta_ids(path)
        out.write(text if text.endswith(NL) or not text else text + NL)

expected_ids = [i for s in included for i in s["fasta_ids"]]
if sorted(written_ids) != sorted(expected_ids) or len(set(written_ids)) != len(written_ids):
    sys.exit("Pathoplexus FASTA IDs do not match the metadata (id / fastaIds), refusing to write a template that would be rejected")
if not multi and any(s["fasta_ids"] != [s["id"]] for s in included):
    sys.exit("Single-sequence entries must use the sample id as FASTA ID")

# ---- TSV -----------------------------------------------------------------------------------
with open(f"{PREFIX}_pathoplexus_metadata.tsv", "w", newline="") as out:
    out.write(TAB.join(f["name"] for f in columns) + NL)
    for row in rows:
        out.write(TAB.join(str(clean(row.get(f["name"], ""))) for f in columns) + NL)

# ---- XLSX ----------------------------------------------------------------------------------
wb = Workbook()
data = wb.active
data.title = "Data"
lists = wb.create_sheet("_lists")
guidance = wb.create_sheet("Guidance")

list_column = 0
last_row = max(len(rows) + 1, 2)
for col_idx, f in enumerate(columns, start=1):
    name, tier = f["name"], f.get("tier", "optional")
    letter = get_column_letter(col_idx)
    is_user = name not in auto_fields

    cell = data.cell(row=1, column=col_idx, value=name)
    cell.font = HEADER_FONT
    cell.fill = fill(USER_HEADER_FILL if is_user else TIER_FILL.get(tier, TIER_FILL["optional"]))
    cell.border = BOX
    note = [f"{f.get('display_name', name)} ({tier})"]
    note.append("PLEASE FILL IN" if is_user else "Filled in automatically by amplicon-nf, please check")
    for key in ("definition", "guidance"):
        if f.get(key):
            note.append(f[key])
    if f.get("example"):
        note.append(f"Example: {f['example']}")
    cell.comment = Comment(NL + NL.join(note), "amplicon-nf", width=400, height=300)

    values = [clean(row.get(name, "")) for row in rows]
    is_int = bool(values) and all(isinstance(v, int) for v in values if v != "") and any(v != "" for v in values)
    number_format = "0" if is_int else "@"
    data.column_dimensions[letter].number_format = number_format
    for row_idx, value in enumerate(values, start=2):
        c = data.cell(row=row_idx, column=col_idx, value=value if value != "" else None)
        c.number_format = number_format
        if is_user and value == "":
            c.fill = fill(USER_EMPTY_FILL.get(tier, USER_EMPTY_FILL["optional"]))
        elif not is_user:
            c.font = AUTO_FONT

    longest = max([len(name)] + [len(str(v)) for v in values])
    data.column_dimensions[letter].width = 24 if is_user else min(max(longest + 2, 12), 60)

    if f.get("options"):
        list_column += 1
        list_letter = get_column_letter(list_column)
        lists.cell(row=1, column=list_column, value=name)
        for i, option in enumerate(f["options"], start=2):
            lists.cell(row=i, column=list_column, value=option)
        end = len(f["options"]) + 1
        formula = f"=_lists!{DOLLAR}{list_letter}{DOLLAR}2:{DOLLAR}{list_letter}{DOLLAR}{end}"
        validation = DataValidation(type="list", formula1=formula, allow_blank=True, showErrorMessage=True)
        validation.error = f"Please choose a value from the {name} list"
        validation.errorTitle = "Invalid value"
        data.add_data_validation(validation)
        validation.add(f"{letter}2:{letter}{VALIDATION_ROWS}")

data.freeze_panes = "B2"
lists.sheet_state = "hidden"
lists.protection.sheet = True

# ---- Guidance sheet ------------------------------------------------------------------------
wrap = Alignment(wrap_text=True, vertical="top")
user_columns = [f for f in columns if f["name"] not in auto_fields]
lines = [
    ("title", f"Pathoplexus submission template: {SCHEME_LABEL}"),
    ("text", f"Generated by {PIPELINE_NAME} v{PIPELINE_VERSION} on {datetime.now().strftime('%Y-%m-%d %H:%M')}."),
    ("blank", ""),
    ("heading", "Please fill in"),
]
for f in user_columns:
    lines.append(("text", f"- {f['name']} ({f.get('display_name', f['name'])}, {f.get('tier', 'optional')})"))
lines += [
    ("text", "Highlighted (orange header) columns on the Data sheet are left for you; sampleCollectionDate, geoLocCountry and authors are required for every organism."),
    ("text", "host is REQUIRED for dengue, yellow fever, zika and chikungunya submissions (processing fails without it) and strongly recommended for all others."),
    ("text", "You can add any other Pathoplexus metadata field as an extra column using its exact field name."),
    ("blank", ""),
    ("heading", "Before uploading"),
    ("text", "- Check the auto-filled values (grey italic text)."),
    ("text", "- Samples that failed QC are included with their QC result filled in. To leave a sample out, delete its row here, its record(s) in the FASTA file and its file(s) in raw_reads/."),
    ("text", "- This template assumes one target organism per primer scheme (segmented genomes and reference-selecting schemes are supported). Primer panels and multi-pathogen schemes are NOT supported."),
    ("blank", ""),
    ("heading", "Uploading"),
    ("text", f"- Website: upload this XLSX (or the TSV) as metadata and {PREFIX}_pathoplexus_sequences.fasta as sequences."),
]
if has_reads:
    lines.append(("text", "- Website raw reads: use 'Upload folder' and select the raw_reads folder next to this file; it holds exactly the files listed in files.rawReads."))
    lines.append(("text", "- API: submit the TSV and FASTA; raw reads must first be uploaded via /files/request-upload and files.rawReads replaced with fileName:fileId pairs."))
else:
    lines.append(("text", "- Raw reads were not prepared (dehosting was skipped or no dehosted reads were available)."))
lines.append(("text", "- Do not add extra sheets: Pathoplexus reads the Data sheet and ignores Guidance and _lists."))
if excluded:
    lines += [("blank", ""), ("heading", "Samples left out of this template")]
    for s in excluded:
        lines.append(("text", f"- {s['id']}: {s.get('exclusion_reason') or 'excluded'}"))
if problems:
    lines += [("blank", ""), ("heading", "Values to check")]
    lines += [("text", f"- {p}") for p in problems]
lines += [("blank", ""), ("heading", "Colour key")]

row_idx = 1
for kind, text in lines:
    c = guidance.cell(row=row_idx, column=1, value=text if text else None)
    if kind == "title":
        c.font = Font(bold=True, size=14)
    elif kind == "heading":
        c.font = Font(bold=True, size=12)
    row_idx += 1

key = [
    (USER_HEADER_FILL, "Header: column for you to fill in"),
    (TIER_FILL["required"], "Header: required by Pathoplexus"),
    (TIER_FILL["desired"], "Header: desired by Pathoplexus"),
    (TIER_FILL["optional"], "Header: optional"),
    (USER_EMPTY_FILL["required"], "Cell: required value missing"),
    (USER_EMPTY_FILL["desired"], "Cell: recommended value missing"),
]
for colour, text in key:
    guidance.cell(row=row_idx, column=1).fill = fill(colour)
    guidance.cell(row=row_idx, column=2, value=text)
    row_idx += 1

row_idx += 1
guidance.cell(row=row_idx, column=1, value="Fields").font = Font(bold=True, size=12)
row_idx += 1
headers = ["Field", "Display name", "Requirement", "Filled by", "Definition", "Guidance", "Example", "Allowed values"]
for col_idx, text in enumerate(headers, start=1):
    c = guidance.cell(row=row_idx, column=col_idx, value=text)
    c.font = HEADER_FONT
    c.border = BOX
row_idx += 1
for f in columns:
    options = f.get("options", [])
    allowed = "; ".join(options) if len(options) <= 20 else (f"{len(options)} values, use the dropdown on the Data sheet" if options else "")
    values = [
        f["name"],
        f.get("display_name", ""),
        f.get("tier", ""),
        "you" if f["name"] not in auto_fields else "amplicon-nf",
        f.get("definition", ""),
        f.get("guidance", ""),
        f.get("example", ""),
        allowed,
    ]
    for col_idx, value in enumerate(values, start=1):
        c = guidance.cell(row=row_idx, column=col_idx, value=value or None)
        c.alignment = wrap
        c.border = BOX
    row_idx += 1

row_idx += 1
guidance.cell(row=row_idx, column=1, value="Provenance").font = Font(bold=True, size=12)
row_idx += 1
for label, value in [
    ("Pipeline", f"{PIPELINE_NAME} v{PIPELINE_VERSION}"),
    ("Primer scheme", SCHEME_LABEL),
    ("Pathoplexus spec (pathoplexus commit)", provenance.get("pathoplexus_sha", "")),
    ("Loculus version", provenance.get("loculus_version", "")),
    ("Spec checked on", provenance.get("checked_on", "")),
]:
    guidance.cell(row=row_idx, column=1, value=label)
    guidance.cell(row=row_idx, column=2, value=str(value))
    row_idx += 1

for letter, width in zip("ABCDEFGH", [34, 28, 12, 12, 60, 60, 24, 40]):
    guidance.column_dimensions[letter].width = width
guidance.protection.sheet = True

wb.save(f"{PREFIX}_pathoplexus_metadata.xlsx")

with open("versions.yml", "w") as fh:
    fh.write(f'"{TASK_PROCESS}":{NL}  python: {platform.python_version()}{NL}  openpyxl: {openpyxl.__version__}{NL}')

print(
    f"[pathoplexus] {PREFIX}: {len(included)} sample(s) written, {len(excluded)} excluded, "
    f"{'multi' if multi else 'single'}-sequence mode, raw reads {'included' if has_reads else 'not included'}"
)
