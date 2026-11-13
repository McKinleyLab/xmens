#!/usr/bin/env python
"""Record how the runs listed in <component>/record/jobs.tsv went, and fill in their blocks in RUNS.md.

For each run (one OUT folder):
  - the SLURM record of every step from sacct: state, exit code, run time, node, start and end;
  - a sha256 manifest of every file in OUT, written to <component>/record/manifests/<run>.sha256 (never overwritten:
    a manifest that no longer matches the folder is kept, and the new one is written next to it with a time stamp);
  - for a rerun (a run named <run>_rerun), every file compared with the run it repeats, by bytes and, where bytes
    differ, by content (time stamps, host names and the run's own folder path are allowed to differ, nothing else);
  - for a run given a --reference folder, the same comparison with that folder, reported but not judged.

It then replaces the block between "<!-- report <run> ..." and "<!-- end report <run> -->" in the component's
RUNS.md and writes the full detail to <component>/record/reports/<run>.<time>.md.

Usage, from the repository root:
  python record/runs_report.py mi samap --reference 2026_11_12_xmens_mi_final=/path/to/reference_run ...
  python record/runs_report.py mi --dry-run      # print the blocks; write nothing
"""

import argparse
import datetime
import gzip
import hashlib
import json
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
NOW = datetime.datetime.now()
STAMP = NOW.strftime("%Y%m%dT%H%M%S")

# Keys (JSON) and dataset names (HDF5) whose values record when, where or how long something ran; a rerun may differ
# in these and nothing else.
VOLATILE = {"written", "built", "started", "finished", "host", "hostname", "node", "date", "timestamp", "elapsed",
            "elapsed_s", "runtime", "duration", "created", "modified", "computed"}
LOG_SUFFIXES = {".log", ".out", ".err"}
TEXT_SUFFIXES = {".csv", ".tsv", ".txt", ".md", ".json", ".yaml", ".yml", ".py", ".sh"}
H5_SUFFIXES = {".h5ad", ".h5", ".hdf5"}
BLOCK = 1 << 24  # bytes read at a time when hashing or comparing


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(BLOCK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def human_size(n_bytes: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n_bytes < 1000 or unit == "TB":
            return f"{n_bytes:,.0f} {unit}" if unit == "B" else f"{n_bytes:,.1f} {unit}"
        n_bytes /= 1000


# ---------------------------------------------------------------------------------------------------------------------
# Jobs and their SLURM record

def read_jobs(component: str) -> dict:
    """Runs of a component, in file order: {run: {"out": Path, "steps": [(step, job ID), ...]}}."""
    runs = {}
    for line in (REPO / component / "record" / "jobs.tsv").read_text().splitlines():
        fields = line.split("\t")
        if not line.strip() or line.startswith("#") or fields[0] == "run":
            continue
        if len(fields) != 4:
            raise ValueError(f"{component}/record/jobs.tsv: expected run, out, step, job: {line!r}")
        run, out, step, job = fields
        runs.setdefault(run, {"out": Path(out), "steps": []})["steps"].append((step, job.strip()))
    return runs


def sacct_records(job_ids: list) -> dict:
    """One summary per job ID (an array job summarized over its tasks), from sacct."""
    command = ["sacct", "-j", ",".join(job_ids), "-X", "-n", "-P",
               "--format=JobID,JobName,State,ExitCode,Elapsed,Start,End,NodeList"]
    try:
        lines = subprocess.run(command, capture_output=True, text=True, check=True).stdout.splitlines()
    except (FileNotFoundError, subprocess.CalledProcessError) as error:
        return {job: {"error": f"sacct failed: {error}"} for job in job_ids}

    tasks = defaultdict(list)
    for line in lines:
        job_id, name, state, exit_code, elapsed, start, end, nodes = line.split("|")
        tasks[job_id.split("_")[0]].append(
            {"name": name, "state": state.split()[0], "exit": exit_code, "elapsed": elapsed, "start": start,
             "end": end, "nodes": nodes})

    records = {}
    for job in job_ids:
        rows = tasks.get(job, [])
        if not rows:
            records[job] = {"error": "not found by sacct"}
            continue
        seconds = [elapsed_seconds(r["elapsed"]) for r in rows]
        starts = [r["start"] for r in rows if r["start"] not in ("Unknown", "None")]
        ends = [r["end"] for r in rows if r["end"] not in ("Unknown", "None")]
        records[job] = {
            "name": rows[0]["name"],
            "n_tasks": len(rows),
            "states": Counter(r["state"] for r in rows),
            "exit_codes": Counter(r["exit"] for r in rows),
            "elapsed_min": min(seconds),
            "elapsed_max": max(seconds),
            "start": min(starts) if starts else "?",
            "end": max(ends) if ends else "?",
            "nodes": sorted({r["nodes"] for r in rows}),
        }
    return records


def elapsed_seconds(text: str) -> int:
    days, _, clock = text.rpartition("-")
    parts = [int(p) for p in clock.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0)
    return int(days or 0) * 86400 + parts[0] * 3600 + parts[1] * 60 + parts[2]


def clock(seconds: int) -> str:
    return f"{seconds // 3600}:{seconds % 3600 // 60:02d}:{seconds % 60:02d}"


def step_ok(record: dict) -> bool:
    return "error" not in record and set(record["states"]) == {"COMPLETED"} and set(record["exit_codes"]) == {"0:0"}


def step_line(step: str, job: str, record: dict) -> str:
    if "error" in record:
        return f"{step} {job}: {record['error']}"
    states = ", ".join(f"{n:,} {s}" if record["n_tasks"] > 1 else s for s, n in record["states"].most_common())
    if record["n_tasks"] > 1:
        timing = f"{clock(record['elapsed_min'])}-{clock(record['elapsed_max'])} per task"
        where = f"{len(record['nodes'])} nodes"
    else:
        timing = clock(record["elapsed_max"])
        where = record["nodes"][0]
    exits = ", ".join(sorted(record["exit_codes"]))
    tasks = f"{record['n_tasks']:,} tasks: " if record["n_tasks"] > 1 else ""
    return (f"{step} {job}: {tasks}{states} (exit {exits}), {timing}, {where}, "
            f"{record['start'].replace('T', ' ')[:16]} to {record['end'].replace('T', ' ')[:16]}")


# ---------------------------------------------------------------------------------------------------------------------
# Manifests

def build_manifest(out: Path) -> dict:
    """{relative path: (sha256, bytes)} for every regular file under out, in sorted order."""
    manifest = {}
    for path in sorted(p for p in out.rglob("*") if p.is_file() and not p.is_symlink()):
        manifest[path.relative_to(out).as_posix()] = (sha256_file(path), path.stat().st_size)
    return manifest


def manifest_text(manifest: dict) -> str:
    return "".join(f"{sha}  {rel}\n" for rel, (sha, _) in manifest.items())


def save_manifest(component: str, run: str, manifest: dict, dry_run: bool) -> tuple:
    """Write the manifest unless an identical one exists; return (path, its sha256, note)."""
    folder = REPO / component / "record" / "manifests"
    text = manifest_text(manifest)
    digest = hashlib.sha256(text.encode()).hexdigest()
    path = folder / f"{run}.sha256"
    note = ""
    if path.exists():
        if path.read_text() == text:
            return path, digest, "unchanged since the last report"
        path = folder / f"{run}.{STAMP}.sha256"
        note = f"the folder changed since {run}.sha256 was written; this manifest is the current state"
    if not dry_run:
        folder.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
    return path, digest, note


# ---------------------------------------------------------------------------------------------------------------------
# Comparing two run folders file by file

class Normalizer:
    """Replaces each run's own folder path and name with one token, so files that record where they were written can
    still be compared."""

    def __init__(self, *folders: Path):
        self.pairs = []
        for folder in folders:
            self.pairs += [(str(folder), "<RUN>"), (folder.name, "<RUN>")]
        self.pairs.sort(key=lambda pair: -len(pair[0]))

    def __call__(self, text: str) -> str:
        for old, new in self.pairs:
            text = text.replace(old, new)
        return text


def compare_runs(a: Path, manifest_a: dict, b: Path, manifest_b: dict) -> list:
    """[(relative path, category, detail)] for the union of the two folders' files."""
    normalize = Normalizer(a, b)
    rows = []
    for rel in sorted(set(manifest_a) | set(manifest_b)):
        if is_log(rel) and (rel not in manifest_a or rel not in manifest_b):
            rows.append((rel, "log", "log in one folder only"))
        elif rel not in manifest_b:
            rows.append((rel, "only in run", ""))
        elif rel not in manifest_a:
            rows.append((rel, "only in other", ""))
        elif manifest_a[rel][0] == manifest_b[rel][0]:
            rows.append((rel, "identical", ""))
        elif is_log(rel):
            rows.append((rel, "log", "logs record times and hosts; not compared"))
        else:
            try:
                category, detail = compare_content(a / rel, b / rel, normalize)
            except Exception as error:  # a comparison that cannot run is reported, never taken as a match
                category, detail = "differs", f"could not compare: {type(error).__name__}: {error}"
            rows.append((rel, category, detail))
    return rows


def is_log(rel: str) -> bool:
    return rel.startswith("logs/") or "/logs/" in rel or Path(rel).suffix in LOG_SUFFIXES


def compare_content(a: Path, b: Path, normalize: Normalizer) -> tuple:
    suffixes = "".join(a.suffixes[-2:])
    suffix = a.suffix
    if suffixes.endswith(".gz"):
        return compare_bytes_normalized(gzip.decompress(a.read_bytes()), gzip.decompress(b.read_bytes()),
                                        normalize, "gzip (header time stamp)")
    if suffix == ".json":
        return compare_json(a, b, normalize)
    if suffix in TEXT_SUFFIXES:
        return compare_bytes_normalized(a.read_bytes(), b.read_bytes(), normalize, "text")
    if suffix in H5_SUFFIXES:
        return compare_h5(a, b, normalize)
    if suffix in (".npz", ".npy"):
        return compare_numpy(a, b)
    if suffix == ".xlsx":
        return compare_xlsx(a, b, normalize)
    if suffix == ".ipynb":
        return compare_notebook(a, b, normalize)
    if suffix == ".svg":
        strip = lambda t: re.sub(r"<dc:date>.*?</dc:date>", "", t, flags=re.S)
        return compare_text(strip(normalize(a.read_text())), strip(normalize(b.read_text())), "svg (date)")
    if suffix == ".pdf":
        strip = lambda d: re.sub(rb"/(CreationDate|ModDate) \([^)]*\)|/ID \[[^\]]*\]", b"", d)
        same = strip(a.read_bytes()) == strip(b.read_bytes())
        return ("same content", "pdf: same apart from its dates and ID") if same else ("differs", "pdf bytes")
    if suffix == ".png":
        from PIL import Image
        with Image.open(a) as image_a, Image.open(b) as image_b:
            same = image_a.size == image_b.size and np.array_equal(np.asarray(image_a), np.asarray(image_b))
        return ("same content", "png: identical pixels") if same else ("differs", "png pixels")
    if suffix == ".parquet":
        import pandas as pd
        same = pd.read_parquet(a).equals(pd.read_parquet(b))
        return ("same content", "parquet: identical table") if same else ("differs", "parquet table")
    if suffix == ".pkl":
        return "not compared", "pickle; not loaded (its contents are compared through the h5ad outputs)"
    return "differs", "binary; no content comparison for this file type"


def compare_bytes_normalized(data_a: bytes, data_b: bytes, normalize: Normalizer, kind: str) -> tuple:
    text_a = normalize(data_a.decode("utf-8", errors="replace"))
    text_b = normalize(data_b.decode("utf-8", errors="replace"))
    return compare_text(text_a, text_b, kind)


def compare_text(text_a: str, text_b: str, kind: str) -> tuple:
    if text_a == text_b:
        return "same content", f"{kind}: identical once each run's folder path is replaced by <RUN>"
    lines_a, lines_b = text_a.splitlines(), text_b.splitlines()
    n_diff = sum(x != y for x, y in zip(lines_a, lines_b)) + abs(len(lines_a) - len(lines_b))
    first = next((i for i, (x, y) in enumerate(zip(lines_a, lines_b)) if x != y), min(len(lines_a), len(lines_b)))
    n_lines = max(len(lines_a), len(lines_b))
    return "differs", f"{kind}: {n_diff:,} of {n_lines:,} lines differ, first at line {first + 1}"


def flatten(value, prefix=""):
    if isinstance(value, dict):
        for key, item in value.items():
            yield from flatten(item, f"{prefix}{key}/")
    else:
        yield prefix.rstrip("/"), value


def compare_json(a: Path, b: Path, normalize: Normalizer) -> tuple:
    flat_a = dict(flatten(json.loads(normalize(a.read_text()))))
    flat_b = dict(flatten(json.loads(normalize(b.read_text()))))
    differing = sorted(k for k in set(flat_a) | set(flat_b) if flat_a.get(k, "<missing>") != flat_b.get(k, "<missing>"))
    if not differing:
        return "same content", "json: identical once each run's folder path is replaced by <RUN>"
    volatile = [k for k in differing if k.split("/")[-1].lower() in VOLATILE]
    if len(volatile) == len(differing):
        return "same content", "json: differs only in " + ", ".join(differing)
    shown = [f"{k}: {short(flat_a.get(k, '<missing>'))} vs {short(flat_b.get(k, '<missing>'))}" for k in differing]
    return "differs", "json: " + "; ".join(shown[:8]) + (f"; ... {len(shown) - 8} more" if len(shown) > 8 else "")


def short(value, width: int = 48) -> str:
    """A value cut to width, keeping both ends (paths and checksums differ at either end)."""
    text = str(value)
    return text if len(text) <= width else f"{text[:18]}...{text[-(width - 21):]}"


def compare_numpy(a: Path, b: Path) -> tuple:
    def load(path):
        if path.suffix == ".npz":
            with np.load(path, allow_pickle=False) as archive:
                return dict(archive)
        return {"": np.load(path, allow_pickle=False)}

    arrays_a, arrays_b = load(a), load(b)
    differing = [k for k in sorted(set(arrays_a) | set(arrays_b))
                 if k not in arrays_a or k not in arrays_b or not arrays_equal(arrays_a[k], arrays_b[k])]
    if not differing:
        return "same content", "numpy: identical arrays"
    return "differs", "numpy arrays differ: " + ", ".join(k or "(array)" for k in differing)


def arrays_equal(x: np.ndarray, y: np.ndarray) -> bool:
    if x.shape != y.shape or x.dtype != y.dtype:
        return False
    if x.dtype.kind in "fc":
        return bool(np.array_equal(x, y, equal_nan=True))
    return bool(np.array_equal(x, y))


def compare_h5(a: Path, b: Path, normalize: Normalizer) -> tuple:
    """Every dataset's values and every attribute, in both files; strings compared after normalization."""
    import h5py

    def entries(handle):
        found = {}
        handle.visititems(lambda name, obj: found.__setitem__(name, obj))
        return found

    differing, missing = [], []
    with h5py.File(a, "r") as file_a, h5py.File(b, "r") as file_b:
        objects_a, objects_b = entries(file_a), entries(file_b)
        missing = sorted(set(objects_a) ^ set(objects_b))
        for name in sorted(set(objects_a) & set(objects_b)):
            obj_a, obj_b = objects_a[name], objects_b[name]
            if not attrs_equal(obj_a.attrs, obj_b.attrs, normalize):
                differing.append(name + " (attributes)")
            if isinstance(obj_a, h5py.Dataset) and not datasets_equal(obj_a, obj_b, normalize):
                differing.append(name)
        if not attrs_equal(file_a.attrs, file_b.attrs, normalize):
            differing.append("/ (attributes)")
    if not differing and not missing:
        return "same content", "hdf5: every dataset and attribute identical"
    problems = missing + differing
    if not missing and all(p.split(" ")[0].split("/")[-1].lower() in VOLATILE for p in differing):
        return "same content", "hdf5: differs only in " + ", ".join(differing)
    shown = ", ".join(problems[:10]) + (f", ... {len(problems) - 10} more" if len(problems) > 10 else "")
    return "differs", f"hdf5: {len(missing)} objects in one file only, {len(differing)} differ: {shown}"


def attrs_equal(attrs_a, attrs_b, normalize: Normalizer) -> bool:
    if set(attrs_a) != set(attrs_b):
        return False
    return all(values_equal(np.asarray(attrs_a[k]), np.asarray(attrs_b[k]), normalize) for k in attrs_a)


def datasets_equal(dataset_a, dataset_b, normalize: Normalizer) -> bool:
    if dataset_a.shape != dataset_b.shape or dataset_a.dtype != dataset_b.dtype:
        return False
    if dataset_a.shape == ():
        return values_equal(np.asarray(dataset_a[()]), np.asarray(dataset_b[()]), normalize)
    rows = max(1, 50_000_000 // max(1, int(np.prod(dataset_a.shape[1:]))))  # about 50 million values at a time
    for start in range(0, dataset_a.shape[0], rows):
        if not values_equal(dataset_a[start:start + rows], dataset_b[start:start + rows], normalize):
            return False
    return True


def values_equal(x, y, normalize: Normalizer) -> bool:
    x, y = np.asarray(x), np.asarray(y)
    if x.dtype.kind in "OSU" or y.dtype.kind in "OSU":
        as_text = lambda v: normalize(v.decode("utf-8", "replace") if isinstance(v, bytes) else str(v))
        if x.shape != y.shape:
            return False
        return all(as_text(u) == as_text(v) for u, v in zip(x.ravel().tolist(), y.ravel().tolist()))
    return arrays_equal(x, y)


def compare_xlsx(a: Path, b: Path, normalize: Normalizer) -> tuple:
    from openpyxl import load_workbook

    def cells(path):
        book = load_workbook(path, read_only=True, data_only=True)
        sheets = {name: [[normalize(str(v)) if isinstance(v, str) else v for v in row]
                         for row in book[name].iter_rows(values_only=True)] for name in book.sheetnames}
        book.close()
        return sheets

    sheets_a, sheets_b = cells(a), cells(b)
    if sheets_a == sheets_b:
        n_cells = sum(len(row) for rows in sheets_a.values() for row in rows)
        return "same content", f"xlsx: {len(sheets_a)} sheets, {n_cells:,} cells identical (file metadata differs)"
    differing = [name for name in sorted(set(sheets_a) | set(sheets_b)) if sheets_a.get(name) != sheets_b.get(name)]
    return "differs", "xlsx sheets differ: " + ", ".join(differing)


def compare_notebook(a: Path, b: Path, normalize: Normalizer) -> tuple:
    code = lambda p: [normalize("".join(c["source"])) for c in json.loads(p.read_text())["cells"]
                      if c["cell_type"] == "code"]
    code_a, code_b = code(a), code(b)
    if code_a == code_b:
        return "same content", f"notebook: same {len(code_a)} code cells (outputs and run times not compared)"
    n_diff = sum(x != y for x, y in zip(code_a, code_b)) + abs(len(code_a) - len(code_b))
    return "differs", f"notebook code: {n_diff} of {max(len(code_a), len(code_b))} code cells differ"


def summarize(rows: list) -> dict:
    counts = Counter(category for _, category, _ in rows)
    by_folder = defaultdict(Counter)
    for rel, category, _ in rows:
        if category not in ("identical",):
            folder = "/".join(rel.split("/")[:-1][:2]) or "."
            by_folder[folder][category] += 1
    return {"counts": counts, "by_folder": by_folder}


def comparison_sentence(rows: list) -> str:
    counts = summarize(rows)["counts"]
    parts = [f"{len(rows):,} files: {counts['identical']:,} identical bytes"]
    if counts["same content"]:
        kinds = Counter(Path(rel).suffix or rel for rel, category, _ in rows if category == "same content")
        parts.append(f"{counts['same content']:,} same content ("
                     + ", ".join(f"{n} {s}" for s, n in kinds.most_common()) + "; time stamps or paths only)")
    for category in ("log", "not compared"):
        if counts[category]:
            parts.append(f"{counts[category]:,} {'logs' if category == 'log' else 'not compared'}")
    parts.append(f"{counts['differs']:,} differ")
    missing = counts["only in run"] + counts["only in other"]
    parts.append(f"{missing:,} in one folder only")
    return "; ".join(parts)


def rerun_passes(rows: list) -> bool:
    counts = summarize(rows)["counts"]
    return counts["differs"] == 0 and counts["only in run"] == 0 and counts["only in other"] == 0


# ---------------------------------------------------------------------------------------------------------------------
# Numbers worth seeing next to each other

def key_numbers(component: str, out: Path) -> dict:
    """A few headline numbers of a run, for the report (not judged)."""
    numbers = {}
    if component == "mi":
        path = out / "overlap" / "methods_overlap.json"
        if path.exists():
            data = dict(flatten(json.loads(path.read_text())))
            numbers = {k: v for k, v in data.items() if k.split("/")[-1] in ("N", "K", "n", "k", "log10_p")}
    elif component == "samap":
        path = out / "2_samap" / "meta.json"
        if path.exists():
            meta = json.loads(path.read_text())
            numbers = {k: meta.get(k) for k in ("mean_alignment", "seed", "hnswlib_threads") if k in meta}
    return numbers


# ---------------------------------------------------------------------------------------------------------------------
# RUNS.md

def entry_code_commit(runs_md: str, run: str):
    """The commit named on the Code line of the entry whose block belongs to run, if any."""
    start = runs_md.find(f"<!-- report {run} ")
    if start < 0:
        return None
    heading = max(runs_md.rfind("\n## ", 0, start), 0)
    match = re.search(r"^Code:.*?@ ([0-9a-f]{7,40})", runs_md[heading:start], flags=re.M)
    return match.group(1) if match else None


def code_changes_since(commit: str, component: str) -> list:
    """Files of the component (records and logs aside) that differ now from the commit the run names."""
    command = ["git", "-C", str(REPO), "diff", "--name-only", commit, "--", component,
               f":(exclude){component}/record", f":(exclude){component}/logs"]
    try:
        return subprocess.run(command, capture_output=True, text=True, check=True).stdout.split()
    except (FileNotFoundError, subprocess.CalledProcessError) as error:
        return [f"(git diff failed: {error})"]


def replace_block(runs_md: str, run: str, block: str) -> str:
    pattern = re.compile(rf"<!-- report {re.escape(run)} .*?-->\n.*?<!-- end report {re.escape(run)} -->\n?", re.S)
    if not pattern.search(runs_md):  # no block to fill: add the report at the end rather than lose it
        return runs_md.rstrip("\n") + f"\n\n## Report for {run} (RUNS.md had no block for this run)\n" + block
    return pattern.sub(lambda _: block, runs_md, count=1)


def wrap(label, text: str, width: int = 120) -> str:
    """'Label:    text' wrapped at width, continuation lines indented under the text; label None: indent only."""
    indent = " " * 10
    lines, line = [], indent if label is None else f"{label + ':':<10}"
    for word in text.split(" "):
        if len(line) + len(word) > width and line.strip() not in ("", f"{label}:"):
            lines.append(line.rstrip())
            line = indent
        line += word + " "
    lines.append(line.rstrip())
    return "\n".join(lines)


# ---------------------------------------------------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("components", nargs="+", help="pipeline folders with record/jobs.tsv, e.g. mi samap")
    parser.add_argument("--reference", action="append", default=[], metavar="RUN=FOLDER",
                        help="also compare RUN with FOLDER (reported, not judged)")
    parser.add_argument("--dry-run", action="store_true", help="print the blocks; write nothing")
    args = parser.parse_args()
    references = dict(item.split("=", 1) for item in args.reference)

    for component in args.components:
        runs = read_jobs(component)
        all_jobs = [job for run in runs.values() for _, job in run["steps"]]
        records = sacct_records(all_jobs)
        runs_md_path = REPO / component / "record" / "RUNS.md"
        runs_md = runs_md_path.read_text()

        manifests, manifest_notes = {}, {}
        for run, info in runs.items():
            print(f"{component} {run}: hashing {info['out']}", flush=True)
            manifests[run] = build_manifest(info["out"]) if info["out"].is_dir() else {}
            manifest_notes[run] = save_manifest(component, run, manifests[run], args.dry_run)

        comparisons = {}  # run -> (other folder, rows, judged)
        for run, info in runs.items():
            original = run.removesuffix("_rerun")
            if run.endswith("_rerun") and original in runs:
                print(f"{component} {run}: comparing with {original}", flush=True)
                comparisons[run] = (runs[original]["out"], compare_runs(
                    info["out"], manifests[run], runs[original]["out"], manifests[original]), True)
            elif run in references:
                folder = Path(references[run])
                print(f"{component} {run}: comparing with the reference {folder}", flush=True)
                reference_manifest = build_manifest(folder)
                save_manifest(component, folder.name, reference_manifest, args.dry_run)
                comparisons[run] = (folder, compare_runs(info["out"], manifests[run], folder, reference_manifest),
                                    False)

        for run, info in runs.items():
            block, detail = render(component, run, info, runs, records, manifests, manifest_notes, comparisons,
                                   references, runs_md)
            report_path = REPO / component / "record" / "reports" / f"{run}.{STAMP}.md"
            block = block.replace("Report:   ____", wrap("Report", str(report_path.relative_to(REPO))))
            if args.dry_run:
                print(block)
                continue
            report_path.parent.mkdir(parents=True, exist_ok=True)
            report_path.write_text(detail)
            runs_md = replace_block(runs_md, run, block)
        if not args.dry_run:
            runs_md_path.write_text(runs_md)
            print(f"wrote {runs_md_path.relative_to(REPO)} and {len(runs)} reports", flush=True)
    return 0


def render(component, run, info, runs, records, manifests, manifest_notes, comparisons, references, runs_md):
    """The RUNS.md block of a run and its detailed report."""
    steps = info["steps"]
    lines_steps = [step_line(step, job, records[job]) for step, job in steps]
    all_ok = all(step_ok(records[job]) for _, job in steps)
    manifest = manifests[run]
    manifest_path, manifest_sha, manifest_note = manifest_notes[run]
    total = sum(size for _, size in manifest.values())

    block = [f"<!-- report {run} (record/runs_report.py fills in this block; written {NOW:%Y-%m-%d %H:%M}) -->"]
    block.append(wrap("Ran", lines_steps[0]))
    block += [wrap(None, line) for line in lines_steps[1:]]
    outputs = (f"{len(manifest):,} files, {human_size(total)}; sha256 of each in "
               f"{manifest_path.relative_to(REPO)} (sha256 {manifest_sha[:16]}...)")
    block.append(wrap("Outputs", outputs + (f"; {manifest_note}" if manifest_note else "")))

    status_parts = ["complete" if all_ok else "INCOMPLETE OR FAILED: " + ", ".join(
        f"{step} {job}" for step, job in steps if not step_ok(records[job]))]
    if run in comparisons:
        folder, rows, judged = comparisons[run]
        label = "Checks" if judged else "Compared"
        what = f"with {folder.name}" + ("" if judged else " (reference; differences expected and not judged)")
        block.append(wrap(label, f"{what}, file by file: {comparison_sentence(rows)}"))
        if judged:
            status_parts.append("identical to the run of record in every output" if rerun_passes(rows)
                                else "DIFFERS from the run of record (see the report)")
    rerun = f"{run}_rerun"
    if rerun in comparisons:
        rows = comparisons[rerun][1]
        block.append(wrap("Rerun", f"{rerun}: " + ("every output reproduced" if rerun_passes(rows)
                                                   else "outputs DIFFER (see that entry)")))
        status_parts.append("reproduced by its rerun" if rerun_passes(rows) else "NOT reproduced by its rerun")

    numbers = key_numbers(component, info["out"])
    if numbers:
        block.append(wrap("Numbers", "; ".join(f"{k} {v}" for k, v in numbers.items())))
    commit = entry_code_commit(runs_md, run)
    if commit:
        changed = code_changes_since(commit, component)
        block.append(wrap("Code now", f"{component}/ as at {commit}" if not changed else
                          f"differs from {commit} in {', '.join(changed)} (records and logs aside)"))
    block.append("Report:   ____")
    block.append(wrap("Status", "; ".join(status_parts)))
    block.append(f"<!-- end report {run} -->")
    block_text = "\n".join(block) + "\n"

    detail = [f"# {component}: {run}", "", f"Written {NOW:%Y-%m-%d %H:%M:%S} by record/runs_report.py", "",
              f"OUT: {info['out']}", "", "## Steps (sacct)", ""]
    detail += [f"- {line}" for line in lines_steps]
    for step, job in steps:
        record = records[job]
        if "error" not in record and record["n_tasks"] > 1:
            detail.append(f"  - {step}: states {dict(record['states'])}, exit codes {dict(record['exit_codes'])}, "
                          f"nodes {', '.join(record['nodes'][:20])}{' ...' if len(record['nodes']) > 20 else ''}")
    detail += ["", "## Outputs", "", outputs, ""]
    if numbers:
        detail += ["## Numbers", ""] + [f"- {k}: {v}" for k, v in numbers.items()] + [""]
    if run in comparisons:
        folder, rows, judged = comparisons[run]
        detail += [f"## Compared with {folder}", "", comparison_sentence(rows), ""]
        groups = summarize(rows)["by_folder"]
        if groups:
            detail += ["Files not identical in bytes, by folder:", ""]
            detail += [f"- {folder_name}: " + ", ".join(f"{n} {c}" for c, n in counts.most_common())
                       for folder_name, counts in sorted(groups.items())]
            detail.append("")
        detail += ["| file | result | detail |", "|---|---|---|"]
        detail += [f"| {rel} | {category} | {text} |" for rel, category, text in rows if category != "identical"]
        detail.append("")
    detail += ["## Manifest", "", f"{manifest_path.relative_to(REPO)} (sha256 {manifest_sha})", ""]
    return block_text, "\n".join(detail)


if __name__ == "__main__":
    sys.exit(main())
