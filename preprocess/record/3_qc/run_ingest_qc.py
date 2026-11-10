"""Run the pipeline's own ingest -> QC -> concat scripts outside Snakemake, on either set of CellBender files.

The scripts are the ones rules ingest (03), qc (04), qc_report (05) and concat (06) run, executed unchanged: each is
handed the `snakemake` object its rule would give it (input, output, params, wildcards, config, log).

  --inputs andy    samples.tsv's `path` (xmens_inputs, Andy's CellBender files): the pipeline's own inputs, so the
                   outputs can be checked against data/interim (compare_qc_outputs.py) to show that the scripts
                   reproduce the pipeline's objects;
  --inputs regen   xmens_inputs_regen/<sample>/<sample>_cellbender_filtered.h5 (the February 2026 re-run's
                   cells-only CellBender output).

Everything goes under --out, never into the repo, so the pipeline's data/interim objects are untouched.

QC runs in the DoubletDetection environment named in config.yaml (envs.doubletdetection), as rule qc does; ingest,
qc_report and concat run in the Python that runs this file (rsc_uv), as in the pipeline. DoubletDetection is random,
so its doublet calls differ slightly from any earlier run's.

Outputs are never overwritten. A step whose output exists is skipped, so an interrupted run can simply be started
again; a step's .h5ad is written under .incomplete/ and moved into place only when the step succeeds, so an
interrupted step never leaves a partial file that a later run would take as finished.

Usage (rsc_uv, on a compute node; run_ingest_qc.sbatch runs one sample per array task, then concat and the check):
  python run_ingest_qc.py --inputs andy  --out /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_05_qc_rerun_andy
  python run_ingest_qc.py --inputs regen --out /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_05_qc_rerun_regen
  python run_ingest_qc.py --inputs andy --out /n/netscratch/eddy_lab/Everyone/nhilgert/2026_11_05_qc_one_sample \\
      --samples 1_dpi_early_mck1 --steps ingest qc

Outputs under --out, mirroring the pipeline's layout:
  raw/<sample>.h5ad          ingest (all CellBender cells, coords attached)
  qc/<sample>.h5ad           QC (threshold QC + DoubletDetection)
  qc_passed_<species>.h5ad   concat (the counterpart of data/interim/qc_passed_<species>.h5ad)
  qc_figures/<sample>/       QC report figures (only with --steps report)
  logs/                      the scripts' own logs, same names as in the pipeline
  manifests/                 one record per invocation: inputs (with sha256), steps, repo commit, Pythons, node, times
"""
import argparse
import hashlib
import json
import os
import socket
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import yaml

REPO = Path("/n/lab_storage/eddy_lab/users/nhilgert/xmens_final")
SCRIPTS = REPO / "workflow/scripts/ingest_qc"
REGEN_DIR = Path("/n/lab_storage/eddy_lab/users/nhilgert/decidualization/xmens_inputs_regen")
CACHE = Path("/n/lab_storage/eddy_lab/users/nhilgert/.cache_scratch")   # keeps caches out of $HOME

# Runs inside the target Python: rebuild the `snakemake` object from a JSON spec and execute the script unchanged.
SHIM = r"""
import json, runpy, sys
from types import SimpleNamespace

class NamedList(list):
    # Snakemake's input/output/log are lists whose items are also attributes (snakemake.input.h5ad).
    def __init__(self, items=(), named=None):
        named = named or {}
        super().__init__(list(items) + list(named.values()))
        for key, value in named.items():
            setattr(self, key, value)

spec = json.load(open(sys.argv[1]))
snakemake = SimpleNamespace(
    input=NamedList(spec["input"].get("items", []), spec["input"].get("named")),
    output=NamedList(named=spec["output"]),
    params=SimpleNamespace(**spec["params"]),
    wildcards=SimpleNamespace(**spec["wildcards"]),
    config=spec["config"],
    log=NamedList([spec["log"]]),
)
runpy.run_path(spec["script"], run_name="__main__", init_globals={"snakemake": snakemake})
"""


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def run_script(python, script, spec, out, label):
    """Execute one pipeline script with its snakemake object; stop with the log tail if it fails."""
    key, target = next(iter(spec["output"].items()))
    target = Path(target)
    if target.exists():
        print(f"  {label}: output exists, skipped ({target})", flush=True)
        return
    log = Path(spec["log"])
    if log.exists():   # left by an earlier attempt that did not finish: keep it under another name
        stamp = time.strftime("%Y%m%d_%H%M%S", time.localtime(log.stat().st_mtime))
        log.rename(log.with_name(f"{log.stem}.unfinished_{stamp}{log.suffix}"))
    written = out / ".incomplete" / f"{os.getpid()}_{target.name}" if target.suffix == ".h5ad" else target
    spec = dict(spec, script=str(script), output={key: str(written)})
    spec_path = out / "specs" / f"{label}_{time.strftime('%Y%m%d_%H%M%S')}.json"
    spec_path.write_text(json.dumps(spec, indent=1, default=str))
    env = dict(os.environ, HDF5_USE_FILE_LOCKING="FALSE", MPLCONFIGDIR=str(CACHE / "mpl"),
               NUMBA_CACHE_DIR=str(CACHE / "numba"), XDG_CACHE_HOME=str(CACHE / "xdg"))
    start = time.time()
    result = subprocess.run([python, "-c", SHIM, str(spec_path)], cwd=out, env=env)
    lines = log.read_text().splitlines() if log.exists() else []
    if result.returncode != 0:
        sys.exit(f"{label} failed (exit {result.returncode}); last log lines:\n" + "\n".join(lines[-15:]))
    if written != target:
        os.replace(written, target)
    print(f"  {label}: done in {time.time() - start:.0f} s", flush=True)
    for line in lines:
        if "WARNING" in line:   # e.g. DoubletDetection missing or failed, which qc.py survives
            print(f"  {label}: {line}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--inputs", required=True, choices=["andy", "regen"],
                        help="andy: samples.tsv paths (the pipeline's inputs); regen: the February 2026 re-run")
    parser.add_argument("--out", required=True, type=Path, help="new output folder (outside the repo)")
    parser.add_argument("--samples", nargs="*", help="sample IDs to run (default: all in samples.tsv)")
    parser.add_argument("--steps", nargs="*", default=["ingest", "qc", "concat"],
                        choices=["ingest", "qc", "report", "concat"], help="steps to run, in pipeline order")
    args = parser.parse_args()

    out = args.out.resolve()
    if REPO in out.parents or out == REPO:
        sys.exit("--out must be outside the repo, so the pipeline's own objects are never touched")
    for sub in ("raw", "qc", "logs", "specs", "qc_figures", "manifests", ".incomplete"):
        (out / sub).mkdir(parents=True, exist_ok=True)

    config = yaml.safe_load((REPO / "config/config.yaml").read_text())
    samples = pd.read_csv(REPO / "config/samples.tsv", sep="\t").set_index("sample_id", drop=False)
    if args.samples:
        unknown = set(args.samples) - set(samples.index)
        if unknown:
            sys.exit(f"unknown samples: {sorted(unknown)}")
        samples = samples.loc[[s for s in samples.index if s in args.samples]]
    if args.inputs == "andy":
        counts = {s: Path(samples.loc[s, "path"]) for s in samples.index}
    else:
        counts = {s: REGEN_DIR / s / f"{s}_cellbender_filtered.h5" for s in samples.index}
    missing = [str(p) for p in counts.values() if not p.exists()]
    if missing:
        sys.exit(f"missing CellBender files: {missing}")

    main_python = sys.executable
    doublet_python = str(Path(config["envs"]["doubletdetection"]) / "bin/python")
    coords = config["coords"]

    def git(*command):
        return subprocess.run(["git", "-C", str(REPO), *command], capture_output=True, text=True).stdout.strip()

    manifest = {
        "started": time.strftime("%Y-%m-%d %H:%M:%S"), "inputs_set": args.inputs, "out": str(out),
        "samples": list(samples.index), "steps": args.steps,
        "repo_commit": git("rev-parse", "HEAD"),
        "modified_since_commit": git("status", "--porcelain", "--", str(SCRIPTS), "config/config.yaml",
                                     "config/samples.tsv") or "none",
        "python_ingest_concat": main_python, "python_qc": doublet_python,
        "host": socket.gethostname(), "cpus_on_node": os.cpu_count(),
        "slurm": {k: os.environ[k] for k in ("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID",
                                             "SLURM_CPUS_PER_TASK") if k in os.environ},
        "inputs": {s: {"path": str(p), "bytes": p.stat().st_size, "sha256": sha256(p),
                       "modified": time.strftime("%Y-%m-%d %H:%M", time.localtime(p.stat().st_mtime))}
                   for s, p in counts.items()},
    }

    for sample, row in samples.iterrows():
        print(f"{sample}", flush=True)
        common = {"wildcards": {"sample": sample}, "config": config}
        if "ingest" in args.steps:
            run_script(main_python, SCRIPTS / "ingest.py", dict(
                common, input={"named": {"counts": str(counts[sample])}},
                output={"h5ad": str(out / "raw" / f"{sample}.h5ad")},
                params={"species": row.species, "timepoint": row.timepoint, "modality": row.modality,
                        "coords_id": row.coords_id, "coords_dir": coords["dir"], "coords_pattern": coords["pattern"],
                        "coords_x": coords["x_col"], "coords_y": coords["y_col"]},
                log=str(out / "logs" / f"00b_ingest_{sample}.log")), out, f"ingest_{sample}")
        if "qc" in args.steps:
            run_script(doublet_python, SCRIPTS / "qc.py", dict(
                common, input={"named": {"h5ad": str(out / "raw" / f"{sample}.h5ad")}},
                output={"h5ad": str(out / "qc" / f"{sample}.h5ad")}, params={"species": row.species},
                log=str(out / "logs" / f"01_qc_{sample}.log")), out, f"qc_{sample}")
        if "report" in args.steps:
            figures = out / "qc_figures" / sample
            figures.mkdir(exist_ok=True)
            run_script(main_python, SCRIPTS / "qc_report.py", dict(
                common, input={"named": {"raw": str(out / "raw" / f"{sample}.h5ad"),
                                         "qcd": str(out / "qc" / f"{sample}.h5ad")}},
                output={"done": str(figures / ".done")}, params={"species": row.species, "outdir": str(figures)},
                log=str(out / "logs" / f"01b_qc_report_{sample}.log")), out, f"report_{sample}")
            (figures / ".done").touch()

    if "concat" in args.steps:
        all_samples = pd.read_csv(REPO / "config/samples.tsv", sep="\t")
        for species, group in samples.groupby("species", sort=False):
            # A partial per-species object would carry the full object's name and block the real one later.
            expected = set(all_samples.loc[all_samples["species"] == species, "sample_id"])
            if set(group.index) != expected:
                print(f"  concat_{species}: skipped (needs all {len(expected)} {species} samples)", flush=True)
                continue
            run_script(main_python, SCRIPTS / "concat.py", {
                "input": {"items": [str(out / "qc" / f"{s}.h5ad") for s in group.index]},
                "output": {"h5ad": str(out / f"qc_passed_{species}.h5ad")}, "params": {"species": species},
                "wildcards": {"species": species}, "config": config,
                "log": str(out / "logs" / f"01c_concat_{species}.log")}, out, f"concat_{species}")

    manifest["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
    label = samples.index[0] if len(samples) == 1 else f"{len(samples)}_samples"
    manifest_path = out / "manifests" / f"{time.strftime('%Y%m%d_%H%M%S')}_{label}_{os.getpid()}.json"
    manifest_path.write_text(json.dumps(manifest, indent=1, default=str))
    print(f"done; manifest {manifest_path}")


if __name__ == "__main__":
    main()
