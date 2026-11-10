"""Replace the obs of final_objects/Samap_Final.h5ad with the six final columns, leaving every other dataset untouched.

New obs (same nuclei, same order): species; sample, timepoint, cell_type from the final fib/dec objects (current names);
transferred_mouse_pseudotime = transferred_mm_dpt_pseudotime_ranknorm (the transfer of record, human rows only);
transferred_mouse_fibdec = transferred_mm_coarse (mouse cell type by SAMap edge-weighted vote; 'unmapped' and the
mouse rows' 'NA' become missing). The obs group is rewritten in a copy; the copy replaces Samap_Final.h5ad only after
every non-obs dataset is shown byte-identical to the pinned original and every new column equals its source.
"""
import datetime
import hashlib
import json
import os
import shutil
import sys

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, write_elem

FO = "/n/netscratch/eddy_lab/Everyone/nhilgert/final_objects"
TARGET = f"{FO}/Samap_Final.h5ad"
TEMP = f"{FO}/.Samap_Final.slim.partial.h5ad"
RECORD = "/n/lab_storage/eddy_lab/users/nhilgert/code/samap-scvi/full_fibdec/runs_raw/vanilla_raw/samap_transferred_pt_harmonized.h5ad"
PIN = "62fb75d7c2ca3e739299be48695476b6c7a8ac48a05c9226c3c0cb47b8a3ad84"
OBJECTS = {"mm": f"{FO}/Cevrim_XMens_Only_Fib_Dec_Cells.h5ad", "hs": f"{FO}/Cevrim_Human_Only_Fib_Dec_Cells.h5ad"}
CELL_TYPE_CHECK = ("/n/lab_storage/eddy_lab/users/nhilgert/tmp/2026_11_02_samap_final_figures/repo/results/figures/"
                   "samap_main/final/transferred_cell_type.csv")
REPORT = "/n/lab_storage/eddy_lab/users/nhilgert/tmp/2026_11_02_samap_final_figures/samap_final_slim_report.json"


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def fail(message):
    if os.path.exists(TEMP):
        os.remove(TEMP)  # our own partial copy, never the target
    sys.exit(f"STOPPED, Samap_Final.h5ad unchanged: {message}")


if os.path.exists(TEMP):
    sys.exit(f"{TEMP} exists; refusing to reuse it")
for path in (TARGET, RECORD):
    observed = sha256_of(path)
    print(f"{path}: sha256 {observed}")
    if observed != PIN:
        sys.exit(f"{path} does not have the pinned sha256; stopping before any change")

# ---- the new obs
with h5py.File(TARGET, "r") as f:
    old = read_elem(f["obs"])
names = old.index.astype(str)
species = old["species"].astype(str)
if species.value_counts().to_dict() != {"hs": 47533, "mm": 43197}:
    fail(f"unexpected species counts {species.value_counts().to_dict()}")

final = {}
for code, path in OBJECTS.items():
    with h5py.File(path, "r") as f:
        obs = read_elem(f["obs"])[["sample", "timepoint", "cell_type"]]
    rows = names[species.to_numpy() == code]
    if not rows.isin(obs.index).all() or len(obs) != len(rows):
        fail(f"{code}: the final object and the SAMap object hold different nuclei")
    final[code] = obs


def combined(column):
    """One categorical over both species, each species' values and category order from its final object."""
    categories = list(dict.fromkeys(list(final["hs"][column].cat.categories) + list(final["mm"][column].cat.categories)))
    values = pd.concat([final["hs"][column].astype(str), final["mm"][column].astype(str)]).reindex(names)
    return pd.Categorical(values.to_numpy(), categories=categories)


cell_types = combined("cell_type").categories
fibdec = old["transferred_mm_coarse"].astype(str).where(lambda s: s.isin(cell_types))
new = pd.DataFrame({
    "species": pd.Categorical(species.to_numpy(), categories=["hs", "mm"]),
    "sample": combined("sample"),
    "timepoint": combined("timepoint"),
    "cell_type": combined("cell_type"),
    "transferred_mouse_pseudotime": old["transferred_mm_dpt_pseudotime_ranknorm"].to_numpy(dtype=float),
    "transferred_mouse_fibdec": pd.Categorical(fibdec.to_numpy(), categories=list(cell_types)),
}, index=pd.Index(names.to_numpy(), dtype=object))

# ---- write it into a copy
shutil.copyfile(TARGET, TEMP)
with h5py.File(TEMP, "r+") as f:
    del f["obs"]
    write_elem(f, "obs", new)

# ---- verify
report = {"written": datetime.datetime.now().isoformat(timespec="seconds"), "source": RECORD, "source_sha256": PIN}
with h5py.File(TEMP, "r") as f:
    back = read_elem(f["obs"])
if list(back.columns) != list(new.columns) or not back.index.astype(str).equals(names):
    fail("obs read back has other columns or another nucleus order")
for column in new.columns:
    a, b = back[column], new[column]
    if a.dtype.kind == "f":
        same = np.array_equal(a.to_numpy(), b.to_numpy(), equal_nan=True)
    else:
        same = a.astype(str).equals(b.astype(str))
    if not same:
        fail(f"obs column {column} read back differs")

hs = back.species.astype(str) == "hs"
counts = {
    "human_with_transferred_pseudotime": int(back.loc[hs, "transferred_mouse_pseudotime"].notna().sum()),
    "human_with_transferred_fibdec": int(back.loc[hs, "transferred_mouse_fibdec"].notna().sum()),
    "transferred_fibdec_counts": {str(k): int(v) for k, v in back.loc[hs, "transferred_mouse_fibdec"].value_counts().items()},
    "mouse_rows_with_any_transfer": int((~hs & (back.transferred_mouse_pseudotime.notna()
                                                | back.transferred_mouse_fibdec.notna())).sum()),
}
expected = {"human_with_transferred_pseudotime": 47352, "human_with_transferred_fibdec": 47509,
            "mouse_rows_with_any_transfer": 0}
if any(counts[k] != v for k, v in expected.items()):
    fail(f"counts {counts} differ from {expected}")

# Independent check: the cell-type vote recomputed from the graph block and the final mouse labels.
check = pd.read_csv(CELL_TYPE_CHECK, index_col="cell")
recomputed = check["transferred_cell_type"].reindex(names[hs.to_numpy()])
if not recomputed.astype(str).equals(back.loc[hs, "transferred_mouse_fibdec"].astype(str).reindex(names[hs.to_numpy()])):
    fail("transferred_mouse_fibdec differs from the graph-block recomputation")
for code in ("hs", "mm"):
    rows = back.species.astype(str) == code
    for column in ("sample", "timepoint", "cell_type"):
        if not back.loc[rows, column].astype(str).equals(final[code][column].astype(str).reindex(back.index[rows])):
            fail(f"{code} {column} differs from the final object")


# Every dataset and attribute outside obs must be byte-identical to the pinned original.
def inventory(path):
    items = {}
    with h5py.File(path, "r") as f:
        items["/"] = ("attrs", sorted((k, str(v)) for k, v in f.attrs.items()))

        def visit(name, obj):
            if name == "obs" or name.startswith("obs/"):
                return
            attrs = sorted((k, str(v)) for k, v in obj.attrs.items())
            if isinstance(obj, h5py.Dataset):
                digest = hashlib.sha256()
                if obj.shape == ():
                    digest.update(np.asarray(obj[()]).tobytes() if obj.dtype.kind != "O" else str(obj[()]).encode())
                else:
                    step = 1 << 22
                    for start in range(0, obj.shape[0], step):
                        block = obj[start:start + step]
                        digest.update(block.tobytes() if obj.dtype.kind != "O" else "\x00".join(map(str, block)).encode())
                items[name] = (str(obj.dtype), obj.shape, attrs, digest.hexdigest())
            else:
                items[name] = ("group", attrs)

        f.visititems(visit)
    return items


before, after = inventory(RECORD), inventory(TEMP)
if before != after:
    changed = sorted(set(before) ^ set(after)) + [k for k in before if k in after and before[k] != after[k]]
    fail(f"non-obs content differs: {changed[:10]}")
report.update({"non_obs_objects_identical": len(before), "obs_columns": list(back.columns), "counts": counts})

# ---- swap in
os.replace(TEMP, TARGET)
report["new_sha256"] = sha256_of(TARGET)
report["bytes"] = os.path.getsize(TARGET)
with open(REPORT, "w") as handle:
    json.dump(report, handle, indent=1)
print(json.dumps(report, indent=1))
