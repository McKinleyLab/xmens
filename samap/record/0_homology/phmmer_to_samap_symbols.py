#!/usr/bin/env python
"""
phmmer_to_samap_symbols.py

Aggregate phmmer .tbl outputs (protein-protein hits) into a SAMap-compatible
homology table at the gene-symbol level. Takes max bit-score per (qsym, tsym)
pair across all isoforms.

Usage:
    python phmmer_to_samap_symbols.py hu_vs_mu [--prot-dir DIR] [--out-dir DIR]
    python phmmer_to_samap_symbols.py mu_vs_hu [--prot-dir DIR] [--out-dir DIR]

Output (SAMap maps convention; <out-dir> defaults to .../homology/hsmm_hmmer/hsmm):
    <out-dir>/hs_to_mm.txt   (hu_vs_mu  →  hs queries vs mm targets)
    <out-dir>/mm_to_hs.txt   (mu_vs_hu  →  mm queries vs hs targets)

2026-09-10: paths made arguments (defaults = live lab_storage locations). The legacy
hardcoded /n/eddy_lab paths are dead, and the legacy default wrote one level above the
hsmm/ dir that SAMap and 09_homology.smk actually use. Parsing/aggregation logic unchanged.
"""
import argparse
import re
import sys
import gzip
from pathlib import Path

GS_RE = re.compile(r'gene_symbol:(\S+)')

DEFAULT_PROT_DIR = '/n/lab_storage/eddy_lab/users/nhilgert/data/proteomes'
DEFAULT_OUT_DIR = ('/n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/'
                   'samap_integrate/homology/hsmm_hmmer/hsmm')
OUT_BASENAME = {'hu_vs_mu': 'hs_to_mm.txt', 'mu_vs_hu': 'mm_to_hs.txt'}


def build_protein_to_symbol(fasta_path):
    p2s = {}
    n_no_symbol = 0
    opener = gzip.open if str(fasta_path).endswith('.gz') else open
    mode = 'rt' if str(fasta_path).endswith('.gz') else 'r'
    with opener(fasta_path, mode) as f:
        for line in f:
            if not line.startswith('>'):
                continue
            pid = line[1:].split(None, 1)[0]
            m = GS_RE.search(line)
            if m:
                p2s[pid] = m.group(1)
            else:
                n_no_symbol += 1
    return p2s, n_no_symbol


def aggregate_chunk(tbl_path, query_p2s, target_p2s):
    pair_max = {}
    n_lines = n_kept = n_dropped_q = n_dropped_t = 0
    with open(tbl_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            n_lines += 1
            cols = line.split()
            if len(cols) < 6:
                continue
            target_pid = cols[0]
            query_pid = cols[2]
            evalue = float(cols[4])
            bitscore = float(cols[5])
            qsym = query_p2s.get(query_pid)
            tsym = target_p2s.get(target_pid)
            if qsym is None:
                n_dropped_q += 1
                continue
            if tsym is None:
                n_dropped_t += 1
                continue
            n_kept += 1
            key = (qsym, tsym)
            prev = pair_max.get(key)
            if prev is None or bitscore > prev[0]:
                pair_max[key] = (bitscore, evalue)
    return pair_max, n_lines, n_kept, n_dropped_q, n_dropped_t


def merge_into(global_max, chunk_max):
    for key, val in chunk_max.items():
        prev = global_max.get(key)
        if prev is None or val[0] > prev[0]:
            global_max[key] = val


def write_samap_format(pair_max, out_path):
    with open(out_path, 'w') as f:
        for (q, t), (bit, ev) in pair_max.items():
            f.write(
                f"{q}\t{t}\t100.0\t100\t0\t0\t1\t100\t1\t100\t{ev}\t{bit}\n"
            )


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('direction', choices=('hu_vs_mu', 'mu_vs_hu'))
    ap.add_argument('--prot-dir', default=DEFAULT_PROT_DIR,
                    help='proteomes root holding fasta/ and hm/<direction>/<prefix>_<i>.tbl')
    ap.add_argument('--out-dir', default=DEFAULT_OUT_DIR, help='directory to write the SAMap map into')
    args = ap.parse_args()

    direction = args.direction
    PROT_DIR = Path(args.prot_dir)
    HMMER_CHUNKS_DIR = PROT_DIR / 'hm'
    OUT_DIR = Path(args.out_dir)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FASTA = {
        'hu': PROT_DIR / 'fasta' / 'Homo_sapiens.GRCh38.pep.all.fa',
        'mu': PROT_DIR / 'fasta' / 'Mus_musculus.GRCm39.pep.all.fa',
    }
    OUT_FILE = {k: OUT_DIR / v for k, v in OUT_BASENAME.items()}
    if direction == 'hu_vs_mu':
        query_species, target_species = 'hu', 'mu'
        chunk_prefix = 'hu'
    else:
        query_species, target_species = 'mu', 'hu'
        chunk_prefix = 'mu'

    print(f"Direction: {direction}")
    print(f"  query  = {query_species}")
    print(f"  target = {target_species}")
    print(f"  output = {OUT_FILE[direction]}")
    print()

    print(f"Building {query_species} protein->symbol map from {FASTA[query_species].name}")
    query_p2s, q_nosym = build_protein_to_symbol(FASTA[query_species])
    print(f"  {len(query_p2s):,} proteins with symbol  ({q_nosym:,} without)")

    print(f"Building {target_species} protein->symbol map from {FASTA[target_species].name}")
    target_p2s, t_nosym = build_protein_to_symbol(FASTA[target_species])
    print(f"  {len(target_p2s):,} proteins with symbol  ({t_nosym:,} without)")
    print()

    chunks_dir = HMMER_CHUNKS_DIR / direction
    if not chunks_dir.is_dir():
        sys.exit(f"missing chunks dir: {chunks_dir}")

    pair_max = {}
    totals = dict(lines=0, kept=0, drop_q=0, drop_t=0)

    for i in range(1, 33):
        tbl = chunks_dir / f"{chunk_prefix}_{i}.tbl"
        if not tbl.is_file():
            print(f"  WARNING missing chunk: {tbl}")
            continue
        chunk_max, n_lines, n_kept, n_dq, n_dt = aggregate_chunk(
            tbl, query_p2s, target_p2s
        )
        merge_into(pair_max, chunk_max)
        totals['lines'] += n_lines
        totals['kept'] += n_kept
        totals['drop_q'] += n_dq
        totals['drop_t'] += n_dt
        if i % 8 == 0 or i == 32:
            print(f"  chunk {i:2d}: {len(pair_max):>8,} unique symbol-pairs accumulated")

    print()
    print(f"Total hit lines parsed: {totals['lines']:,}")
    print(f"  kept  (both ends symbolized): {totals['kept']:,}")
    print(f"  dropped (query had no symbol): {totals['drop_q']:,}")
    print(f"  dropped (target had no symbol): {totals['drop_t']:,}")
    print(f"Unique (qsym, tsym) pairs: {len(pair_max):,}")
    print()

    out_path = OUT_FILE[direction]
    print(f"Writing {out_path}")
    write_samap_format(pair_max, out_path)
    print(f"Wrote {len(pair_max):,} lines.")
    print("Done.")


if __name__ == '__main__':
    main()
