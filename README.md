# Coset lifted product (CLP) quantum LDPC codes, search and certification release

This repository holds every script, pair list, results file and per code certificate
behind the codes reported in the paper. Nothing here is run from this folder. It is a
frozen snapshot with checksums in MANIFEST.md5.

## Certified codes

| Code | w | Group | H | Folder |
|---|---|---|---|---|
| [[168,10,14]] | 6 | SmallGroup(252,21) | <f2 f5, f4> = C6 | codes/168_10_14 |
| [[240,8,20]]  | 6 | SmallGroup(360,99) | <f2, f3> = C6    | codes/240_8_20  |
| [[360,20,20]] | 6 | SmallGroup(270,26) | <f2 f4> = C3     | codes/360_20_20 |
| [[252,16,18]] | 9 | SmallGroup(189,8)  | <f2> = C3        | codes/252_16_18 |

Each folder holds HX and HZ as 0/1 text, the construction record (G, H, a, B), the
solver logs for every sector on both sides, the distance witnesses, the output of
python/verify.py, and a PROVENANCE.txt saying which log proves which sector.

## Layout

    gap/        CLP_core.g (construction), NOVA_lib.g, rescreen.g, recoverD3.g (rebuild check)
    python/     certify.py, sectors.py, witness.py, verify.py, collect.py
    pairs/      (G,H) pair lists used by the Phase 4 drivers
    drivers/    W7 to W10 drivers and F9, plus G1 and G2 from Phase 3
    results/    raw results files from every Phase 3 and 4 driver
    phase1/     Phase 1 (D1 to D4) results and the certified D3 matrices, run on a personal machine
    phase2/     Phase 2 wide search scripts and results, including the H10 recovery chain
    hope/       F1 to F5 drivers between Phases 2 and 3
    negative/   layouts G3 to G9, archived, negative results
    codes/      one folder per certified or partially certified code
    ORIGINS.md  which driver and results file produced every code in the paper

## Versions

GAP 4.15.1 with the SmallGroups library and QDistRnd. Python 3.13.9, NumPy 2.3.5,
Gurobi 13.0.3 through gurobipy. Runs on Ubuntu, 64 cores, 125 GB RAM.

## Reproducing a code from its record

Every record gives G as a SmallGroup id, H by generators written in the polycyclic
generating sequence f1, f2, ... that GAP returns for that group, and a and B as sums of
such words. gap/recoverD3.g shows the full procedure. It rebuilds HX and HZ from
(G, H, a, B) with CLP_core.g and compares them with the stored files. For H10 the
same check is in phase2/recoverH10.g and its output is codes/360_20_20/H10_package/H10_elements.txt.

## Rerunning a certification

    python3 python/certify.py <id> <n> <d_ub> --dir <folder> --threads 24 --ram 40

runs every sector on both sides with a cutoff at d_ub minus one half. A sector that
reaches the cutoff proves no logical operator of weight below d_ub anticommutes with
that basis element. python/sectors.py runs a chosen range of sectors on one side, and
python/witness.py extracts a weight d_ub logical operator so the distance is exact.

## Notes

Phase 1 matrix files for n >= 216 contain GAP backslash line continuations. Strip them
before parsing with a strict reader. The D3_results.txt kept on the server mixes two runs
with restarted counters. phase1/D3_results_original.txt is the clean single run file.
