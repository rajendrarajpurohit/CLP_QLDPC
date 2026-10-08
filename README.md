# Coset lifted product (CLP) quantum LDPC codes, search and certification release

This repository holds every script, pair list, results file and per code certificate
behind the codes reported in the paper. It is a frozen snapshot with checksums in
MANIFEST.md5.

## Certified codes

| Code          | w | Group               | H                    | Folder                    |
|---------------|---|---------------------|----------------------|---------------------------|
| [[168,10,14]] | 6 | SmallGroup(252,21)  | `<f2 f5, f4>` = C6   | codes/168_10_14           |
| [[240,8,20]]  | 6 | SmallGroup(360,99)  | `<f2, f3>` = C6      | codes/240_8_20            |
| [[360,20,20]] | 6 | SmallGroup(270,26)  | `<f2 f4>` = C3       | codes/360_20_20           |
| [[360,24,20]] | 7 | SmallGroup(810,84)  | `<f6, f2 f4>` = C9   | codes/360_24_20           |
| [[216,8,18]]  | 7 | SmallGroup(162,3)   | `<f3 f5>` = C3       | codes/phase1_extras/D2_10 |
| [[252,16,18]] | 9 | SmallGroup(189,8)   | `<f2>` = C3          | codes/252_16_18           |
| [[288,16,20]] | 9 | SmallGroup(144,164) | `<f1 f5>` = C2       | codes/288_16_20           |

Each folder holds HX and HZ as 0/1 text, the construction record (G, H, A, B), the
solver logs for every sector on both sides, the distance witnesses, the output of
python/verify.py, and a PROVENANCE.txt saying which log proves which sector.
Ten further certified codes with smaller parameters are in codes/phase1_extras.

## Layout

    gap/        CLP_core.g (construction), rebuild.g (build a code from its printed description),
                NOVA_lib.g, rescreen.g, recoverD3.g (internal rebuild check for the Phase 1 codes)
    python/     certify.py, sectors.py, witness.py, verify.py, collect.py
    pairs/      (G,H) pair lists used by the Phase 4 drivers
    drivers/    W7 to W10 drivers and F9, plus G1 and G2 from Phase 3
    results/    raw results files from every Phase 3 and 4 driver
    phase1/     Phase 1 (D1 to D4) results and the certified D3 matrices, run on a personal machine
    phase2/     Phase 2 wide search scripts and results, including the H10 recovery chain
    hope/       F1 to F5 drivers between Phases 2 and 3
    negative/   layouts G3 to G9, archived, negative results
    codes/      one folder per certified code
    ORIGINS.md  which driver and results file produced every code in the paper

## Versions

GAP 4.15.1 with the SmallGroups library and QDistRnd. Python 3.13.9, NumPy 2.3.5,
Gurobi 13.0.3 through gurobipy. Runs on Ubuntu, 64 cores, 125 GB RAM.

## Reproducing a code from the paper

Every code in the paper is printed as a SmallGroup identifier, the generators of H,
and the entries of A and B as words in the pc generators f1, f2, ... of that group.
`gap/rebuild.g` turns such a block into parity check matrices using the same
construction functions the search used (`CLP_core.g`). Nothing else is needed.

### 1. Transcribe the block into a call

    Rebuild( l, i, Hw, Aw, Bw, sA, sB, name );

| argument | meaning |
|----------|---------|
| `l, i`   | the SmallGroup identifier |
| `Hw`     | `f -> [ generators of H ]` |
| `Aw`     | `f -> A` as a list of rows, each entry a list of its terms |
| `Bw`     | `f -> B` in the same form |
| `sA, sB` | `[1,1]` and `[2,2]` for the 1x1 over 2x2 layout, `[2,2]` and `[1,1]` for 2x2 over 1x1 |
| `name`   | prefix for the output files `name_HX.txt`, `name_HZ.txt` |

An entry equal to 1 is written `One(f[1])`. A sum such as 1 + f6 + f4^4 f6^2 is the
list `[ One(f[1]), f[6], f[4]^4*f[6]^2 ]`. The pc generators f1, f2, ... are those of
GAP's `SmallGroup(l,i)`, obtained inside the call as `Pcgs(G)`.

A worked example, the [[120,8,8]] code of the paper, is `gap/example_120_8_8.g`.

### 2. Build the matrices

    cd gap
    gap -q example_120_8_8.g

The script prints n, k, the check weight and a QDistRnd upper bound, and writes
`example_HX.txt` and `example_HZ.txt` in the current directory. If a transcription
error puts an entry of B outside N_G(H), or an entry of A in the core, it stops with
a message naming which matrix is wrong.

### 3. Check the structure and certify the distance

    python3 ../python/verify.py  example 120 8 --dir .
    python3 ../python/certify.py example 120 8 --dir .

`verify.py` checks orthogonality, dimension, check weights, qubit degrees and the
logical bases from the matrix files alone. `certify.py` proves d >= 8 by solving one
integer program per logical sector with Gurobi, cutoff 7.5, then finds a weight 8
logical operator, which gives d = 8. The 120 qubit code takes seconds. The length
360 codes with d = 20 take days on a workstation. `python/sectors.py` runs a chosen
range of sectors on one side so the work can be split across machines, and
`python/witness.py` extracts the weight d logical operator separately.

## Rerunning a stored certification

    python3 python/certify.py <id> <n> <d_ub> --dir <folder> --threads 24 --ram 40

runs every sector on both sides with a cutoff at d_ub minus one half. A sector that
reaches the cutoff proves no logical operator of weight below d_ub anticommutes with
that basis element.

## Notes

Phase 1 matrix files for n >= 216 contain GAP backslash line continuations. Strip them
before parsing with a strict reader. The D3_results.txt kept on the server mixes two runs
with restarted counters. phase1/D3_results_original.txt is the clean single run file.
For the Phase 1 codes, `gap/recoverD3.g` rebuilds each one from its printed (G, H, A, B)
and confirms the result is identical to the stored matrices. For H10 the same check is
`phase2/recoverH10.g`.
