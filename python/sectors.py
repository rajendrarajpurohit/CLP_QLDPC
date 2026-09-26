#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sectors.py - prove a range of sectors for one side of a CSS code.
Same Gurobi configuration as certify.py.

usage:
    python sectors.py b336_12_23_ii 336 23 --side Z --first 3 --last 5 --threads 8 --ram 5
"""
import os
import sys
import json
import time
import argparse
import numpy as np
import gurobipy as gp
from gurobipy import GRB

MEMLIMIT_GB = 60.0


def parse_matrix(path, N):
    with open(path) as f:
        d = [int(c) for c in f.read() if c in '01']
    if len(d) % N:
        sys.exit("%s: %d digits, not a multiple of %d" % (path, len(d), N))
    return np.array(d, dtype=np.int64).reshape(-1, N)


def gf2_rref(A):
    A = A.copy() % 2
    m, n = A.shape
    r = 0
    piv = []
    for c in range(n):
        rows = np.nonzero(A[r:, c])[0]
        if len(rows) == 0:
            continue
        p = r + rows[0]
        if p != r:
            A[[r, p]] = A[[p, r]]
        piv.append(c)
        mask = A[:, c].astype(bool)
        mask[r] = False
        A[mask] ^= A[r]
        r += 1
        if r == m:
            break
    return r, A, piv


def logical_basis(H_same, H_other):
    n = H_other.shape[1]
    _, R, piv = gf2_rref(H_other)
    free = [c for c in range(n) if c not in piv]
    ker = []
    for f in free:
        v = np.zeros(n, dtype=np.int64)
        v[f] = 1
        for i, p in enumerate(piv):
            if R[i, f]:
                v[p] = 1
        ker.append(v)
    base = [x for x in (H_same % 2) if x.any()]
    rank = gf2_rref(np.array(base))[0] if base else 0
    out = []
    for v in ker:
        r2 = gf2_rref(np.array(base + [v]))[0]
        if r2 > rank:
            out.append(v)
            base.append(v)
            rank = r2
    return out


def log(logfile, msg):
    line = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), msg)
    print(line, flush=True)
    with open(logfile, "a") as f:
        f.write(line + "\n")


def prove_sector(H, lg, name, d_ub, threads, ramGB, gurobi_log):
    N = H.shape[1]
    rw = H.sum(axis=1)
    md = gp.Model(name)
    md.Params.OutputFlag = 1
    md.Params.LogFile = gurobi_log
    md.Params.Threads = threads
    md.Params.NodefileStart = ramGB
    md.Params.NodefileDir = "."
    md.Params.SoftMemLimit = MEMLIMIT_GB
    md.Params.Cutoff = d_ub - 0.5
    md.Params.MIPFocus = 3
    md.Params.Symmetry = 2
    e = md.addMVar(N, vtype=GRB.BINARY, name="e")
    z = md.addMVar(H.shape[0], vtype=GRB.INTEGER, lb=0, ub=(rw // 2))
    w = md.addVar(vtype=GRB.INTEGER, lb=0, ub=int(lg.sum()) // 2)
    md.addConstr(H @ e == 2 * z)
    md.addConstr(lg @ e == 1 + 2 * w)
    md.setObjective(e.sum(), GRB.MINIMIZE)
    md.optimize()
    st = md.Status
    if st == GRB.CUTOFF:
        res = {"status": "CUTOFF", "min_weight": ">=%d" % d_ub,
               "bound": md.ObjBound}
        wit = None
    elif st == GRB.OPTIMAL:
        v = int(round(md.ObjVal))
        wit = np.round(e.X).astype(np.int64)
        res = {"status": "OPTIMAL", "min_weight": v}
    elif st == GRB.INFEASIBLE:
        res = {"status": "INFEASIBLE", "min_weight": None}
        wit = None
    elif st == GRB.MEM_LIMIT:
        b = md.ObjBound
        md.dispose()
        sys.exit("%s: hit %.0f GB soft limit (bound %.2f)"
                 % (name, MEMLIMIT_GB, b))
    elif st == GRB.INTERRUPTED:
        b = md.ObjBound
        md.dispose()
        sys.exit("%s: interrupted (bound %.2f); finished sectors saved"
                 % (name, b))
    else:
        md.dispose()
        sys.exit("%s: unexpected Gurobi status %d" % (name, st))
    res["nodes"] = md.NodeCount
    res["seconds"] = md.Runtime
    res["threads"] = threads
    md.dispose()
    return res, wit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("id")
    ap.add_argument("N", type=int)
    ap.add_argument("d_ub", type=int)
    ap.add_argument("--side", choices=["Z", "X"], default="Z")
    ap.add_argument("--first", type=int, default=1)
    ap.add_argument("--last", type=int, default=1)
    ap.add_argument("--dir", default=".")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--ram", type=float, default=5.0)
    a = ap.parse_args()

    base = os.path.join(a.dir, a.id)
    tag = "d" + a.side
    stem = "%s_%s_%d-%d" % (a.id, tag, a.first, a.last)
    logfile = stem + ".log"
    gurobi_log = stem + "_gurobi.log"
    out_f = stem + ".json"

    HX = parse_matrix(base + "_HX.txt", a.N)
    HZ = parse_matrix(base + "_HZ.txt", a.N)
    if ((HX @ HZ.T) % 2).any():
        sys.exit("CSS orthogonality fails")

    if a.side == "Z":
        H, Ho = HX, HZ
    else:
        H, Ho = HZ, HX
    L = logical_basis(H, Ho)
    k = a.N - gf2_rref(HX)[0] - gf2_rref(HZ)[0]
    if len(L) != k:
        sys.exit("basis size %d != k=%d" % (len(L), k))
    if a.last > k:
        sys.exit("--last %d exceeds k=%d" % (a.last, k))

    log(logfile, "=== %s [[%d,%d]] %s sectors %d-%d of %d  cutoff=%.1f "
                 "threads=%d nodefile=%.0fGB"
        % (a.id, a.N, k, tag, a.first, a.last, k,
           a.d_ub - 0.5, a.threads, a.ram))

    if os.path.exists(out_f):
        results = json.load(open(out_f))
    else:
        results = {}

    for i in range(a.first, a.last + 1):
        if str(i) in results:
            log(logfile, "%s sector %d already done: %s"
                % (tag, i, results[str(i)]["min_weight"]))
            continue
        t0 = time.time()
        res, wit = prove_sector(H, L[i - 1], "%s_%s_%d" % (a.id, tag, i),
                                a.d_ub, a.threads, a.ram, gurobi_log)
        if wit is not None:
            if ((H @ wit) % 2).any():
                sys.exit("witness not in kernel")
            if (int(lg_dot(L[i - 1], wit))) % 2 != 1:
                sys.exit("witness does not anticommute")
            fn = "%s_%s_%d_witness_w%d.txt" % (a.id, tag, i, res["min_weight"])
            np.savetxt(fn, wit, fmt="%d")
            res["witness_file"] = fn
            log(logfile, "%s sector %d: FOUND weight %d  <-- below cutoff"
                % (tag, i, res["min_weight"]))
        else:
            log(logfile, "%s sector %d: %s  (%.0fs, %.0f nodes)"
                % (tag, i, res["min_weight"], res["seconds"], res["nodes"]))
        results[str(i)] = res
        json.dump(results, open(out_f, "w"), indent=1)

    log(logfile, "range complete: %s"
        % {i: r["min_weight"] for i, r in results.items()})


def lg_dot(lg, v):
    return int(np.dot(lg, v))


if __name__ == "__main__":
    main()