#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
certify.py - exact CSS distance + witness, with logging and resume.

usage:
    python certify.py a168_20_10 168 10 --dir .
    python certify.py b336_12_23_ii 336 23 --dir . --threads 24 --ram 8
"""
import os
import sys
import json
import time
import argparse
import numpy as np
import gurobipy as gp
from gurobipy import GRB

MEMLIMIT_GB = 1000


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
        res = (d_ub, None)
    elif st == GRB.OPTIMAL:
        res = (int(round(md.ObjVal)), np.round(e.X).astype(np.int64))
    elif st == GRB.INFEASIBLE:
        res = (10 ** 9, None)
    elif st == GRB.MEM_LIMIT:
        b = md.ObjBound
        md.dispose()
        sys.exit("%s: hit %.0f GB soft limit (bound %.2f); "
                 "rerun with fewer --threads" % (name, MEMLIMIT_GB, b))
    elif st == GRB.INTERRUPTED:
        b = md.ObjBound
        md.dispose()
        sys.exit("%s: interrupted (bound %.2f); finished sectors are saved"
                 % (name, b))
    else:
        md.dispose()
        sys.exit("%s: unexpected Gurobi status %d" % (name, st))
    md.dispose()
    return res


def find_witness(H, logicals, target, threads):
    N = H.shape[1]
    rw = H.sum(axis=1)
    for i, lg in enumerate(logicals):
        md = gp.Model("wit%d" % i)
        md.Params.OutputFlag = 0
        md.Params.Threads = threads
        md.Params.SoftMemLimit = MEMLIMIT_GB
        md.Params.BestObjStop = target + 0.5
        e = md.addMVar(N, vtype=GRB.BINARY)
        z = md.addMVar(H.shape[0], vtype=GRB.INTEGER, lb=0, ub=(rw // 2))
        w = md.addVar(vtype=GRB.INTEGER, lb=0, ub=int(lg.sum()) // 2)
        md.addConstr(H @ e == 2 * z)
        md.addConstr(lg @ e == 1 + 2 * w)
        md.setObjective(e.sum(), GRB.MINIMIZE)
        md.optimize()
        got = md.SolCount > 0
        v = np.round(e.X).astype(np.int64) if got else None
        md.dispose()
        if got and int(v.sum()) == target:
            return v, i
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("id")
    ap.add_argument("N", type=int)
    ap.add_argument("d_ub", type=int)
    ap.add_argument("--dir", default=".")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--ram", type=float, default=8.0)
    a = ap.parse_args()

    base = os.path.join(a.dir, a.id)
    logfile = "%s_certify.log" % a.id
    gurobi_log = "%s_gurobi.log" % a.id
    prog_f = "%s_progress.json" % a.id

    HX = parse_matrix(base + "_HX.txt", a.N)
    HZ = parse_matrix(base + "_HZ.txt", a.N)
    if ((HX @ HZ.T) % 2).any():
        sys.exit("CSS orthogonality fails")

    LX = logical_basis(HX, HZ)
    LZ = logical_basis(HZ, HX)
    k = a.N - gf2_rref(HX)[0] - gf2_rref(HZ)[0]
    if not (len(LX) == len(LZ) == k):
        sys.exit("basis size %d/%d != k=%d" % (len(LX), len(LZ), k))

    np.savetxt("%s_LX_basis.txt" % a.id, np.array(LX), fmt="%d")
    np.savetxt("%s_LZ_basis.txt" % a.id, np.array(LZ), fmt="%d")

    wt = int(max(HX.sum(axis=1).max(), HZ.sum(axis=1).max(),
                 (HX.sum(axis=0) + HZ.sum(axis=0)).max()))

    log(logfile, "=== %s [[%d,%d]] w=%d d_ub=%d Xlog=%d Zlog=%d "
                 "threads=%d nodefile=%.0fGB memlimit=%.0fGB"
        % (a.id, a.N, k, wt, a.d_ub, len(LX), len(LZ),
           a.threads, a.ram, MEMLIMIT_GB))

    if os.path.exists(prog_f):
        prog = json.load(open(prog_f))
    else:
        prog = {"id": a.id, "N": a.N, "k": k, "w": wt,
                "d_ub": a.d_ub, "dZ": [], "dX": []}

    for tag, H, L in (("dZ", HX, LX), ("dX", HZ, LZ)):
        for i in range(len(prog[tag]), len(L)):
            t0 = time.time()
            val, wit = prove_sector(H, L[i], "%s_%d" % (tag, i + 1),
                                    prog["d_ub"], a.threads, a.ram, gurobi_log)
            if wit is not None:
                if ((H @ wit) % 2).any():
                    sys.exit("witness not in kernel")
                np.savetxt("%s_%s_improved_w%d.txt" % (a.id, tag, val),
                           wit, fmt="%d")
                prog["d_ub"] = val
                log(logfile, "%s sector %d: FOUND weight %d, new bound"
                    % (tag, i + 1, val))
            else:
                log(logfile, "%s sector %d/%d: >= %d  (%.0fs)"
                    % (tag, i + 1, len(L), prog["d_ub"], time.time() - t0))
            prog[tag].append(val)
            json.dump(prog, open(prog_f, "w"), indent=1)

    d = min(min(prog["dZ"]), min(prog["dX"]))
    log(logfile, "PROVEN d = %d" % d)

    if min(prog["dZ"]) <= min(prog["dX"]):
        H, L, tag = HX, LX, "dZ"
    else:
        H, L, tag = HZ, LZ, "dX"

    v, idx = find_witness(H, L, d, a.threads)
    if v is None:
        log(logfile, "WARNING: no witness of weight %d found" % d)
    else:
        fn = "%s_witness_%s_w%d.txt" % (a.id, tag, d)
        np.savetxt(fn, v, fmt="%d")
        supp = np.nonzero(v)[0].tolist()
        log(logfile, "witness %s sector %d weight %d support %s"
            % (tag, idx + 1, int(v.sum()), supp))
        prog["witness_file"] = fn
        prog["witness_support"] = supp

    prog["d_exact"] = d
    prog["Q"] = float(k * d * d) / a.N
    json.dump(prog, open(prog_f, "w"), indent=1)
    log(logfile, "=== FINAL [[%d,%d,%d]] w=%d Q=%.3f"
        % (a.N, k, d, wt, prog["Q"]))


if __name__ == "__main__":
    main()