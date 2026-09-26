import json, argparse, numpy as np, gurobipy as gp
from gurobipy import GRB


def parse_matrix(path, N):
    d = [int(c) for c in open(path).read() if c in '01']
    return np.array(d, dtype=int).reshape(-1, N)


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
    r, R, piv = gf2_rref(H_other)
    free = [c for c in range(n) if c not in piv]
    ker = []
    for f in free:
        v = np.zeros(n, dtype=int)
        v[f] = 1
        for i, p in enumerate(piv):
            if R[i, f]:
                v[p] = 1
        ker.append(v)
    base = [x for x in H_same % 2 if x.any()]
    rank = gf2_rref(np.array(base))[0]
    logs = []
    for v in ker:
        r2 = gf2_rref(np.array(base + [v]))[0]
        if r2 > rank:
            logs.append(v)
            base.append(v)
            rank = r2
    return logs


def sector(H, lg, name, d_ub, threads, ramGB):
    N = H.shape[1]
    rw = H.sum(axis=1)
    md = gp.Model(name)
    P = md.Params
    P.OutputFlag = 1
    P.LogFile = name + ".log"
    P.Threads = threads
    P.NodefileStart = ramGB
    P.NodefileDir = "."
    P.Cutoff = d_ub - 0.5
    P.MIPGap = 0
    P.MIPFocus = 3
    P.Symmetry = 2
    e = md.addMVar(N, vtype=GRB.BINARY)
    z = md.addMVar(H.shape[0], vtype=GRB.INTEGER, lb=0, ub=rw // 2)
    w = md.addVar(vtype=GRB.INTEGER, lb=0, ub=int(lg.sum()) // 2)
    md.addConstr(H @ e == 2 * z)
    md.addConstr(lg @ e == 1 + 2 * w)
    md.setObjective(e.sum(), GRB.MINIMIZE)
    md.optimize()
    if md.Status == GRB.CUTOFF:
        res = {"status": "CUTOFF", "min_weight": ">=" + str(d_ub), "bound": md.ObjBound}
    elif md.Status == GRB.OPTIMAL:
        wit = np.round(e.X).astype(int)
        assert not ((H @ wit) % 2).any()
        assert (lg @ wit) % 2 == 1
        np.save(name + "_witness.npy", wit)
        res = {"status": "OPTIMAL", "min_weight": int(round(md.ObjVal))}
    else:
        raise RuntimeError(name + ": status " + str(md.Status))
    res["nodes"] = md.NodeCount
    res["seconds"] = md.Runtime
    res["gurobi"] = ".".join(str(x) for x in gp.gurobi.version())
    md.dispose()
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("id")
    ap.add_argument("N", type=int)
    ap.add_argument("d_ub", type=int)
    ap.add_argument("--first", type=int, default=1)
    ap.add_argument("--last", type=int, default=5)
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--ram", type=float, default=35.0)
    a = ap.parse_args()

    HX = parse_matrix(a.id + "_HX.txt", a.N)
    HZ = parse_matrix(a.id + "_HZ.txt", a.N)
    assert not ((HX @ HZ.T) % 2).any(), "HX HZ^T != 0"
    LX = logical_basis(HX, HZ)
    k = a.N - gf2_rref(HX)[0] - gf2_rref(HZ)[0]
    assert len(LX) == k, "basis size mismatch"
    print("[[%d,%d]]  Z-sectors %d..%d of %d, cutoff %.1f" % (a.N, k, a.first, a.last, k, a.d_ub - 0.5))

    out = "%s_dZ_sectors_%d-%d.json" % (a.id, a.first, a.last)
    try:
        results = json.load(open(out))
    except FileNotFoundError:
        results = {}
    for i in range(a.first, a.last + 1):
        if str(i) in results:
            print("sector %d already done: %s" % (i, results[str(i)]))
            continue
        results[str(i)] = sector(HX, LX[i - 1], "%s_dZ_%d" % (a.id, i), a.d_ub, a.threads, a.ram)
        print(">>> dZ sector %d/%d: %s" % (i, k, results[str(i)]))
        json.dump(results, open(out, "w"), indent=1)
    print("all sectors:", {i: r["min_weight"] for i, r in results.items()})


if __name__ == "__main__":
    main()