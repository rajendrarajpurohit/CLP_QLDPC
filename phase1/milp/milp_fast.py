import sys, json, argparse, numpy as np, gurobipy as gp
from gurobipy import GRB

def parse_matrix(path, N):
    d = [int(c) for c in open(path).read() if c in '01']
    return np.array(d, dtype=int).reshape(-1, N)

def gf2_rref(A):
    A = A.copy() % 2; m, n = A.shape; r = 0; piv = []
    for c in range(n):
        rows = np.nonzero(A[r:, c])[0]
        if len(rows) == 0: continue
        p = r + rows[0]; A[[r, p]] = A[[p, r]]; piv.append(c)
        mask = A[:, c].astype(bool); mask[r] = False
        A[mask] ^= A[r]; r += 1
        if r == m: break
    return r, A, piv

def logical_basis(H_same, H_other):
    n = H_other.shape[1]; r, R, piv = gf2_rref(H_other)
    free = [c for c in range(n) if c not in piv]
    ker = []
    for f in free:
        v = np.zeros(n, dtype=int); v[f] = 1
        for i, p in enumerate(piv):
            if R[i, f]: v[p] = 1
        ker.append(v)
    base = [x for x in H_same % 2 if x.any()]; rank = gf2_rref(np.array(base))[0]; logs = []
    for v in ker:
        r2 = gf2_rref(np.array(base + [v]))[0]
        if r2 > rank: logs.append(v); base.append(v); rank = r2
    return logs

def sector(H, lg, name, d_ub, threads, ramGB):
    N = H.shape[1]; rw = H.sum(axis=1)
    md = gp.Model(name); P = md.Params
    P.OutputFlag = 1; P.Threads = threads; P.NodefileStart = ramGB; P.NodefileDir = "."
    P.Cutoff = d_ub - 0.5; P.MIPFocus = 3; P.Symmetry = 2
    e = md.addMVar(N, vtype=GRB.BINARY)
    z = md.addMVar(H.shape[0], vtype=GRB.INTEGER, lb=0, ub=rw // 2)
    w = md.addVar(vtype=GRB.INTEGER, lb=0, ub=int(lg.sum()) // 2)
    md.addConstr(H @ e == 2 * z); md.addConstr(lg @ e == 1 + 2 * w)
    md.setObjective(e.sum(), GRB.MINIMIZE); md.optimize()
    if md.Status == GRB.CUTOFF:   val, wit = d_ub, None            
    elif md.Status == GRB.OPTIMAL: val, wit = int(round(md.ObjVal)), np.round(e.X).astype(int)
    else: raise RuntimeError(f"{name}: status {md.Status}")
    md.dispose(); return val, wit

def witness(H, lg, target, threads):
    N = H.shape[1]; rw = H.sum(axis=1)
    md = gp.Model("wit"); P = md.Params; P.OutputFlag = 0; P.Threads = threads; P.BestObjStop = target
    e = md.addMVar(N, vtype=GRB.BINARY); z = md.addMVar(H.shape[0], vtype=GRB.INTEGER, lb=0, ub=rw // 2)
    w = md.addVar(vtype=GRB.INTEGER, lb=0, ub=int(lg.sum()) // 2)
    md.addConstr(H @ e == 2 * z); md.addConstr(lg @ e == 1 + 2 * w)
    md.setObjective(e.sum(), GRB.MINIMIZE); md.optimize()
    v = np.round(e.X).astype(int); md.dispose(); return v

def run_side(H, logs, tag, d_ub, prog, threads, ramGB):
    done = prog[tag]
    for i in range(len(done), len(logs)):
        val, wit = sector(H, logs[i], f"{tag}_{i+1}", d_ub, threads, ramGB)
        print(f">>> {tag} sector {i+1}/{len(logs)}: {'>= '+str(d_ub) if wit is None else '= '+str(val)}")
        if wit is not None:
            assert not ((H @ wit) % 2).any(); np.save(f"{prog['id']}_{tag}_witness.npy", wit); d_ub = val
        done.append(val); prog["d_ub"] = d_ub; json.dump(prog, open(prog["file"], "w"))
    return min(done + [d_ub])

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("id"); ap.add_argument("N", type=int)
    ap.add_argument("d_ub", type=int); ap.add_argument("--done-z", type=int, default=0)
    ap.add_argument("--threads", type=int, default=8); ap.add_argument("--ram", type=float, default=35.0)
    a = ap.parse_args()
    HX = parse_matrix(f"{a.id}_HX.txt", a.N); HZ = parse_matrix(f"{a.id}_HZ.txt", a.N)
    assert not ((HX @ HZ.T) % 2).any(), "HX HZ^T != 0"
    LX, LZ = logical_basis(HX, HZ), logical_basis(HZ, HX)
    k = a.N - gf2_rref(HX)[0] - gf2_rref(HZ)[0]
    assert len(LX) == len(LZ) == k, f"basis size {len(LX)},{len(LZ)} != k={k}"
    print(f"[[{a.N},{k}]]  d_ub={a.d_ub}  X-logicals={len(LX)}  Z-logicals={len(LZ)}")
    pf = f"progress_{a.id}.json"
    try: prog = json.load(open(pf))
    except FileNotFoundError:
        prog = {"id": a.id, "file": pf, "d_ub": a.d_ub, "dZ": [a.d_ub] * a.done_z, "dX": []}
    dZ = run_side(HX, LX, "dZ", prog["d_ub"], prog, a.threads, a.ram)   
    dX = run_side(HZ, LZ, "dX", min(prog["d_ub"], dZ), prog, a.threads, a.ram)
    d = min(dZ, dX); print(f"\nEXACT: [[{a.N},{k},{d}]]   dZ={dZ} dX={dX}")
    H, L, tag = (HX, LX, "dZ") if dZ <= dX else (HZ, LZ, "dX")
    for lg in L:
        v = witness(H, lg, d, a.threads)
        if v.sum() == d: np.save(f"{a.id}_{tag}_witness_w{d}.npy", v); print(f"witness saved ({tag}, weight {d})"); break

if __name__ == "__main__": main()
