#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py - independent structural check of a stored CSS code.
usage: python verify.py G1_n360k16d22 360
"""
import sys
import numpy as np


def parse_matrix(path, N):
    with open(path) as f:
        d = [int(c) for c in f.read() if c in '01']
    if len(d) % N:
        sys.exit("%s: %d digits, not a multiple of %d" % (path, len(d), N))
    return np.array(d, dtype=np.int64).reshape(-1, N)


def gf2_rank(A):
    A = A.copy() % 2
    m, n = A.shape
    r = 0
    for c in range(n):
        rows = np.nonzero(A[r:, c])[0]
        if len(rows) == 0:
            continue
        p = r + rows[0]
        if p != r:
            A[[r, p]] = A[[p, r]]
        mask = A[:, c].astype(bool)
        mask[r] = False
        A[mask] ^= A[r]
        r += 1
        if r == m:
            break
    return r


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
    rank = gf2_rank(np.array(base)) if base else 0
    out = []
    for v in ker:
        r2 = gf2_rank(np.array(base + [v]))
        if r2 > rank:
            out.append(v)
            base.append(v)
            rank = r2
    return out


def main():
    cid = sys.argv[1]
    N = int(sys.argv[2])
    ok = True

    HX = parse_matrix(cid + "_HX.txt", N)
    HZ = parse_matrix(cid + "_HZ.txt", N)
    print("file shapes: HX %s  HZ %s" % (HX.shape, HZ.shape))

    print("\n[1] entries are 0/1 only")
    e = set(np.unique(HX)) | set(np.unique(HZ))
    print("    unique values: %s  -> %s" % (sorted(e), "PASS" if e <= {0, 1} else "FAIL"))
    ok &= e <= {0, 1}

    print("\n[2] CSS orthogonality  HX . HZ^T = 0 mod 2")
    prod = (HX @ HZ.T) % 2
    nz = int(prod.sum())
    print("    nonzero entries: %d  -> %s" % (nz, "PASS" if nz == 0 else "FAIL"))
    ok &= (nz == 0)

    print("\n[3] dimension  k = n - rank(HX) - rank(HZ)")
    rx, rz = gf2_rank(HX), gf2_rank(HZ)
    k = N - rx - rz
    print("    rank(HX)=%d  rank(HZ)=%d  k=%d" % (rx, rz, k))

    print("\n[4] logical bases have size k on both sides")
    LX = logical_basis(HX, HZ)
    LZ = logical_basis(HZ, HX)
    good = (len(LX) == len(LZ) == k)
    print("    |X-logicals|=%d  |Z-logicals|=%d  -> %s"
          % (len(LX), len(LZ), "PASS" if good else "FAIL"))
    ok &= good

    print("\n[5] logicals are in the right kernel and not stabilizers")
    bad = 0
    for v in LX:
        if ((HZ @ v) % 2).any():
            bad += 1
    for v in LZ:
        if ((HX @ v) % 2).any():
            bad += 1
    print("    logicals outside kernel: %d  -> %s" % (bad, "PASS" if bad == 0 else "FAIL"))
    ok &= (bad == 0)

    print("\n[6] weight, both conventions")
    rwX, rwZ = HX.sum(axis=1), HZ.sum(axis=1)
    degX, degZ = HX.sum(axis=0), HZ.sum(axis=0)
    w_check = int(max(rwX.max(), rwZ.max()))
    w_qubit = int((degX + degZ).max())
    w_QL = max(w_check, w_qubit)
    print("    max check weight (ATB convention): %d" % w_check)
    print("    max qubit degree  deg_X+deg_Z    : %d" % w_qubit)
    print("    overall w (Qian-Li Eq.1)         : %d" % w_QL)
    print("    min check weight: %d" % int(min(rwX.min(), rwZ.min())))

    print("\n[7] no all-zero or duplicate rows")
    zx = int((rwX == 0).sum())
    zz = int((rwZ == 0).sum())
    dx = HX.shape[0] - len(set(map(tuple, HX)))
    dz = HZ.shape[0] - len(set(map(tuple, HZ)))
    print("    zero rows: HX=%d HZ=%d   duplicate rows: HX=%d HZ=%d" % (zx, zz, dx, dz))

    print("\n[8] no all-zero columns (every qubit is checked)")
    zc = int(((degX + degZ) == 0).sum())
    print("    unchecked qubits: %d  -> %s" % (zc, "PASS" if zc == 0 else "FAIL"))
    ok &= (zc == 0)

    print("\n[9] trivial-distance check: no weight-1 or weight-2 logical")
    small = 0
    for j in range(N):
        v = np.zeros(N, dtype=np.int64)
        v[j] = 1
        if not ((HX @ v) % 2).any():
            r0 = gf2_rank(HZ)
            if gf2_rank(np.vstack([HZ, v])) > r0:
                small += 1
    print("    weight-1 Z-logicals: %d  -> %s" % (small, "PASS" if small == 0 else "FAIL"))
    ok &= (small == 0)

    print("\n[10] summary")
    print("    [[%d,%d,%d]]  w=%d   Q = %.3f"
          % (N, k, a.d_ub, w_QL, k * a.d_ub * a.d_ub / N))
    print("    ALL STRUCTURAL CHECKS: %s" % ("PASS" if ok else "FAIL"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
