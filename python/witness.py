#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""usage: python witness.py H10_n360k20d20 360 20"""
import sys, json
import numpy as np
import gurobipy as gp
from gurobipy import GRB
exec(open("certify.py").read().split("def log(")[0].split('"""')[2])

cid, N, d = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
HX = parse_matrix(cid + "_HX.txt", N)
HZ = parse_matrix(cid + "_HZ.txt", N)
out = {}
for tag, H, Ho in (("dZ", HX, HZ), ("dX", HZ, HX)):
    L = logical_basis(H, Ho)
    found = False
    for i, lg in enumerate(L):
        m = gp.Model()
        m.Params.OutputFlag = 0
        m.Params.Threads = 4
        m.Params.BestObjStop = d + 0.5
        e = m.addMVar(N, vtype=GRB.BINARY)
        z = m.addMVar(H.shape[0], vtype=GRB.INTEGER, lb=0, ub=(H.sum(axis=1) // 2))
        w = m.addVar(vtype=GRB.INTEGER, lb=0, ub=int(lg.sum()) // 2)
        m.addConstr(H @ e == 2 * z)
        m.addConstr(lg @ e == 1 + 2 * w)
        m.setObjective(e.sum(), GRB.MINIMIZE)
        m.optimize()
        if m.SolCount > 0:
            v = np.round(e.X).astype(int)
            if int(v.sum()) == d:
                assert not ((H @ v) % 2).any(), "not in kernel"
                assert int(lg @ v) % 2 == 1, "does not anticommute"
                fn = "%s_witness_%s_sector%d_w%d.txt" % (cid, tag, i + 1, d)
                np.savetxt(fn, v, fmt="%d")
                out[tag] = {"sector": i + 1, "weight": d,
                            "support": np.nonzero(v)[0].tolist(), "file": fn}
                print("%s: weight-%d logical in sector %d" % (tag, d, i + 1))
                print("   support %s" % np.nonzero(v)[0].tolist())
                m.dispose()
                found = True
                break
        m.dispose()
    if not found:
        print("%s: NO logical of weight %d found" % (tag, d))
json.dump(out, open("%s_witnesses.json" % cid, "w"), indent=1)
print("wrote %s_witnesses.json" % cid)
