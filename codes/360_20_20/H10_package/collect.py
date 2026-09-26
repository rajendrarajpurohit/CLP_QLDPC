#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""usage: python collect.py H10_n360k20d20 20"""
import sys, glob, json, os
cid = sys.argv[1]
k = int(sys.argv[2])
res = {}
src = {}
for f in sorted(glob.glob(cid + "_progress.json") + glob.glob(cid + "_d*_*.json")):
    try:
        j = json.load(open(f))
    except Exception:
        continue
    if "dZ" in j and isinstance(j.get("dZ"), list):
        for i, v in enumerate(j["dZ"]):
            res[("dZ", i + 1)] = v;  src[("dZ", i + 1)] = f
        for i, v in enumerate(j.get("dX", [])):
            res[("dX", i + 1)] = v;  src[("dX", i + 1)] = f
    else:
        side = "dZ" if "_dZ_" in f else "dX"
        for i, r in j.items():
            mw = r.get("min_weight")
            val = int(str(mw).lstrip(">=")) if mw else None
            res[(side, int(i))] = val;  src[(side, int(i))] = f
print("%-4s %-7s %-8s %s" % ("side", "sector", "result", "source file"))
missing = []
for side in ("dZ", "dX"):
    for i in range(1, k + 1):
        if (side, i) in res:
            print("%-4s %-7d >=%-6s %s" % (side, i, res[(side, i)], src[(side, i)]))
        else:
            print("%-4s %-7d %-8s -" % (side, i, "MISSING"))
            missing.append((side, i))
done = len(res)
print("\n%d of %d sectors proven.  missing: %s" % (done, 2 * k, missing if missing else "none"))
if not missing:
    print("d = %d  ->  EXACT" % min(res.values()))
