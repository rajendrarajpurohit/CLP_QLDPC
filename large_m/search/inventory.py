import re, glob
txt = ""
for f in glob.glob("W6*_results.txt") + glob.glob("mut960_s*_results.txt") + glob.glob("shape_m216/*_results.txt"):
    txt += open(f).read().replace("\n", " ") + "\n"
recs = {}
for m in re.finditer(r"#### (\S+)\s+\[\[(\d+),(\d+),<=(\d+)\]\]\s+w=(\d+)\s+dZ<=\s*(\d+)\s+dX<=\s*(\d+)\s+Q=([\d.]+)\s+group := \"([^\"]+)\"", txt):
    cid, n, k, d, w, dz, dx, q, g = m.groups()
    n, k, d, w, dz, dx = map(int, (n, k, d, w, dz, dx))
    key = (n, k)
    if key not in recs or d > recs[key][0]:
        recs[key] = (d, w, dz, dx, g, cid)
rows = sorted([(n, k) + v for (n, k), v in recs.items()])
def dominated(r):
    n, k, d = r[0], r[1], r[2]
    return any((o[0] <= n and o[1] >= k and o[2] >= d and o != r) for o in rows)
print(f"{'n':>5} {'k':>4} {'d<=':>4} {'w':>2} {'dZ':>4} {'dX':>4} {'Q':>6}  pareto  group / id")
for r in rows:
    n, k, d, w, dz, dx, g, cid = r
    print(f"{n:>5} {k:>4} {d:>4} {w:>2} {dz:>4} {dx:>4} {k*d*d/n:>6.1f}  {'*' if not dominated(r) else ' '}       {g}  {cid}")
print(len(rows), "distinct (n,k), all bounds are 200000 round unless rescreened")
