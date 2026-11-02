import re, glob, json, os, shutil
BIG = os.path.expanduser("~/Rajendra/Big"); QC = os.path.expanduser("~/Rajendra/qldpc-challenge")
recs = ""
for f in glob.glob(BIG+"/W6*_results.txt") + glob.glob(BIG+"/mut960_s*_results.txt"):
    recs += open(f).read().replace("\n", " ") + "\n"
res = open(BIG+"/screen1e6/rescreen_results.txt").read() + open(BIG+"/rescreen/rescreen_results.txt").read()
ids = ["W6focus_c7_1_n576k16d30w6","W6focus_c16_1_n576k20d24w6","W6focus_c10_2_n588k18d30w6","W6focus_c7_2_n640k16d34w6",
       "W6focus_c4_1_n640k24d24w6","W6focus_c10_4_n672k24d26w6","W6rest_c20_1_n720k16d40w6","W6rest_c11_1_n840k16d42w6",
       "W6rest_c18_1_n864k16d48w6","W6rest_c11_2_n864k20d42w6","W6rest_c9_2_n960k20d58w6","mut960_s21_2_n960k20d64w6",
       "W6rest_c11_3_n960k16d54w6","W6rest_c29_3_n960k16d54w6","W6rest_c22_1_n960k16d48w6","W6rest_c16_1_n960k16d62w6","W6rest_c20_2_n1000k16d62w6"]
def word(s): return s.replace("<identity> of ...","1").replace("*"," ")
table = []
for cid in ids:
    m = re.search(r"#### "+re.escape(cid)+r"\s+\[\[(\d+),(\d+),<=(\d+)\]\].*?Hgens := (\[.*?\]);\s*A := (\[.*?\]);\s*B := (\[.*?\]);\s*Aidx := (\[.*?\]);\s*Bidx :=\s*(\[.*?\]);", recs)
    if not m: print("no record", cid); continue
    n, k = int(m.group(1)), int(m.group(2))
    g = re.search(r"SmallGroup\((\d+),(\d+)\) H=(\w+)", recs[m.start():m.end()]); l, i, Ht = g.group(1), g.group(2), g.group(3)
    lines = [x for x in res.splitlines() if cid in x]
    if not lines: print("no screen", cid); continue
    dz = min(int(re.search(r"dZ<=(\d+)",x).group(1)) for x in lines); dx = min(int(re.search(r"dX<=(\d+)",x).group(1)) for x in lines)
    d = min(dz, dx); rounds = sum(int(re.search(r"rounds=(\d+)",x).group(1)) for x in lines)
    bw = None; dj = QC+f"/{cid}.deep.json"
    if os.path.exists(dj) and os.path.getsize(dj)>0:
        w = json.load(open(dj))
        if "doc" in w: bw = w["doc"]["distance"]; d = min(d, bw["X"]["value"], bw["Z"]["value"])
    folder = f"{n}_{k}_{d}"
    if os.path.exists(folder): folder = f"{n}_{k}_{d}_G{l}_{i}"
    os.makedirs(folder, exist_ok=True)
    for s in ("HX","HZ"):
        for src in (BIG+"/screen1e6/", BIG+"/rescreen/", BIG+"/"):
            if os.path.exists(src+cid+f"_{s}.txt"): shutil.copy(src+cid+f"_{s}.txt", folder+f"/{s}.txt"); break
    open(folder+"/RECORD.txt","w").write(recs[m.start():m.end()]+"\n")
    open(folder+"/screen.txt","w").write(f"# QDistRnd two sided screen, {len(lines)} run(s), total rounds {rounds}, min over runs dZ<={dz} dX<={dx}, d<={d}\n"+"\n".join(lines)+"\n")
    if bw: json.dump(bw, open(folder+"/board_witnesses.json","w"), indent=1); open(folder+"/screen.txt","a").write(f"# board tool RIS witness: dX<={bw['X']['value']} dZ<={bw['Z']['value']}, see board_witnesses.json\n")
    open(folder+"/REBUILD.txt","w").write(f"# G = SmallGroup({l},{i}), H = {Ht}, Hgens as exponent vectors in Pcgs(G): {m.group(4)}\n"
        f"# a = {word(m.group(5))}\n# B = {word(m.group(6))}\n"
        f"# index form for CLP_core.g: Aidx = {m.group(7)}  Bidx = {m.group(8)}\n"
        f"# rebuild: gap/rebuild.g, Rebuild({l},{i}, Hw, Aw, Bw, [1,1], [2,2], \"{folder}\") with Hw, Aw, Bw transcribed from the words above\n")
    table.append((n,k,d,k*d*d/n,l,i,Ht,rounds,cid,folder))
table.sort()
with open("TABLE.txt","w") as f:
    f.write(f"{'n':>5} {'k':>3} {'d<=':>4} {'Q':>6}  {'group':<20} {'H':<4} {'rounds':>9}  id\n")
    for n,k,d,q,l,i,Ht,r,cid,fo in table: f.write(f"{n:>5} {k:>3} {d:>4} {q:>6.1f}  SmallGroup({l},{i}){'':<{20-len(f'SmallGroup({l},{i})')}} {Ht:<4} {r:>9}  {fo:<22} {cid}\n")
print(open("TABLE.txt").read())
