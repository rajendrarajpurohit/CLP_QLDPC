cd C:\Users\rajen\Documents\release
@'
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