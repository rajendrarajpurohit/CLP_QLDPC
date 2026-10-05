# Where each code came from

Phase 1  ~/Rajendra/CLP        Sep 7 to 8    drivers D1 to D4 run inline, not saved
Phase 2  ~/Rajendra/Gold       Sep 9 to 14   wide search, records prefixed G in gold360k, called H here
Hope     ~/Rajendra/Hope       Sep 14 to 15  F1 to F5 and screen4
Phase 3  ~/Rajendra/Nova       Sep 15 to 16  layouts G1 to G9, G3 to G9 archived as negative
Phase 4  ~/Rajendra/Nova       Sep 16 on     W7 to W10 drivers, F9

Code               w   Phase  Record                          Results file
[[168,10,14]]      6   1      D3_6                            CLP/D3_results.txt
[[240,8,20]]       6   1      D3_4                            CLP/D3_results.txt
[[360,20,20]]      6   2      H10 (gold360k G10, line 69)     Gold/gold360k_results.txt
[[252,16,18]]      9   4      W10_m56-100_12                  Nova/W10_m56-100_results.txt
[[224,12,<=21]]    9   4      W9_m56-100_15                   Nova/W9_m56-100_results.txt
[[224,12,<=20]]    9   4      W9_m56-100_16                   Nova/W9_m56-100_results.txt
[[288,16,<=20]]    9   4      W10_m56-100_1                   Nova/W10_m56-100_results.txt
[[336,18,<=20]]    9   4      W9_m56-100_5                    Nova/W9_m56-100_results.txt
[[360,24,<=20]]    7   4      W7_m56-100_38                   Nova/W7_m56-100_results.txt
[[336,16,<=20]]    7   4      W7_m56-100_23                   Nova/W7_m56-100_results.txt
[[288,12,<=23]]    9   4      W9_m56-100_13                   Nova/W9_m56-100_results.txt
[[288,12,<=24]]    9   4      W9_m56-100_1                    Nova/W9_m56-100_results.txt
[[336,16,<=24]]    9   4      W9_m56-100_4                    Nova/W9_m56-100_results.txt
[[360,16,<=25]]    7   4      W7_m56-100_53                   Nova/W7_m56-100_results.txt
[[360,16,<=26]]    9   4      W10_m56-100_32                  Nova/W10_m56-100_results.txt
[[360,12,<=30]]    9   4      W9_m56-100_31                   Nova/W9_m56-100_results.txt
[[360,12,<=31]]    9   4      W9_m56-100_36                   Nova/W9_m56-100_results.txt
[[372,10,<=32]]    9   4      W9_m56-100_8                    Nova/W9_m56-100_results.txt
[[360,20,<=27]]   10   4      W10_m56-100_28                  Nova/W10_m56-100_results.txt
[[384,14,<=24]]    9   3      G2_13                           Nova/G2_results.txt
[[360,16,<=24]]    9   3      G1_94                           Nova/G1_results.txt
[[360,20,<=22]]    7   2      H75 (gold360k G75, line 459)    Gold/gold360k_results.txt
[[360,16,<=26]]    7   2      V287                            Gold/wide_results.txt
[[360,20,<=22]]    7   2      V177                            Gold/wide_results.txt

Certification trails
[[240,8,20]]    CLP/Rajendra/D3_4_dZ_*.log, D3_4_dZ_sectors_1-5.json, D3_4_dZ_witness_w20.npy
[[168,10,14]]   CLP/keep/a168_10_14_certify.log
[[360,20,20]]   CLP/keep/H10_*.log, h10*.log, H10_package/
[[252,16,18]]   Nova/frontier/252_package/
[[360,20,20]]   codes/360_20_20/  H10 logs and H10_package
[[252,16,18]]   codes/252_16_18/  252_package

Phase 1 layouts  D1 1x1/1x1 (ATB shape, validation only)  D2 2x2/1x1  D3 1x1/2x2  D4 2x2/2x2 (empty)
Note  Phase 1 matrix files use GAP backslash line continuation for n >= 216, unwrap before strict parsing

Phase 1 extras, certified 27 Sep 2026, codes/phase1_extras/<id>/, H pinned by gap/recoverD3.g
[[120,8,8]]     w=5  D2_15  SG(360,99)  H=<f2,f3f6,f5>=C6xC2
[[180,12,8]]    w=5  D2_23  SG(135,3)   H=<f2>=C3
[[240,16,8]]    w=5  D2_4   SG(360,99)  H=<f2,f3>=C6
[[112,12,10]]   w=7  D3_41  SG(336,188) H=<f6,f1,f2>=D12   (second class rebuilds identically)
[[168,14,12]]   w=7  D3_39  SG(336,188) H=<f1f3,f2,f5>=D8  (second class rebuilds identically)
[[168,20,10]]   w=7  D3_56  SG(336,188) H=<f1f2,f2f3,f5>=D8 (second class rebuilds identically)
[[180,16,12]]   w=7  D3_17  SG(135,3)   H=<f1f2>=C3
[[224,24,9]]    w=7  D3_5   SG(336,188) H=<f6,f1f2>=S3
[[240,16,13]]   w=7  D3_3   SG(360,99)  H=<f3f6,f5>=C6
[[252,22,10]]   w=7  D3_8   SG(126,8)   H=<f1>=C2

Not reported
[[224,12,<=21]] w=9  W9_m56-100_15   MILP started 25 Sep, stopped for memory, upper bound only
[[216,8,<=18]]  w=7  D2_10           rescreen held at 10^6, MILP queued
D1 validation codes [[56,6,8]] [[36,4,6]] and the gross code [[144,12,12]]  pipeline checks only

Queued when compute is free
[[360,24,<=20]] w=7  W7_m56-100_38   then  [[216,8,<=18]] w=7  then  [[288,16,<=20]] w=9
[[216,8,18]]    w=7  D2_10  SG(162,3)   H=<f3f5>=C3   certified 2 Oct 2026, moved from Not reported
[[360,24,20]]   w=7  W7_m56-100_38  SG(810,84) H=<f6,f2f4>=C9   certified 2 Oct 2026, codes/360_24_20/
[[216,8,18]]    w=7  D2_10  SG(162,3)   H=<f3f5>=C3   certified 2 Oct 2026, moved from Not reported
[[360,24,20]]   w=7  W7_m56-100_38  SG(810,84) H=<f6,f2f4>=C9   certified 2 Oct 2026, codes/360_24_20/
[[288,16,20]]   w=9  W10_m56-100_1  SG(144,164) H=<f1f5>=C2   certified 5 Oct 2026, codes/288_16_20/, w=9 leader Q=22.22
