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
