Read("CLP_core.g");;
Read("CLP_patch.g");;
Read("pairs_W10.g");;
SEED := 901;;  Reset( GlobalMersenneTwister, SEED );;
OUT_PREFIX := "W9_";;  RECORD_FILE := "W9_results.txt";;  PrintTo( RECORD_FILE, "" );;
CODE_COUNTER := 0;;  SEEN_PARAMS := [];;
FULL_ROUNDS := 200000;;  CAND_PER_PAIR := 40000;;  RATE_FLOOR := 0.06;;
T := rec( label := "W9  3x4/5x3", sA := [3,4], sB := [5,3], wMax := 9, Kmin := 1,
          Nfrom := m -> 29 * m,
          pairOK := pr -> pr.m = 12,
          genA := P -> GaugeBase( P.nA, 3, 4 ),
          genB := P -> GaugeBase( P.nL, 5, 3 ) );
RunDriver2( T );
QUIT;
