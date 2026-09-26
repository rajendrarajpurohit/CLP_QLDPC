Read("CLP_core.g");;
Read("CLP_patch.g");;
Read("pairs_W10.g");;
SEED := 1001;;  Reset( GlobalMersenneTwister, SEED );;
OUT_PREFIX := "W10_";;  RECORD_FILE := "W10_results.txt";;  PrintTo( RECORD_FILE, "" );;
CODE_COUNTER := 0;;  SEEN_PARAMS := [];;
FULL_ROUNDS := 200000;;  CAND_PER_PAIR := 40000;;  RATE_FLOOR := 0.10;;
T := rec( label := "W10  3x5/5x3", sA := [3,5], sB := [5,3], wMax := 10, Kmin := 1,
          Nfrom := m -> 34 * m,
          pairOK := pr -> pr.m = 12,
          genA := P -> GaugeBase( P.nA, 3, 5 ),
          genB := P -> GaugeBase( P.nL, 5, 3 ) );
RunDriver2( T );
QUIT;
