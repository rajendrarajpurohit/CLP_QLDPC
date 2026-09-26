Read("CLP_core.g");;
Read("CLP_patch.g");;
Read("pairs_bigm.g");;
SEED := 701;;  Reset( GlobalMersenneTwister, SEED );;
OUT_PREFIX := "W7_";;  RECORD_FILE := "W7_results.txt";;  PrintTo( RECORD_FILE, "" );;
CODE_COUNTER := 0;;  SEEN_PARAMS := [];;
FULL_ROUNDS := 200000;;  CAND_PER_PAIR := 60000;;  RATE_FLOOR := 0.035;;
T := rec( label := "W7  1x1/2x2 large-m", sA := [1,1], sB := [2,2], wMax := 7, Kmin := 12,
          Nfrom := m -> 4 * m,
          pairOK := pr -> pr.m >= 60 and pr.m <= 100 and 4*pr.m <= 400 and pr.r = 2,
          genA := P -> [ [ SupportWithOne( P.nA, 3 ) ] ],
          genB := P -> [ [ SupportWithOne( P.nL, Random([1,2]) ), SupportWithOne( P.nL, Random([1,2]) ) ],
                         [ SupportWithOne( P.nL, Random([1,2]) ), RandomEntry( P.nL, 1, 2 ) ] ] );
RunDriver2( T );
QUIT;
