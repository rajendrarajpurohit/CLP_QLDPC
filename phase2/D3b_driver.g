Read( "CLP_core.g" );;
Read( "pairs.g" );;
SEED := 304;;  Reset( GlobalMersenneTwister, SEED );;
OUT_PREFIX := "D3b_";;  RECORD_FILE := "D3b_results.txt";;
PrintTo( RECORD_FILE, "" );;
NMAX := 400;;  CAND_PER_PAIR := 60000;;
MIN_D_LOG := 14;;  FOM_FLOOR := 13.0;;
CODE_COUNTER := 0;;  SEEN_PARAMS := [];;

T := rec( label := "D3b  1x1/2x2 k>=10",
          sA := [1,1], sB := [2,2], wMax := 7, Kmin := 10,
          Nfrom := m -> 4 * m,
          pairOK := pr -> 4 * pr.m <= NMAX and 4 * pr.m >= 240,
          genA := P -> [ [ SupportWithOne( P.nA, 3 ) ] ],
          genB := P -> [ [ SupportWithOne( P.nL, Random( [1,2] ) ), SupportWithOne( P.nL, Random( [1,2] ) ) ],
                         [ SupportWithOne( P.nL, Random( [1,2] ) ), RandomEntry( P.nL, 1, 2 ) ] ] );
RunDriver( T );
QUIT;
