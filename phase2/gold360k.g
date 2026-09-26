Read("CLP_core.g");;
Read("pairs_bigm.g");;
SEED := 2402;;  Reset( GlobalMersenneTwister, SEED );;
RECORD := "gold360k_results.txt";;
PrintTo( RECORD, "# k>=12 then d; log only if Q beats frontier\n" );;
COUNT := 0;;

# frontier: w=6 -> 15.3, w=7 -> 18.0.  d capped at 22 (certifiable).
Frontier := function( w )
    if w = 6 then return 15.3; elif w = 7 then return 18.0; else return 1000.0; fi;
end;;

sel := Filtered( PAIRS, p -> p.l > 0 and p.m = 90 and p.r = 2 );;
Print( "pairs: ", Length(sel), "\n" );

for pr in sel do
  P := LoadPair( pr );;  PrepareActions( P );;
  nK := 0;;  nD := 0;;  bestK := 0;;  bestQ := 0.0;;
  Print( "== ", P.desc, " m=", P.m, " lift=", P.lift, " stab=", P.stab,
         " n=", 4*P.m, "\n" );
  for t in [ 1 .. 20000 ] do
    A := [ [ SupportWithOne( P.nA, 3 ) ] ];;
    B := [ [ SupportWithOne(P.nL,Random([1,2])), SupportWithOne(P.nL,Random([1,2])) ],
           [ SupportWithOne(P.nL,Random([1,2])), RandomEntry(P.nL,1,2) ] ];;
    if ForAny( Flat(A), x -> not P.AinN[x] ) then
      c := CLP_Build( P, A, B, [1,1], [2,2] );;
      if not HasCancellation( c, A, B, [1,1], [2,2], P.m ) then
        K := c.N - RankMat(c.HX) - RankMat(c.HZ);;
        if K > bestK then bestK := K; fi;
        # ---- GATE 1: dimension ----
        if K >= 12 then
          nK := nK + 1;
          wts := CodeWeights( c.HX, c.HZ );;
          if wts.wmax <= 7 then
            # smallest d that would beat the frontier, capped at 22
            need := 99;;
            for dd in [ 10 .. 22 ] do
              if Float( K*dd^2/c.N ) > Frontier( wts.wmax ) then need := dd; break; fi;
            od;
            # ---- GATE 2: distance, only for candidates that could win ----
            if need <= 22 then
              nD := nD + 1;
              d := DistUB( c.HX, c.HZ, 500, need - 1 );;
              if d >= need then
                d := DistUB( c.HX, c.HZ, 20000, need - 1 );;
                if d >= need then
                  COUNT := COUNT + 1;
                  bestQ := Float( K*d^2/c.N );
                  Print( "\n>>> GOLD [[", c.N, ",", K, ",<=", d, "]] w=", wts.wmax,
                         " Q=", bestQ, "  ", P.desc, "\n" );
                  AppendTo( RECORD, "\n#### G", COUNT, "  [[", c.N, ",", K, ",<=", d,
                    "]]  w=", wts.wmax, "  Q=", bestQ, "\n",
                    "group := \"", P.desc, "\";  l := ", P.l, ";  i := ", P.i,
                    ";  m := ", P.m, ";  lift := ", P.lift, ";  stab := ", P.stab, ";\n",
                    "Aidx := ", A, ";  Bidx := ", B, ";\n" );
                  WriteBinaryMatrix( Concatenation( "H", String(COUNT), "_n", String(c.N),
                     "k", String(K), "d", String(d), "_HX.txt" ), c.HX );
                  WriteBinaryMatrix( Concatenation( "H", String(COUNT), "_n", String(c.N),
                     "k", String(K), "d", String(d), "_HZ.txt" ), c.HZ );
                fi;
              fi;
            fi;
          fi;
        fi;
      fi;
    fi;
  od;
  Print( "   bestK=", bestK, "  K>=12: ", nK, "  reached distance stage: ", nD,
         "  logged so far: ", COUNT, "\n" );
od;
Print( "\ndone, ", COUNT, " codes\n" );
QUIT;
