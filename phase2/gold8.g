Read("CLP_core.g");;
Read("pairs_bigm.g");;
SEED := 8001;;  Reset( GlobalMersenneTwister, SEED );;
RECORD := "gold8_results.txt";;
PrintTo( RECORD, "# w=8: A has 4 terms.  log only if Q > 30.6\n" );;
COUNT := 0;;
sel := Filtered( PAIRS, p -> p.l > 0 and p.m = 90 and p.r = 2 );;
Print( "pairs: ", Length(sel), "\n" );
for pr in sel do
  P := LoadPair( pr );;  PrepareActions( P );;
  nK := 0;;  bestK := 0;;
  Print( "== ", P.desc, " m=", P.m, " lift=", P.lift, " n=", 4*P.m, "\n" );
  for t in [ 1 .. 20000 ] do
    A := [ [ SupportWithOne( P.nA, 4 ) ] ];;
    B := [ [ SupportWithOne(P.nL,Random([1,2])), SupportWithOne(P.nL,Random([1,2])) ],
           [ SupportWithOne(P.nL,Random([1,2])), RandomEntry(P.nL,1,2) ] ];;
    if ForAny( Flat(A), x -> not P.AinN[x] ) then
      c := CLP_Build( P, A, B, [1,1], [2,2] );;
      if not HasCancellation( c, A, B, [1,1], [2,2], P.m ) then
        K := c.N - RankMat(c.HX) - RankMat(c.HZ);;
        if K > bestK then bestK := K; fi;
        if K >= 12 then
          nK := nK + 1;
          wts := CodeWeights( c.HX, c.HZ );;
          if wts.wmax >= 8 and wts.wmax <= 10 then
            need := 99;;
            for dd in [ 12 .. 26 ] do
              if Float( K*dd^2/c.N ) > 30.6 then need := dd; break; fi;
            od;
            if need <= 26 then
              d := DistUB( c.HX, c.HZ, 500, need - 1 );;
              if d >= need then
                d := DistUB( c.HX, c.HZ, 20000, need - 1 );;
                if d >= need then
                  COUNT := COUNT + 1;
                  Print( "\n>>> W8GOLD [[", c.N, ",", K, ",<=", d, "]] w=", wts.wmax,
                         " Q=", Float(K*d^2/c.N), "  ", P.desc, "\n" );
                  AppendTo( RECORD, "\n#### W", COUNT, "  [[", c.N, ",", K, ",<=", d,
                    "]]  w=", wts.wmax, "  Q=", Float(K*d^2/c.N), "\n",
                    "group := \"", P.desc, "\";  l := ", P.l, ";  i := ", P.i,
                    ";  m := ", P.m, ";  lift := ", P.lift, ";  stab := ", P.stab, ";\n",
                    "A := ", List(A, row -> List(row, e -> List(e, k -> P.Areps[k]))), "\n",
                    "B := ", List(B, row -> List(row, e -> List(e, k -> P.Lreps[k]))), "\n",
                    "Aidx := ", A, ";  Bidx := ", B, ";\n" );
                  WriteBinaryMatrix( Concatenation("W",String(COUNT),"_n",String(c.N),
                     "k",String(K),"d",String(d),"w",String(wts.wmax),"_HX.txt"), c.HX );
                  WriteBinaryMatrix( Concatenation("W",String(COUNT),"_n",String(c.N),
                     "k",String(K),"d",String(d),"w",String(wts.wmax),"_HZ.txt"), c.HZ );
                fi;
              fi;
            fi;
          fi;
        fi;
      fi;
    fi;
  od;
  Print( "   bestK=", bestK, "  K>=12: ", nK, "  logged: ", COUNT, "\n" );
od;
Print( "\ndone, ", COUNT, " codes\n" );
QUIT;
