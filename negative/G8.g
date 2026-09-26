Read("CLP_core.g");;
Read("pairs_mid.g");;  P1 := PAIRS;;
Read("pairs_bigm.g");;  PAIRS := Concatenation( P1, PAIRS );;
SEED := 2608;;  Reset( GlobalMersenneTwister, SEED );;
RECORD := "G8_results.txt";;  NT := 20000;;  COUNT := 0;;
AppendTo( RECORD, "# G8 shape 3x3/1x1 n=6m  both sides screened  d=min(dZ,dX)\n" );;
Frontier := function( w )
  if w = 6 then return 15.3; elif w = 7 then return 18.0; elif w = 8 then return 30.6;
  elif w = 9 then return 20.3; elif w = 10 then return 38.8; else return 1000.0; fi;
end;;
NeedD := function( K, N, w )
  local dd;
  for dd in [ 8 .. 20 ] do if Float( K*dd^2/N ) > Frontier( w ) then return dd; fi; od;
  return 99;
end;;
sel := Filtered( PAIRS, p -> p.l > 0 and p.m >= 28 and p.m <= 64 and p.lift >= 3 );;
Print( "G8: pairs ", Length(sel), "\n" );
for pr in sel do
  P := LoadPair( pr );;  PrepareActions( P );;
  nK := 0;;  nD := 0;;  bestK := 0;;
  Print( "== ", P.desc, " m=", P.m, " lift=", P.lift, " stab=", P.stab, " n=", 6*P.m, "\n" );
  for t in [ 1 .. NT ] do
    A := List( [1..3], ii -> List( [1..3], jj -> [ 1 ] ) );;
    for ii in [2,3] do for jj in [2,3] do A[ii][jj] := [ Random( [ 2 .. P.nA ] ) ]; od; od;
    B := [ [ SupportWithOne( P.nL, 3 ) ] ];;
    if ForAny( Flat(A), x -> not P.AinN[x] ) then
      c := CLP_Build( P, A, B, [3,3], [1,1] );;
      if not HasCancellation( c, A, B, [3,3], [1,1], P.m ) then
        K := c.N - RankMat(c.HX) - RankMat(c.HZ);;
        if K > bestK then bestK := K; fi;
        if K >= 8 then
          nK := nK + 1;
          wts := CodeWeights( c.HX, c.HZ );;
          if wts.wmax <= 8 then
            need := NeedD( K, c.N, wts.wmax );;
            if need <= 20 then
              nD := nD + 1;
              if DistUB( c.HX, c.HZ, 50, 4 ) > 0 then
                d := DistUB( c.HX, c.HZ, 500, need-1 );;
                if d >= need then d := DistUB( c.HX, c.HZ, 20000, need-1 ); fi;
                if d >= need then
                  dZ := DistRandCSS( c.HX, c.HZ, 200000, need-1 : field := GF(2) );;
                  dX := DistRandCSS( c.HZ, c.HX, 200000, need-1 : field := GF(2) );;
                  if dZ >= need and dX >= need then
                    COUNT := COUNT + 1;  d := Minimum( dZ, dX );;
                    id := Concatenation( "G8_", String(COUNT), "_n", String(c.N), "k", String(K), "d", String(d), "w", String(wts.wmax) );;
                    Print( "\n>>> [[", c.N, ",", K, ",<=", d, "]] w=", wts.wmax, " dZ<=", dZ, " dX<=", dX, " Q=", Float(K*d^2/c.N), "  ", P.desc, "\n" );
                    AppendTo( RECORD, "\n#### ", id, "  [[", c.N, ",", K, ",<=", d, "]]  w=", wts.wmax, "  dZ<=", dZ, "  dX<=", dX, "  Q=", Float(K*d^2/c.N), "\n",
                      "group := \"", P.desc, "\";  l := ", P.l, ";  i := ", P.i, ";  m := ", P.m, ";  lift := ", P.lift, ";  stab := ", P.stab, ";\n",
                      "Hgens := ", pr.Hgens, ";\n",
                      "A := ", List(A, r -> List(r, e -> List(e, x -> P.Areps[x]))), ";\n",
                      "B := ", List(B, r -> List(r, e -> List(e, x -> P.Lreps[x]))), ";\n",
                      "Aidx := ", A, ";  Bidx := ", B, ";\n" );
                    WriteBinaryMatrix( Concatenation( id, "_HX.txt" ), c.HX );
                    WriteBinaryMatrix( Concatenation( id, "_HZ.txt" ), c.HZ );
                  fi;
                fi;
              fi;
            fi;
          fi;
        fi;
      fi;
    fi;
  od;
  Print( "   bestK=", bestK, "  K>=8: ", nK, "  distance stage: ", nD, "  logged: ", COUNT, "\n" );
od;
Print( "\nG8 done, ", COUNT, " codes\n" );
QUIT;
