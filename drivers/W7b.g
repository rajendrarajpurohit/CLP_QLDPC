W := 7;;  MLO := 56;;  MHI := 100;;
Read("CLP_core.g");;
Read("pairs_final.g");;  P1 := PAIRS;;
Read("pairs_mid.g");;    P2 := PAIRS;;
Read("pairs_bigm.g");;   PAIRS := Concatenation( P1, P2, PAIRS );;
SEED := 2700 + W + 100*MLO;;  Reset( GlobalMersenneTwister, SEED );;
TAG := Concatenation( "W", String(W), "_m", String(MLO), "-", String(MHI) );;
RECORD := Concatenation( TAG, "_results.txt" );;  NT := 5000;;  COUNT := 0;;
AppendTo( RECORD, "# ", TAG, " shape 1x1/2x2 n=4m  w<=", W, "  both sides  d=min(dZ,dX)\n" );;
Frontier := function( w )
  if w = 6 then return 15.3; elif w = 7 then return 18.0; elif w = 8 then return 30.6;
  elif w = 9 then return 20.3; elif w = 10 then return 38.8; else return 1000.0; fi;
end;;
NeedD := function( K, N, w )
  local dd;
  for dd in [ 8 .. 28 ] do if Float( K*dd^2/N ) > Frontier( w ) then return dd; fi; od;
  return 99;
end;;
Bsum := function( B )
  return Maximum( Concatenation( List( B, r -> Sum( List( r, Length ) ) ),
                                 List( [1,2], j -> Sum( List( B, r -> Length( r[j] ) ) ) ) ) );
end;;
sel := Filtered( PAIRS, p -> p.l > 0 and p.m >= MLO and p.m <= MHI and p.lift >= 4 );;
Print( TAG, ": pairs ", Length(sel), "\n" );
for pr in sel do
  P := LoadPair( pr );;  PrepareActions( P );;
  nK := 0;;  nD := 0;;  bestK := 0;;
  Print( "== ", P.desc, " m=", P.m, " lift=", P.lift, " stab=", P.stab, " n=", 4*P.m, "\n" );
  for t in [ 1 .. NT ] do
    wa := Random( [ 3 .. W-3 ] );;
    A := [ [ SupportWithOne( P.nA, wa ) ] ];;
    B := [ [ SupportWithOne(P.nL,Random([1,2,3])), SupportWithOne(P.nL,Random([1,2,3])) ],
           [ SupportWithOne(P.nL,Random([1,2,3])), RandomEntry(P.nL,1,3) ] ];;
    if Bsum( B ) <= W - wa and ForAny( Flat(A), x -> not P.AinN[x] ) then
      c := CLP_Build( P, A, B, [1,1], [2,2] );;
      if not IsZero( c.HX * TransposedMat( c.HZ ) ) then Error( "not CSS" ); fi;
      if not HasCancellation( c, A, B, [1,1], [2,2], P.m ) then
        K := c.N - RankMat(c.HX) - RankMat(c.HZ);;
        if K > bestK then bestK := K; fi;
        if K >= 8 then
          nK := nK + 1;
          wts := CodeWeights( c.HX, c.HZ );;
          if wts.wmax <= W then
            need := NeedD( K, c.N, wts.wmax );;
            if need <= 28 then
              nD := nD + 1;
              if DistUB( c.HX, c.HZ, 50, 4 ) > 0 then
                d := DistUB( c.HX, c.HZ, 500, need-1 );;
                if d >= need then d := DistUB( c.HX, c.HZ, 20000, need-1 ); fi;
                if d >= need then
                  dZ := DistRandCSS( c.HX, c.HZ, 200000, need-1 : field := GF(2) );;
                  dX := DistRandCSS( c.HZ, c.HX, 200000, need-1 : field := GF(2) );;
                  if dZ >= need and dX >= need then
                    COUNT := COUNT + 1;  d := Minimum( dZ, dX );;
                    id := Concatenation( TAG, "_", String(COUNT), "_n", String(c.N), "k", String(K), "d", String(d), "w", String(wts.wmax) );;
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
Print( "\n", TAG, " done, ", COUNT, " codes\n" );
QUIT;
