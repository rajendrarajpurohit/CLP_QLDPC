Read("CLP_core.g");;  Read("CLP_patch.g");;  Read("pairs_W10.g");;
best := 0.0;;  n := 0;;
for pr in Filtered( PAIRS, p -> p.l > 0 and p.m = 12 ) do
  P := LoadPair( pr );;  PrepareActionsCapped( P, 3000 );;
  for t in [ 1 .. 500 ] do
    A := GaugeBase( P.nA, 3, 5 );;  B := GaugeBase( P.nL, 5, 3 );;
    if ForAny( Flat(A), x -> not P.AinN[x] ) then
      code := CLP_Build( P, A, B, [3,5], [5,3] );;
      if not HasCancellation( code, A, B, [3,5], [5,3], P.m ) then
        K := code.N - RankMat(code.HX) - RankMat(code.HZ);;
        if K >= 60 then
          d := AbsInt( DistRandCSS( code.HX, code.HZ, 2000, 2 : field := GF(2) ) );;
          n := n + 1;
          if Float(K*d^2/code.N) > best then
            best := Float(K*d^2/code.N);
            Print( "  best so far: [[", code.N, ",", K, ",<=", d, "]] Q=", best,
                   "  ", P.desc, "\n" );
          fi;
        fi;
      fi;
    fi;
  od;
  Print( P.desc, ": ", n, " samples, best Q=", best, "\n" );
od;
QUIT;
