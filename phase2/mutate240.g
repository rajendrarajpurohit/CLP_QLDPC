Read("CLP_core.g");;
Read("pairs.g");;
SEED := 8020;;  Reset( GlobalMersenneTwister, SEED );;

pr := First( PAIRS, p -> p.l = 360 and p.i = 99 and p.m = 60 );;
if pr = fail then Error("pair SG(360,99) m=60 not found in pairs.g"); fi;
P := LoadPair( pr );;  PrepareActions( P );;
Print( P.desc, " m=", P.m, " nA=", P.nA, " nL=", P.nL, "\n" );

A0 := [ [ [ 1, 114, 116 ] ] ];;
B0 := [ [ [ 1 ], [ 1, 26 ] ], [ [ 1, 29 ], [ 28 ] ] ];;
c := CLP_Build( P, A0, B0, [1,1], [2,2] );;
Print( "PARENT: n=", c.N, " K=", c.N - RankMat(c.HX) - RankMat(c.HZ),
       " w=", CodeWeights(c.HX,c.HZ).wmax, "  (expect K=8 w=6)\n" );

PrintTo( "mutate240_results.txt", "# neighbours of a240_8_20\n" );;
best := 8;;  COUNT := 0;;
for t in [ 1 .. 300000 ] do
  A := StructuralCopy( A0 );;  B := StructuralCopy( B0 );;
  for r in [ 1 .. Random([1,2]) ] do
    if Random([1,2]) = 1 then
      A[1][1] := SupportWithOne( P.nA, 3 );
    else
      i := Random([1,2]);;  j := Random([1,2]);;
      if i = 1 and j = 1 then B[1][1] := SupportWithOne( P.nL, Random([1,2]) );
      elif i = 1 then B[1][2] := SupportWithOne( P.nL, Random([1,2]) );
      elif j = 1 then B[2][1] := SupportWithOne( P.nL, Random([1,2]) );
      else B[2][2] := RandomEntry( P.nL, 1, 2 ); fi;
    fi;
  od;
  if ForAny( Flat(A), x -> not P.AinN[x] ) then
    c := CLP_Build( P, A, B, [1,1], [2,2] );;
    if not HasCancellation( c, A, B, [1,1], [2,2], P.m ) then
      K := c.N - RankMat(c.HX) - RankMat(c.HZ);;
      if K > best then best := K; Print( "   new bestK=", K, " at t=", t, "\n" ); fi;
      if K >= 10 then
        wts := CodeWeights( c.HX, c.HZ );;
        if wts.wmax <= 7 then
          need := 15.3;;  if wts.wmax = 7 then need := 18.0; fi;
          d := AbsInt( DistUB( c.HX, c.HZ, 2000, 2 ) );;
          if Float( K*d^2/c.N ) > need then
            COUNT := COUNT + 1;
            Print( ">>> [[", c.N, ",", K, ",<=", d, "]] w=", wts.wmax,
                   " Q=", Float(K*d^2/c.N), "\n" );
            AppendTo( "mutate240_results.txt", "[[", c.N, ",", K, ",<=", d,
              "]] w=", wts.wmax, " Q=", Float(K*d^2/c.N),
              "  Aidx:=", A, " Bidx:=", B, "\n" );
            WriteBinaryMatrix( Concatenation("M",String(COUNT),"_k",String(K),
              "d",String(d),"_HX.txt"), c.HX );
            WriteBinaryMatrix( Concatenation("M",String(COUNT),"_k",String(K),
              "d",String(d),"_HZ.txt"), c.HZ );
          fi;
        fi;
      fi;
    fi;
  fi;
  if t mod 20000 = 0 then
    Print( "  ", t, " tried, bestK=", best, ", logged=", COUNT, "\n" );
  fi;
od;
Print( "done: ", COUNT, " codes\n" );
QUIT;
