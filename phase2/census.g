Read("CLP_core.g");;
Read("pairs_bigm.g");;
SEED := 5001;;  Reset( GlobalMersenneTwister, SEED );;
OUT := "census.txt";;
PrintTo( OUT, "# l i m n bestK nK6 nK7 hits\n" );;
sel := Filtered( PAIRS, p -> p.l > 0 and p.m >= 56 and p.m <= 96 and p.r = 2 );;
Print( "pairs: ", Length(sel), "\n" );
cnt := 0;;
for pr in sel do
  cnt := cnt + 1;
  P := LoadPair( pr );;  PrepareActions( P );;
  bestK := 0;;  nK6 := 0;;  nK7 := 0;;  hits := 0;;
  for t in [ 1 .. 500 ] do
    A := [ [ SupportWithOne( P.nA, 3 ) ] ];;
    B := [ [ SupportWithOne(P.nL,Random([1,2])), SupportWithOne(P.nL,Random([1,2])) ],
           [ SupportWithOne(P.nL,Random([1,2])), RandomEntry(P.nL,1,2) ] ];;
    if ForAny( Flat(A), x -> not P.AinN[x] ) then
      c := CLP_Build( P, A, B, [1,1], [2,2] );;
      if not HasCancellation( c, A, B, [1,1], [2,2], P.m ) then
        K := c.N - RankMat(c.HX) - RankMat(c.HZ);;
        if K > bestK then bestK := K; fi;
        if K >= 10 then
          wts := CodeWeights( c.HX, c.HZ );;
          if wts.wmax = 6 and Float( K*400/c.N ) > 15.3 then hits := hits + 1; fi;
          if wts.wmax = 7 and Float( K*400/c.N ) > 18.0 then hits := hits + 1; fi;
          if wts.wmax = 6 then nK6 := nK6 + 1; fi;
          if wts.wmax = 7 then nK7 := nK7 + 1; fi;
        fi;
      fi;
    fi;
  od;
  AppendTo( OUT, pr.l, " ", pr.i, " ", pr.m, " ", 4*pr.m, " ",
            bestK, " ", nK6, " ", nK7, " ", hits, "\n" );
  if cnt mod 25 = 0 then Print( cnt, "/", Length(sel), "\n" ); fi;
od;
Print( "census done\n" );
QUIT;
