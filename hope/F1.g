Read("CLP_core.g");;
SEED := 1101;;  Reset( GlobalMersenneTwister, SEED );;
OUT := "F1_results.txt";;
PrintTo( OUT, "# l i m n k w d Q alpha beta shape\n" );;
NCAND := 600;;  NPAIR := 16;;

dupcols := function( M )
    local T;
    T := TransposedMat( M );
    return Length( Set( T ) ) < Length( T );
end;;

thincol := function( HX, HZ )
    local cx, cz, i;
    cx := List( TransposedMat(HX), c -> Number(c, y -> not IsZero(y)) );
    cz := List( TransposedMat(HZ), c -> Number(c, y -> not IsZero(y)) );
    return ForAny( [1..Length(cx)], i -> cx[i]+cz[i] < 2 );
end;;

kfloor := function( w )
    if w = 5 then return 6; elif w = 6 then return 10; elif w = 7 then return 12;
    elif w = 8 then return 16; elif w = 9 then return 12; else return 20; fi;
end;;

qbar := function( w )
    if w = 5 then return 6.0; elif w = 6 then return 12.0; elif w = 7 then return 14.0;
    elif w = 8 then return 20.0; elif w = 9 then return 16.0; else return 25.0; fi;
end;;

entry := function( nreps, a )
    if a = 1 then return [ Random([1..nreps]) ]; fi;
    return RandomSubset( [1..nreps], a );
end;;

Read("pairs_final.g");;  Q1 := Filtered(PAIRS, p->p.l>0);;
Read("pairs_mid.g");;    Q2 := Filtered(PAIRS, p->p.l>0);;
Read("pairs_bigm.g");;   Q3 := Filtered(PAIRS, p->p.l>0);;
ALL := Concatenation( Q1, Q2, Q3 );;
Print( "pairs available: ", Length(ALL), "\n" );

NPM := 4;;  A0 := 1;;  A1 := 1;;  B0 := 2;;  B1 := 2;;
ALPHA := [3];;  BETA := [1,2];;
COUNT := 0;;

for mm in [ 15 .. 96 ] do
  sel := Filtered( ALL, p -> p.m = mm and mm*NPM >= 60 and mm*NPM <= 400 );
  if Length(sel) = 0 then continue; fi;
  if Length(sel) > NPAIR then sel := RandomSubset( sel, NPAIR ); fi;
  Print( "== m=", mm, " n=", mm*NPM, " pairs=", Length(sel), "\n" );
  for pr in sel do
    P := LoadPair( pr );;  PrepareActions( P );;
    best := 0.0;;
    for t in [1..NCAND] do
      al := Random( ALPHA );;  be := Random( BETA );;
      A := [ [ Concatenation([1], Difference(entry(P.nA,al),[1])) ] ];;
      B := [ [ entry(P.nL,be), entry(P.nL,be) ],
             [ entry(P.nL,be), entry(P.nL,be) ] ];;
      B[1][1] := Concatenation([1], Difference(B[1][1],[1]));;
      if ForAny( Flat(A), x -> not P.AinN[x] ) then
        c := CLP_Build( P, A, B, [A0,A1], [B0,B1] );;
        if not dupcols(c.HX) and not dupcols(c.HZ) and not thincol(c.HX,c.HZ) then
          wts := CodeWeights( c.HX, c.HZ );;
          if wts.wmax >= 5 and wts.wmax <= 10 then
            K := c.N - RankMat(c.HX) - RankMat(c.HZ);;
            if K >= kfloor( wts.wmax ) then
              d := AbsInt( DistRandCSS( c.HX, c.HZ, 500, 2 : field := GF(2) ) );;
              q := Float( K*d^2/c.N );;
              if q > best then best := q; fi;
              if q > qbar( wts.wmax ) then
                COUNT := COUNT + 1;
                AppendTo( OUT, pr.l," ",pr.i," ",P.m," ",c.N," ",K," ",wts.wmax,
                          " ",d," ",q," ",al," ",be," 1x1/2x2\n" );
                Print( "  >>> [[", c.N, ",", K, ",<=", d, "]] w=", wts.wmax,
                       " Q=", q, "  SG(", pr.l, ",", pr.i, ")\n" );
                WriteBinaryMatrix( Concatenation("F1_",String(COUNT),"_n",String(c.N),
                  "k",String(K),"d",String(d),"w",String(wts.wmax),"_HX.txt"), c.HX );
                WriteBinaryMatrix( Concatenation("F1_",String(COUNT),"_n",String(c.N),
                  "k",String(K),"d",String(d),"w",String(wts.wmax),"_HZ.txt"), c.HZ );
                AppendTo( OUT, "#   A := ", List(A, r -> List(r, e -> List(e, x -> P.Areps[x]))),
                          "\n#   B := ", List(B, r -> List(r, e -> List(e, x -> P.Lreps[x]))),
                          "\n#   Aidx := ", A, "  Bidx := ", B, "\n" );
              fi;
            fi;
          fi;
        fi;
      fi;
    od;
    Print( "   SG(", pr.l, ",", pr.i, ") bestQ=", best, "\n" );
  od;
od;
Print( "\nF1 done, ", COUNT, " codes\n" );
QUIT;
