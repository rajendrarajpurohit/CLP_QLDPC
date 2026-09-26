Read("CLP_core.g");;
SEED := 9001;;  Reset( GlobalMersenneTwister, SEED );;
OUT := "hope_scan.txt";;
PrintTo( OUT, "# a0 a1 b0 b1 alpha beta l i m n k w d Q\n" );;

kfloor := function( n )
    if n < 150 then return 6; elif n < 250 then return 8; else return 10; fi;
end;;

entry := function( nreps, alpha )
    if alpha = 1 then return [ Random([1..nreps]) ]; fi;
    return RandomSubset( [1..nreps], alpha );
end;;

genarr := function( nreps, r, c, alpha )
    local M, i, j;
    M := [];
    for i in [1..r] do
        M[i] := [];
        for j in [1..c] do M[i][j] := entry( nreps, alpha ); od;
    od;
    M[1][1] := Concatenation( [1], Difference( M[1][1], [1] ) );
    return M;
end;;

SHAPES := [];;
for a0 in [1,2,3] do for a1 in [1,2,3] do
for b0 in [1,2,3] do for b1 in [1,2,3] do
  if not ( a0=1 and a1=1 and b0=1 and b1=1 ) then
    for al in [1..5] do for be in [1..5] do
      if a0*al >= 3 and b1*be >= 3 then
        Add( SHAPES, [a0,a1,b0,b1,al,be] );
      fi;
    od; od;
  fi;
od; od; od; od;
Print( "shape/entry combinations: ", Length(SHAPES), "\n" );

Read("pairs_final.g");;  P1 := PAIRS;;
Read("pairs_mid.g");;    P2 := PAIRS;;
Read("pairs_bigm.g");;   P3 := PAIRS;;
ALL := Concatenation( Filtered(P1,p->p.l>0), Filtered(P2,p->p.l>0), Filtered(P3,p->p.l>0) );;
Print( "total pairs available: ", Length(ALL), "\n" );

for sh in SHAPES do
  a0 := sh[1];; a1 := sh[2];; b0 := sh[3];; b1 := sh[4];; al := sh[5];; be := sh[6];;
  npm := a1*b0 + a0*b1;;
  cand := Filtered( ALL, p -> p.m*npm >= 60 and p.m*npm <= 400 );;
  if Length(cand) = 0 then continue; fi;
  if Length(cand) > 6 then cand := RandomSubset( cand, 6 ); fi;
  Print( "== ", sh, "  n/m=", npm, "  pairs=", Length(cand), "\n" );
  for pr in cand do
    P := LoadPair( pr );;  PrepareActions( P );;
    n := P.m * npm;;
    kf := kfloor( n );;
    best := 0.0;;
    for t in [1..300] do
      A := genarr( P.nA, a0, a1, al );;
      B := genarr( P.nL, b0, b1, be );;
      if ForAny( Flat(A), x -> not P.AinN[x] ) then
        c := CLP_Build( P, A, B, [a0,a1], [b0,b1] );;
        K := c.N - RankMat(c.HX) - RankMat(c.HZ);;
        if K >= kf then
          wts := CodeWeights( c.HX, c.HZ );;
          if wts.wmax >= 5 and wts.wmax <= 10 then
            d := AbsInt( DistRandCSS( c.HX, c.HZ, 500, 2 : field := GF(2) ) );;
            q := Float( K*d^2/c.N );;
            if q > best then best := q; fi;
            AppendTo( OUT, a0," ",a1," ",b0," ",b1," ",al," ",be," ",
                      pr.l," ",pr.i," ",P.m," ",c.N," ",K," ",wts.wmax," ",d," ",q,"\n" );
          fi;
        fi;
      fi;
    od;
    Print( "   SG(", pr.l, ",", pr.i, ") m=", P.m, " n=", n, " bestQ=", best, "\n" );
  od;
od;
Print( "hope scan done\n" );
QUIT;
