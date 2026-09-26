Read("CLP_core.g");;  Read("pairs_mid.g");;  P1 := PAIRS;;
Read("pairs_bigm.g");;  PAIRS := Concatenation( P1, PAIRS );;
pr := First( PAIRS, p -> p.l = 84 and p.i = 8 and p.m = 42 and p.lift = 14 );;
P := LoadPair( pr );;  PrepareActions( P );;  Print( P.desc, "\n" );
for t in [ 1 .. 400 ] do
A := [ [ SupportWithOne( P.nA, 3 ) ] ];;
B := List( [1..3], ii -> List( [1..3], jj -> [ 1 ] ) );;
for ii in [2,3] do for jj in [2,3] do B[ii][jj] := [ Random( [ 2 .. P.nL ] ) ]; od; od;
if ForAny( Flat(A), x -> not P.AinN[x] ) then
c := CLP_Build( P, A, B, [1,1], [3,3] );;
if not HasCancellation( c, A, B, [1,1], [3,3], P.m ) then
K := c.N - RankMat(c.HX) - RankMat(c.HZ);;
if K >= 12 then
dZ := DistRandCSS( c.HX, c.HZ, 300, 0 : field := GF(2) );;
dX := DistRandCSS( c.HZ, c.HX, 300, 0 : field := GF(2) );;
Print( "K=", K, "  dZ<=", dZ, "  dX<=", dX, "\n" );
fi; fi; fi;
od;
QUIT;
