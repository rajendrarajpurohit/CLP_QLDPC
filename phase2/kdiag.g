Read("CLP_core.g");;
Read("pairs_bigm.g");;
sel := Filtered( PAIRS, p -> p.l > 0 and p.m in [60,72,84] and p.r = 2 );;
Print( "pairs found: ", Length(sel), "\n" );
if Length(sel) > 5 then sel := sel{[1..5]}; fi;
for pr in sel do
  P := LoadPair( pr );;
  PrepareActions( P );;
  Print( "== ", P.desc, " m=", P.m, " lift=", P.lift,
         " stab=", P.stab, " nA=", P.nA, " nL=", P.nL, "\n" );
  for t in [1..20] do
    A := [ [ SupportWithOne( P.nA, 3 ) ] ];;
    B := [ [ SupportWithOne(P.nL,2), SupportWithOne(P.nL,2) ],
           [ SupportWithOne(P.nL,2), RandomEntry(P.nL,1,2) ] ];;
    if ForAny( Flat(A), x -> not P.AinN[x] ) then
      c := CLP_Build( P, A, B, [1,1], [2,2] );;
      if not HasCancellation( c, A, B, [1,1], [2,2], P.m ) then
        Print( "   n=", c.N, "  K=", c.N - RankMat(c.HX) - RankMat(c.HZ), "\n" );
      fi;
    fi;
  od;
od;
QUIT;
