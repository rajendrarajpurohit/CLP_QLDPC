Read("CLP_core.g");;
Read("pairs_bigm.g");;
pr := First( PAIRS, p -> p.l=270 and p.i=25 and p.m=90 );;
P := LoadPair(pr);;  PrepareActions(P);;
wh := [];;  kh := [];;
for t in [1..300] do
  A := [ [ SupportWithOne( P.nA, 4 ) ] ];;
  B := [ [ SupportWithOne(P.nL,Random([1,2])), SupportWithOne(P.nL,Random([1,2])) ],
         [ SupportWithOne(P.nL,Random([1,2])), RandomEntry(P.nL,1,2) ] ];;
  if ForAny( Flat(A), x -> not P.AinN[x] ) then
    c := CLP_Build( P, A, B, [1,1], [2,2] );;
    if not HasCancellation( c, A, B, [1,1], [2,2], P.m ) then
      K := c.N - RankMat(c.HX) - RankMat(c.HZ);;
      if K >= 12 then
        Add( wh, CodeWeights(c.HX,c.HZ).wmax );  Add( kh, K );
      fi;
    fi;
  fi;
od;
Print("samples with K>=12: ", Length(wh), "\n");
Print("weight histogram: ", Collected(wh), "\n");
Print("K histogram: ", Collected(kh), "\n");
QUIT;
