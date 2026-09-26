Read("CLP_core.g");;
Read("pairs_bigm.g");;
rdm := function( f, N )
    local s, d, i, M;
    s := StringFile(f);  d := Filtered(s, c -> c='0' or c='1');  M := [];
    for i in [0..Length(d)/N-1] do
        Add(M, List(d{[i*N+1..i*N+N]}, c -> (IntChar(c)-48)*One(GF(2))));
    od;
    return M;
end;;
pr := First( PAIRS, p -> p.l=270 and p.i=26 and p.m=90 );;
P := LoadPair( pr );;
cos := RightCosets( P.G, P.H );;
hom := ActionHomomorphism( P.G, cos, OnRight );;
m := P.m;;
Lm := g -> PermutationMat( Image( hom, g^-1 )^-1, m, GF(2) );;
HX := rdm( "H10_n360k20d20_HX.txt", 360 );;
tgt := HX{[1..m]}{[1..m]};;
els := AsList( P.G );;
mats := List( els, Lm );;
supp1 := Positions( tgt[1], One(GF(2)) );;
Print("row-1 support of target: ", supp1, "\n");
cand := Filtered( [1..Length(els)], i -> Position(mats[i][1],One(GF(2))) in supp1 );;
Print("candidates: ", Length(cand), "\n");
sol := fail;;
for s in Combinations( cand, 3 ) do
  if Sum( List(s, i -> mats[i]) ) = tgt then sol := s; break; fi;
od;
if sol <> fail then
  Print("\nFOUND a = ", List(sol, i -> els[i]), "\n");
  Print("orders ", List(sol, i -> Order(els[i])), "\n");
  Print("identity present: ", ForAny(sol, i -> els[i] = One(P.G)), "\n");
  Print("verify (1,1): ", Sum(List(sol,i->mats[i])) = tgt, "\n");
  Print("verify (2,2): ", Sum(List(sol,i->mats[i])) = HX{[m+1..2*m]}{[m+1..2*m]}, "\n");
  Print("some outside N: ", ForAny(sol, i -> not (els[i] in P.N)), "\n");
else
  Print("\nno 3-subset works.  checking whether tgt is a sum of ANY L-matrices:\n");
  Print("  is tgt a permutation matrix? ",
        ForAll(tgt, r -> Number(r, y -> not IsZero(y)) = 1), "\n");
  Print("  row weights present: ",
        Set(List(tgt, r -> Number(r, y -> not IsZero(y)))), "\n");
  Print("  col weights present: ",
        Set(List(TransposedMat(tgt), c -> Number(c, y -> not IsZero(y)))), "\n");
fi;
QUIT;
