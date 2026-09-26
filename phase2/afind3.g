Read("CLP_core.g");;
Read("pairs_bigm.g");;
readmat := function( f, N )
    local s, d, i, M;
    s := StringFile(f);  d := Filtered(s, c -> c='0' or c='1');  M := [];
    for i in [0..Length(d)/N-1] do
        Add(M, List(d{[i*N+1..i*N+N]}, c -> (IntChar(c)-48)*One(GF(2))));
    od;
    return M;
end;;
pr := First( PAIRS, p -> p.l=270 and p.i=26 and p.m=90 );;
P := LoadPair( pr );;  PrepareActions( P );;
m := P.m;;
HX := readmat( "H10_n360k20d20_HX.txt", 360 );;
tgt := HX{[1..m]}{[1..m]};;
# every element of G, not the transversal
all := AsList(P.G);;
Print("testing all ", Length(all), " elements of G\n");
mats := List( all, g -> PermutationMat( Image( ActionHomomorphism(P.G,
          RightCosets(P.G,P.H), OnRight), g^-1 )^-1, m, GF(2) ) );;
cand := Filtered( [1..Length(all)],
          i -> Position(mats[i][1], One(GF(2))) in Positions(tgt[1], One(GF(2))) );;
Print("candidates: ", Length(cand), "\n");
sol := fail;;
for s in Combinations( cand, 3 ) do
  if Sum( List(s, i -> mats[i]) ) = tgt then sol := s; break; fi;
od;
if sol = fail then
  Print("FAILED\n");
else
  Print("\na elements = ", List(sol, i -> all[i]), "\n");
  Print("a orders = ", List(sol, i -> Order(all[i])), "\n");
  Print("verify (1,1): ", Sum(List(sol,i->mats[i])) = tgt, "\n");
  Print("verify (2,2): ", Sum(List(sol,i->mats[i])) = HX{[m+1..2*m]}{[m+1..2*m]}, "\n");
  Print("some outside N: ", ForAny(sol, i -> not (all[i] in P.N)), "\n");
fi;
QUIT;
