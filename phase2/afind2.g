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
Print("target row 1 support: ", Positions(tgt[1], One(GF(2))), "\n");
# a permutation matrix belongs to the sum iff its row-1 entry lands in the support
cand := Filtered( [1..P.nA],
          i -> Position(P.Lmats[i][1], One(GF(2))) in Positions(tgt[1], One(GF(2))) );;
Print("candidates matching row 1: ", Length(cand), " -> ", cand, "\n");
# now test every 3-subset of those candidates
sol := fail;;
for s in Combinations( cand, 3 ) do
  if Sum( List(s, i -> P.Lmats[i]) ) = tgt then sol := s; break; fi;
od;
if sol = fail then
  Print("still not found -- printing first 3 rows of target for inspection\n");
  Print(Positions(tgt[1],One(GF(2))), "\n");
  Print(Positions(tgt[2],One(GF(2))), "\n");
  Print(Positions(tgt[3],One(GF(2))), "\n");
else
  Print("\na indices = ", sol, "\n");
  Print("a elements = ", List(sol, i -> P.Areps[i]), "\n");
  Print("a orders = ", List(sol, i -> Order(P.Areps[i])), "\n");
  Print("verify (1,1): ", Sum(List(sol,i->P.Lmats[i])) = tgt, "\n");
  Print("verify (2,2): ", Sum(List(sol,i->P.Lmats[i])) = HX{[m+1..2*m]}{[m+1..2*m]}, "\n");
  Print("some a outside N: ", ForAny(sol, i -> not (P.Areps[i] in P.N)), "\n");
fi;
QUIT;
