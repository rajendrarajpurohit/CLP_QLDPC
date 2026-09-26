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
res := tgt - P.Lmats[1];;
Print("searching pairs from ", P.nA, " elements...\n");
found := fail;;
for i in [2..P.nA] do
  for j in [i+1..P.nA] do
    if P.Lmats[i] + P.Lmats[j] = res then
      found := [1,i,j];  break;
    fi;
  od;
  if found <> fail then break; fi;
od;
if found = fail then
  Print("no pair found; trying all triples without assuming identity\n");
else
  Print("\na indices (this session) = ", found, "\n");
  Print("a elements = ", List(found, i -> P.Areps[i]), "\n");
  Print("a orders   = ", List(found, i -> Order(P.Areps[i])), "\n");
  Print("verify block(1,1): ", Sum(List(found, i -> P.Lmats[i])) = tgt, "\n");
  Print("verify block(2,2): ",
        Sum(List(found, i -> P.Lmats[i])) = HX{[m+1..2*m]}{[m+1..2*m]}, "\n");
  Print("some a outside N: ", ForAny(found, i -> not (P.Areps[i] in P.N)), "\n");
  Print("a generates ", StructureDescription(Subgroup(P.G, List(found, i -> P.Areps[i]))),
        " (index ", Index(P.G, Subgroup(P.G, List(found, i -> P.Areps[i]))), ")\n");
fi;
QUIT;
