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
blk := function( M, br, bc ) return M{[(br-1)*m+1..br*m]}{[(bc-1)*m+1..bc*m]}; end;;
Print("nA = ", P.nA, "   nL = ", P.nL, "\n");
for i in [1..2] do
  for j in [1..4] do
    Print("HX block(",i,",",j,") rowwt = ",
          Number( blk(HX,i,j)[1], y -> not IsZero(y) ), "\n");
  od;
od;
Print("\nA from results file, Aidx = [1,32,262]:\n");
cand := P.Lmats[1] + P.Lmats[32] + P.Lmats[262];;
Print("  matches block(1,1): ", cand = blk(HX,1,1), "\n");
Print("  matches block(2,2): ", cand = blk(HX,2,2), "\n");
Print("  elements: ", [P.Areps[1],P.Areps[32],P.Areps[262]], "\n");
Print("  orders:   ", List([1,32,262], i -> Order(P.Areps[i])), "\n");
QUIT;
