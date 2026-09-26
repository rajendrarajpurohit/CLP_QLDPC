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
HX := rdm( "H10_n360k20d20_HX.txt", 360 );;
sel := Filtered( PAIRS, p -> p.l=270 and p.i=26 and p.m=90 );;
Print("trying ", Length(sel), " candidate pairs\n");
for idx in [1..Length(sel)] do
  P := LoadPair( sel[idx] );;  PrepareActions( P );;
  m := P.m;;
  tgt := HX{[1..m]}{[1..m]};;
  supp1 := Positions( tgt[1], One(GF(2)) );;
  cand := Filtered( [1..P.nA],
            i -> Position(P.Lmats[i][1],One(GF(2))) in supp1 );;
  sol := fail;;
  for s in Combinations( cand, 3 ) do
    if Sum( List(s, i -> P.Lmats[i]) ) = tgt then sol := s; break; fi;
  od;
  Print("pair ", idx, " (Hgens=", sel[idx].Hgens, "): ");
  if sol = fail then
    Print("no\n");
  else
    Print("MATCH  a = ", List(sol, i -> P.Areps[i]),
          "  orders ", List(sol, i -> Order(P.Areps[i])), "\n");
    Print("   verify (2,2): ",
          Sum(List(sol,i->P.Lmats[i])) = HX{[m+1..2*m]}{[m+1..2*m]}, "\n");
    Print("   b11 = ", Positions(HX{[1..m]}{[2*m+1..3*m]}[1], One(GF(2))), "\n");
  fi;
od;
QUIT;
