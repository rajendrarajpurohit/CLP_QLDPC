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
sel := Filtered( PAIRS, p -> p.l=270 and p.i=26 and p.m=90 );;
pr := sel[3];;
P := LoadPair( pr );;  PrepareActions( P );;
m := P.m;;
HX := rdm( "H10_n360k20d20_HX.txt", 360 );;
HZ := rdm( "H10_n360k20d20_HZ.txt", 360 );;
bl := function( M, br, bc ) return M{[(br-1)*m+1..br*m]}{[(bc-1)*m+1..bc*m]}; end;;
match := function( tgt, mats, nmax )
    local supp, cand, t, s;
    supp := Positions( tgt[1], One(GF(2)) );
    cand := Filtered( [1..Length(mats)],
              i -> Position(mats[i][1],One(GF(2))) in supp );
    for t in [1..nmax] do
      for s in Combinations( cand, t ) do
        if Sum( List(s, i -> mats[i]) ) = tgt then return s; fi;
      od;
    od;
    return fail;
end;;
Ai  := match( bl(HX,1,1), P.Lmats, 3 );;
B11 := match( bl(HX,1,3), P.Rmats, 2 );;
B12 := match( bl(HX,1,4), P.Rmats, 2 );;
B21 := match( bl(HX,2,3), P.Rmats, 2 );;
B22 := match( bl(HX,2,4), P.Rmats, 2 );;
Print("========== H10  [[360,20,20]]  w=6  EXACT ==========\n\n");
Print("G = SmallGroup(270,26) = ", StructureDescription(P.G), "\n");
Print("  direct factors: ", List(DirectFactorsOfGroup(P.G), StructureDescription), "\n");
Print("  centre ", StructureDescription(Centre(P.G)),
      ", derived ", StructureDescription(DerivedSubgroup(P.G)), "\n");
Print("H = C3, Hgens exponents = ", pr.Hgens, ", normal = ", IsNormal(P.G,P.H), "\n");
Print("N_G(H) = ", StructureDescription(P.N), ", Core = ", Size(P.core), "\n");
Print("m=", P.m, " r=", P.r, " lift=", P.lift, " stab=", P.stab, "\n\n");
Print("a   = ", List(Ai,  i -> P.Areps[i]), "   orders ", List(Ai, i->Order(P.Areps[i])), "\n");
Print("b11 = ", List(B11, i -> P.Lreps[i]), "   orders ", List(B11,i->Order(P.Lreps[i])), "\n");
Print("b12 = ", List(B12, i -> P.Lreps[i]), "   orders ", List(B12,i->Order(P.Lreps[i])), "\n");
Print("b21 = ", List(B21, i -> P.Lreps[i]), "   orders ", List(B21,i->Order(P.Lreps[i])), "\n");
Print("b22 = ", List(B22, i -> P.Lreps[i]), "   orders ", List(B22,i->Order(P.Lreps[i])), "\n\n");
Print("Aidx = ", [[Ai]], "\n");
Print("Bidx = ", [[B11,B12],[B21,B22]], "\n\n");
Print("some a outside N_G(H): ", ForAny(Ai, i -> not (P.Areps[i] in P.N)), "\n");
Print("all b in N_G(H): ",
      ForAll(Concatenation(B11,B12,B21,B22), i -> P.Lreps[i] in P.N), "\n");
Print("rebuild check: ");
c := CLP_Build( P, [[Ai]], [[B11,B12],[B21,B22]], [1,1], [2,2] );;
Print("HX identical = ", c.HX = HX, ",  HZ identical = ", c.HZ = HZ, "\n");
Print("k = ", c.N - RankMat(c.HX) - RankMat(c.HZ),
      ",  w = ", CodeWeights(c.HX,c.HZ).wmax, "\n");
QUIT;
