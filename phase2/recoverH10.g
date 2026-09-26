Read("CLP_core.g");;
Read("pairs_bigm.g");;

readmat := function( f, N )
    local s, d, i, M;
    s := StringFile(f);
    if s = fail then return fail; fi;
    d := Filtered(s, c -> c='0' or c='1');  M := [];
    for i in [0..Length(d)/N-1] do
        Add(M, List(d{[i*N+1..i*N+N]}, c -> (IntChar(c)-48)*One(GF(2))));
    od;
    return M;
end;;

# exact: which subset of MATS sums to BLK, by GF(2) linear algebra
MatchExact := function( blk, mats )
    local vecs, tgt, sol, hits, i;
    vecs := List( mats, M -> Flat(M) );
    tgt := Flat( blk );
    sol := SolutionMat( vecs, tgt );
    if sol = fail then return fail; fi;
    hits := Filtered( [1..Length(sol)], i -> not IsZero(sol[i]) );
    if Sum( List( hits, i -> mats[i] ) ) <> blk then return fail; fi;
    return hits;
end;;

pr := First( PAIRS, p -> p.l=270 and p.i=26 and p.m=90 );;
if pr = fail then Error("pair (270,26) m=90 not found"); fi;
P := LoadPair( pr );;  PrepareActions( P );;
m := P.m;;  off := 2*m;;

HX := readmat( "H10_n360k20d20_HX.txt", 360 );;
HZ := readmat( "H10_n360k20d20_HZ.txt", 360 );;

Print("========== H10  [[360,20,20]]  w=6  EXACT ==========\n\n");
Print("-- group --\n");
Print("G = SmallGroup(270,26)\n");
Print("  StructureDescription : ", StructureDescription(P.G), "\n");
Print("  |G| = ", Size(P.G), ",  abelian = ", IsAbelian(P.G),
      ",  nilpotent = ", IsNilpotent(P.G), ",  solvable = ", IsSolvable(P.G), "\n");
Print("  exponent = ", Exponent(P.G), ",  derived length = ", DerivedLength(P.G), "\n");
Print("  direct factors : ", List( DirectFactorsOfGroup(P.G), StructureDescription ), "\n");
Print("  centre = ", StructureDescription(Centre(P.G)), ", |Z| = ", Size(Centre(P.G)), "\n");
Print("  derived subgroup = ", StructureDescription(DerivedSubgroup(P.G)), "\n");
Print("H = ", StructureDescription(P.H), ",  |H| = ", Size(P.H),
      ",  normal = ", IsNormal(P.G,P.H), "\n");
Print("N_G(H) = ", StructureDescription(P.N), ",  |N| = ", Size(P.N), "\n");
Print("Core_G(H) = ", StructureDescription(P.core), ",  |Core| = ", Size(P.core), "\n");
Print("m = [G:H] = ", P.m, ",  r = [G:N] = ", P.r,
      ",  lift = [N:H] = ", P.lift, ",  stab = [H:Core] = ", P.stab, "\n");
Print("|Im L| = ", Size(P.G)/Size(P.core), ",  |Im R| = ", P.lift, "\n\n");

blk := function( M, br, bc ) return M{[(br-1)*m+1..br*m]}{[(bc-1)*m+1..bc*m]}; end;;
Ai  := MatchExact( blk(HX,1,1), P.Lmats );;
B11 := MatchExact( blk(HX,1,3), P.Rmats );;
B12 := MatchExact( blk(HX,1,4), P.Rmats );;
B21 := MatchExact( blk(HX,2,3), P.Rmats );;
B22 := MatchExact( blk(HX,2,4), P.Rmats );;

Print("-- group algebra elements --\n");
if fail in [Ai,B11,B12,B21,B22] then
    Print("FAILED: Ai=",Ai," B11=",B11," B12=",B12," B21=",B21," B22=",B22,"\n");
else
    Print("a  = ", List(Ai, i -> P.Areps[i]), "\n");
    Print("     orders ", List(Ai, i -> Order(P.Areps[i])), "\n");
    Print("     generates ", StructureDescription(Subgroup(P.G, List(Ai,i->P.Areps[i]))),
          "  (index ", Index(P.G, Subgroup(P.G, List(Ai,i->P.Areps[i]))), " in G)\n");
    Print("b11 = ", List(B11, i -> P.Lreps[i]), "   orders ", List(B11, i -> Order(P.Lreps[i])), "\n");
    Print("b12 = ", List(B12, i -> P.Lreps[i]), "   orders ", List(B12, i -> Order(P.Lreps[i])), "\n");
    Print("b21 = ", List(B21, i -> P.Lreps[i]), "   orders ", List(B21, i -> Order(P.Lreps[i])), "\n");
    Print("b22 = ", List(B22, i -> P.Lreps[i]), "   orders ", List(B22, i -> Order(P.Lreps[i])), "\n");
    Print("\n-- checks --\n");
    Print("some a-term outside N_G(H) (genuinely coset-based): ",
          ForAny(Ai, i -> not (P.Areps[i] in P.N)), "\n");
    Print("all b-terms in N_G(H): ",
          ForAll(Concatenation(B11,B12,B21,B22), i -> P.Lreps[i] in P.N), "\n");
fi;
Print("\n-- code --\n");
Print("n = 360,  k = ", 360 - RankMat(HX) - RankMat(HZ),
      ",  CSS ok = ", IsZero(HX*TransposedMat(HZ)), "\n");
QUIT;
