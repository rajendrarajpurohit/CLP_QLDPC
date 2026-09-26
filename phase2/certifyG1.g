Read("pairs_bigm.g");;

readmat := function( f, N )
    local s, d, i, M;
    s := StringFile(f);
    d := Filtered(s, c -> c='0' or c='1');
    M := [];
    for i in [0..Length(d)/N-1] do
        Add(M, List(d{[i*N+1..i*N+N]}, c -> (IntChar(c)-48)*One(GF(2))));
    od;
    return M;
end;;

pr := First( PAIRS, p -> p.l=270 and p.i=25 and p.m=90 );;
grp := SmallGroup(270,25);;
sub := Subgroup( grp, List( pr.Hgens, v -> PcElementByExponents( Pcgs(grp), v ) ) );;
nor := Normalizer(grp,sub);;
cor := Core(grp,sub);;
cos := RightCosets(grp,sub);;
m := Length(cos);;
hom := ActionHomomorphism( grp, cos, OnRight );;

Print("=== STEP A: group data ===\n");
Print("G = ", StructureDescription(grp), ", |G| = ", Size(grp), "\n");
Print("H = ", StructureDescription(sub), ", |H| = ", Size(sub), "\n");
Print("H normal in G: ", IsNormal(grp,sub), "   (must be false)\n");
Print("N_G(H) = ", StructureDescription(nor), ", |N| = ", Size(nor), "\n");
Print("Core_G(H) = ", StructureDescription(cor), ", |Core| = ", Size(cor), "\n");
Print("m = [G:H] = ", m, ", r = [G:N] = ", Index(grp,nor), ", lift = [N:H] = ", Index(nor,sub), "\n");

Lm := g -> PermutationMat( Image(hom,g^-1)^-1, m, GF(2) );;
Rm := function( u )
    local p;
    p := PermList( List( cos, c -> Position( cos, RightCoset( sub, u^-1 * Representative(c) ) ) ) );
    return PermutationMat( p^-1, m, GF(2) );
end;;

gg := GeneratorsOfGroup(grp);;
f1 := gg[1];; f2 := gg[2];; f3 := gg[3];; f4 := gg[4];; f5 := gg[5];;
at := [ One(grp), f2^2*f3^2*f5, f1*f2*f5^3 ];;
b11 := [ One(grp), f4^2*f5^4 ];;
b12 := [ One(grp), f2^2*f5 ];;
b21 := [ One(grp), f2*f4^2*f5^3 ];;
b22 := [ f2*f5^3 ];;
ball := Concatenation(b11,b12,b21,b22);;

Print("\n=== STEP B: the elements ===\n");
Print("a = ", at, "\n    orders ", List(at,Order), "\n");
Print("all a in G: ", ForAll(at, x -> x in grp), "\n");
Print("some a outside N (else it is 2BGA): ", ForAny(at, x -> not (x in nor)), "\n");
Print("all b in N_G(H): ", ForAll(ball, x -> x in nor), "\n");

Print("\n=== STEP C: representation properties ===\n");
tst := List([1..6], i -> Random(grp));;
Print("L is a homomorphism: ", ForAll(tst, x -> ForAll(tst, y -> Lm(x*y) = Lm(x)*Lm(y))), "\n");
tsn := List([1..6], i -> Random(nor));;
Print("R is an anti-homomorphism: ", ForAll(tsn, x -> ForAll(tsn, y -> Rm(x*y) = Rm(y)*Rm(x))), "\n");
Print("L(e) = I: ", Lm(One(grp)) = IdentityMat(m,GF(2)), "\n");
Print("R(e) = I: ", Rm(One(grp)) = IdentityMat(m,GF(2)), "\n");
Print("ker L = Core (|Im L| = |G/Core|): ", Length(Set(List(AsList(grp),Lm))) = Size(grp)/Size(cor), "\n");
Print("ker R = H (|Im R| = |N/H|): ", Length(Set(List(AsList(nor),Rm))) = Index(nor,sub), "\n");
Print("L(g)R(u) = R(u)L(g) for all g,u: ", ForAll(tst, x -> ForAll(tsn, y -> Lm(x)*Rm(y) = Rm(y)*Lm(x))), "\n");

amat := Sum(List(at,Lm));;
m11 := Sum(List(b11,Rm));;
m12 := Sum(List(b12,Rm));;
m21 := Sum(List(b21,Rm));;
m22 := Sum(List(b22,Rm));;

HXs := readmat("G1_n360k16d22_HX.txt",360);;
HZs := readmat("G1_n360k16d22_HZ.txt",360);;

Print("\n=== STEP D: do the stored blocks equal the algebra? ===\n");
blk := function( M, br, bc ) return M{[(br-1)*m+1..br*m]}{[(bc-1)*m+1..bc*m]}; end;;
Print("HX is 2x4 blocks of size ", m, "\n");
for i in [1..2] do
  for j in [1..4] do
    Print("  HXblock(",i,",",j,") = ",
      ["A","B11","B12","B21","B22","ZERO","other"][
        PositionProperty([amat,m11,m12,m21,m22,NullMat(m,m,GF(2))],
                         X -> X = blk(HXs,i,j)) ], "\n");
  od;
od;
for i in [1..2] do
  for j in [1..4] do
    Print("  HZblock(",i,",",j,") = ",
      ["A^T","B11^T","B12^T","B21^T","B22^T","ZERO","other"][
        PositionProperty([TransposedMat(amat),TransposedMat(m11),TransposedMat(m12),
                          TransposedMat(m21),TransposedMat(m22),NullMat(m,m,GF(2))],
                         X -> X = blk(HZs,i,j)) ], "\n");
  od;
od;

Print("\n=== STEP E: code properties from the stored matrix ===\n");
Print("CSS: HX.HZ^T = 0: ", IsZero(HXs*TransposedMat(HZs)), "\n");
Print("k = ", 360 - RankMat(HXs) - RankMat(HZs), "\n");
Print("max HX row wt = ", Maximum(List(HXs, r -> Number(r, y -> not IsZero(y)))), "\n");
Print("max HZ row wt = ", Maximum(List(HZs, r -> Number(r, y -> not IsZero(y)))), "\n");
QUIT;
