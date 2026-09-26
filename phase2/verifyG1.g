Read("pairs_bigm.g");;

CLPReadN := function( f, N )
    local s, d, i, M;
    s := StringFile(f);  d := Filtered(s, c -> c='0' or c='1');  M := [];
    for i in [0..Length(d)/N-1] do
        Add(M, List(d{[i*N+1..i*N+N]}, c -> (IntChar(c)-48)*One(GF(2))));
    od;
    return M;
end;;

pr := First( PAIRS, p -> p.l=270 and p.i=25 and p.m=90 );;
G := SmallGroup(270,25);;
H := Subgroup( G, List( pr.Hgens, v -> PcElementByExponents( Pcgs(G), v ) ) );;
N := Normalizer(G,H);;
C := Core(G,H);;
cos := RightCosets(G,H);;
m := Length(cos);;
Print("|G|=",Size(G)," |H|=",Size(H)," |N|=",Size(N)," |Core|=",Size(C)," m=",m,"\n");
Print("H normal: ", IsNormal(G,H), "\n");

Lmat := function( g )
    local M, j, t;
    M := NullMat(m,m,GF(2));
    for j in [1..m] do
        t := PositionProperty( cos, c -> g*Representative(cos[j]) in c );
        M[t][j] := One(GF(2));
    od;
    return M;
end;;

Rmat := function( u )
    local M, j, t;
    M := NullMat(m,m,GF(2));
    for j in [1..m] do
        t := PositionProperty( cos, c -> Representative(cos[j])*u in c );
        M[t][j] := One(GF(2));
    od;
    return M;
end;;

gg := GeneratorsOfGroup(G);;
f1 := gg[1];; f2 := gg[2];; f3 := gg[3];; f4 := gg[4];; f5 := gg[5];;
aterms := [ One(G), f2^2*f3^2*f5, f1*f2*f5^3 ];;
b11 := [ One(G), f4^2*f5^4 ];;
b12 := [ One(G), f2^2*f5 ];;
b21 := [ One(G), f2*f4^2*f5^3 ];;
b22 := [ f2*f5^3 ];;

Print("\n-- membership --\n");
Print("a in G: ", List(aterms, x -> x in G), "  orders: ", List(aterms,Order), "\n");
Print("some a outside N (needed): ", ForAny(aterms, x -> not (x in N)), "\n");
Print("b in N: ", List([b11,b12,b21,b22], e -> List(e, x -> x in N)), "\n");

La := List( aterms, Lmat );;
Ball := Concatenation( b11, b12, b21, b22 );;
Rb := List( Ball, Rmat );;
Print("\n-- commutation --\n");
Print("L(a) R(b) = R(b) L(a) for all pairs: ", ForAll(La, X -> ForAll(Rb, Y -> X*Y = Y*X)), "\n");
Print("L homomorphism: ", Lmat(aterms[2]*aterms[3]) = Lmat(aterms[2])*Lmat(aterms[3]), "\n");
Print("R anti-homomorphism: ", Rmat(b11[2]*b21[2]) = Rmat(b21[2])*Rmat(b11[2]), "\n");

A := Sum( La );;
B11 := Sum( List(b11,Rmat) );;
B12 := Sum( List(b12,Rmat) );;
B21 := Sum( List(b21,Rmat) );;
B22 := Sum( List(b22,Rmat) );;
nullblk := NullMat(m,m,GF(2));;
tA := TransposedMat(A);;
t11 := TransposedMat(B11);;
t12 := TransposedMat(B12);;
t21 := TransposedMat(B21);;
t22 := TransposedMat(B22);;

RowX1 := r -> Concatenation( A[r], nullblk[r], B11[r], B12[r] );;
RowX2 := r -> Concatenation( nullblk[r], A[r], B21[r], B22[r] );;
RowZ1 := r -> Concatenation( t11[r], t21[r], tA[r], nullblk[r] );;
RowZ2 := r -> Concatenation( t12[r], t22[r], nullblk[r], tA[r] );;
HXn := Concatenation( List([1..m],RowX1), List([1..m],RowX2) );;
HZn := Concatenation( List([1..m],RowZ1), List([1..m],RowZ2) );;

Print("\n-- hand-built --\n");
Print("n=", Length(HXn[1]), " HXrows=", Length(HXn), " HZrows=", Length(HZn), "\n");
Print("CSS ok: ", IsZero( HXn*TransposedMat(HZn) ), "\n");
Print("k = ", Length(HXn[1]) - RankMat(HXn) - RankMat(HZn), "\n");
Print("max HX rowwt=", Maximum(List(HXn, r -> Number(r, y->not IsZero(y)))), "\n");
Print("max HZ rowwt=", Maximum(List(HZn, r -> Number(r, y->not IsZero(y)))), "\n");

HXs := CLPReadN("G1_n360k16d22_HX.txt", 360);;
HZs := CLPReadN("G1_n360k16d22_HZ.txt", 360);;
Print("\n-- stored --\n");
Print("stored CSS ok: ", IsZero(HXs*TransposedMat(HZs)), "\n");
Print("stored k = ", 360 - RankMat(HXs) - RankMat(HZs), "\n");
Print("HXn = HXs exactly: ", HXn = HXs, "\n");
rx := RankMat(HXs);;
rz := RankMat(HZs);;
Print("same X rowspace: ", RankMat(HXn) = rx and RankMat(Concatenation(HXn,HXs)) = rx, "\n");
Print("same Z rowspace: ", RankMat(HZn) = rz and RankMat(Concatenation(HZn,HZs)) = rz, "\n");
QUIT;
