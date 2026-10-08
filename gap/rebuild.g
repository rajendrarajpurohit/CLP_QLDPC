# rebuild.g   build a CLP code from the description printed in the paper
#
#   Rebuild( l, i, Hw, Aw, Bw, sA, sB, name )
#     l, i    SmallGroup(l,i)
#     Hw      f -> list of generators of H as words in the pc generators f[1..]
#     Aw      f -> A as a list of rows, each entry a list of group elements (the terms of the sum)
#     Bw      f -> B in the same form
#     sA, sB  block shapes, [1,1] and [2,2] for the 1x1 over 2x2 layout, [2,2] and [1,1] for 2x2 over 1x1
#     name    output prefix, writes name_HX.txt and name_HZ.txt
#   then   python3 verify.py name N d   and   python3 certify.py name N d
Read( "CLP_core.g" );
Rebuild := function( l, i, Hw, Aw, Bw, sA, sB, name )
    local G, f, H, pr, P, idxA, idxB, A, B, code, k, w;
    G := SmallGroup( l, i );  f := Pcgs( G );
    H := Subgroup( G, Hw( f ) );
    pr := rec( l := l, i := i, pc := true, Hdesc := StructureDescription( H ),
               Hgens := List( GeneratorsOfGroup( H ), h -> ExponentsOfPcElement( f, h ) ) );
    P := LoadPair( pr );  PrepareActions( P );
    f := Pcgs( P.G );
    idxA := g -> First( [ 1 .. P.nA ], k -> P.Areps[k]^-1 * g in P.core );
    idxB := u -> First( [ 1 .. P.nL ], k -> P.Lreps[k]^-1 * u in P.H );
    A := List( Aw( f ), row -> List( row, ent -> List( ent, idxA ) ) );
    B := List( Bw( f ), row -> List( row, ent -> List( ent, idxB ) ) );
    if fail in Flat( A ) then Error( "an entry of A is not in G/Core" ); fi;
    if fail in Flat( B ) then Error( "an entry of B is not in N_G(H)" ); fi;
    code := CLP_Build( P, A, B, sA, sB );
    k := code.N - RankMat( code.HX ) - RankMat( code.HZ );
    w := CodeWeights( code.HX, code.HZ );
    Print( name, "  ", P.desc, "  m=", P.m, "  [[", code.N, ",", k, "]]  check weight ", w.wmax,
           "  QDistRnd 20000 rounds d<=", DistUB( code.HX, code.HZ, 20000, 1 ), "\n" );
    WriteBinaryMatrix( Concatenation( name, "_HX.txt" ), code.HX );
    WriteBinaryMatrix( Concatenation( name, "_HZ.txt" ), code.HZ );
end;
