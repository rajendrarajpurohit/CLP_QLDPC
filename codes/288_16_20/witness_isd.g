Read("rescreen.g");;
FindWitness := function( Hcheck, Hstab, w, rounds, outfile )
    local n, K, r, perm, Kp, M, row, v;
    n := Length( Hcheck[1] );
    K := NullspaceMat( TransposedMat( Hcheck ) );
    ConvertToMatrixRep( K, GF(2) );
    for r in [ 1 .. rounds ] do
        perm := Random( SymmetricGroup( n ) );
        Kp := List( K, x -> Permuted( x, perm ) );
        M := TriangulizedMat( Kp );
        for row in M do
            if WeightVecFFE( row ) <= w then
                v := Permuted( row, perm^-1 );
                if SolutionMat( Hstab, v ) = fail then
                    PrintTo( outfile, Concatenation( List( v, x -> String( IntFFE( x ) ) ) ), "\n" );
                    Print( "found weight ", WeightVecFFE( v ), " at round ", r, "\n" );
                    return v;
                fi;
            fi;
        od;
        if r mod 20000 = 0 then Print( "round ", r, "\n" ); fi;
    od;
    Print( "not found in ", rounds, " rounds\n" );
    return fail;
end;;
id := "W10_m56-100_1_n288k16d20w9";;
HX := ReadBin( Concatenation( id, "_HX.txt" ) );;
HZ := ReadBin( Concatenation( id, "_HZ.txt" ) );;
ConvertToMatrixRep( HX, GF(2) );;  ConvertToMatrixRep( HZ, GF(2) );;
Print( "Z side\n" );  FindWitness( HX, HZ, 20, 400000, Concatenation( id, "_witness_Z_w20.txt" ) );
Print( "X side\n" );  FindWitness( HZ, HX, 20, 400000, Concatenation( id, "_witness_X_w20.txt" ) );
QUIT;
