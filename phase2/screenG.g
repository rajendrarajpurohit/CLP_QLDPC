if LoadPackage( "QDistRnd" ) = fail then Error("no QDistRnd"); fi;
RowWeights := function( M ) return List( M, r -> Number( r, y -> not IsZero(y) ) ); end;
CodeWeights := function( HX, HZ )
    local rwX, rwZ, dX, dZ;
    rwX := RowWeights(HX);  rwZ := RowWeights(HZ);
    dX := RowWeights(TransposedMat(HX));  dZ := RowWeights(TransposedMat(HZ));
    return rec( wmax := Maximum(Maximum(rwX),Maximum(rwZ),Maximum(dX+dZ)) );
end;
CLPReadN := function( fname, N )
    local s, d, i, M;
    s := StringFile( fname );
    if s = fail then Print("MISSING ",fname,"\n"); return fail; fi;
    d := Filtered( s, c -> c = '0' or c = '1' );
    M := [];
    for i in [ 0 .. Length(d)/N - 1 ] do
        Add( M, List( d{[i*N+1..i*N+N]}, c -> (IntChar(c)-48)*One(GF(2)) ) );
    od;
    return M;
end;
G := [ "G1_n360k16d22", "G2_n360k16d21", "G3_n360k12d22", "G4_n360k12d22",
       "G5_n360k12d22", "G6_n360k16d20", "G7_n360k12d22" ];
PrintTo( "screenG_1e6.txt", "# 10^6 both sides, mindist=2\n" );
for nm in G do
    HX := CLPReadN( Concatenation(nm,"_HX.txt"), 360 );
    HZ := CLPReadN( Concatenation(nm,"_HZ.txt"), 360 );
    if HX <> fail and HZ <> fail then
        K := 360 - RankMat(HX) - RankMat(HZ);
        w := CodeWeights(HX,HZ).wmax;
        dZ := AbsInt( DistRandCSS( HX, HZ, 10^6, 2 : field := GF(2) ) );
        dX := AbsInt( DistRandCSS( HZ, HX, 10^6, 2 : field := GF(2) ) );
        d := Minimum(dZ,dX);
        AppendTo( "screenG_1e6.txt", nm, "  [[360,", K, ",<=", d, "]] w=", w,
                  "  Q=", Float(K*d^2/360), "  dZ=", dZ, " dX=", dX, "\n" );
        Print( nm, ": [[360,", K, ",<=", d, "]] w=", w, " Q=", Float(K*d^2/360), "\n" );
    fi;
od;
QUIT;
