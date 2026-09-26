if LoadPackage( "QDistRnd" ) = fail then Error("no QDistRnd"); fi;
rw := function( M ) return List( M, r -> Number( r, y -> not IsZero(y) ) ); end;
cw := function( HX, HZ )
    local a, b, c, d;
    a := rw(HX);  b := rw(HZ);
    c := rw(TransposedMat(HX));  d := rw(TransposedMat(HZ));
    return Maximum( Maximum(a), Maximum(b), Maximum(c+d) );
end;
rd := function( f, N )
    local s, x, i, M;
    s := StringFile(f);
    if s = fail then Print("MISSING ",f,"\n"); return fail; fi;
    x := Filtered(s, c -> c='0' or c='1');  M := [];
    for i in [0..Length(x)/N-1] do
        Add(M, List(x{[i*N+1..i*N+N]}, c -> (IntChar(c)-48)*One(GF(2))));
    od;
    return M;
end;
nm := "W1_n360k18d26w8";
HX := rd( Concatenation(nm,"_HX.txt"), 360 );
HZ := rd( Concatenation(nm,"_HZ.txt"), 360 );
Print("file: ", nm, "\n");
Print("CSS ok: ", IsZero(HX*TransposedMat(HZ)), "\n");
K := 360 - RankMat(HX) - RankMat(HZ);
Print("k = ", K, "   w = ", cw(HX,HZ), "\n");
Print("10^6 rounds, Z side...\n");
dZ := AbsInt( DistRandCSS( HX, HZ, 10^6, 2 : field := GF(2) ) );
Print("dZ = ", dZ, "\n");
Print("10^6 rounds, X side...\n");
dX := AbsInt( DistRandCSS( HZ, HX, 10^6, 2 : field := GF(2) ) );
Print("dX = ", dX, "\n");
d := Minimum(dZ,dX);
Print("RESULT [[360,", K, ",<=", d, "]] w=", cw(HX,HZ),
      "  Q=", Float(K*d^2/360), "  d/sqrt(n)=", Float(d/Sqrt(360.0)), "\n");
QUIT;
