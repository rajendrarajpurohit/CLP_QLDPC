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
    if s = fail then return fail; fi;
    x := Filtered(s, c -> c='0' or c='1');  M := [];
    for i in [0..Length(x)/N-1] do
        Add(M, List(x{[i*N+1..i*N+N]}, c -> (IntChar(c)-48)*One(GF(2))));
    od;
    return M;
end;
TOP := [ "H75_n360k20d22", "H162_n360k16d24", "H243_n360k16d24",
         "H34_n360k16d23", "H164_n360k16d23", "H83_n360k20d20",
         "H221_n360k20d20", "H152_n360k20d19" ];
PrintTo( "screenTop_1e6.txt", "# 10^6 rounds each side, mindist=2, no early exit\n" );
for nm in TOP do
    HX := rd( Concatenation(nm,"_HX.txt"), 360 );
    HZ := rd( Concatenation(nm,"_HZ.txt"), 360 );
    if HX = fail or HZ = fail then
        Print( nm, ": MISSING\n" );
        AppendTo( "screenTop_1e6.txt", nm, "  MISSING\n" );
    elif not IsZero( HX * TransposedMat(HZ) ) then
        Print( nm, ": CSS FAILS\n" );
        AppendTo( "screenTop_1e6.txt", nm, "  CSS FAILS\n" );
    else
        K := 360 - RankMat(HX) - RankMat(HZ);
        w := cw(HX,HZ);
        dZ := AbsInt( DistRandCSS( HX, HZ, 10^6, 2 : field := GF(2) ) );
        dX := AbsInt( DistRandCSS( HZ, HX, 10^6, 2 : field := GF(2) ) );
        d := Minimum(dZ,dX);
        AppendTo( "screenTop_1e6.txt", nm, "  [[360,", K, ",<=", d, "]]  w=", w,
            "  Q=", Float(K*d^2/360), "  dZ=", dZ, " dX=", dX,
            "  d/sqrt(n)=", Float(d/Sqrt(360.0)), "\n" );
        Print( nm, ": [[360,", K, ",<=", d, "]] w=", w,
               "  Q=", Float(K*d^2/360), "  (dZ=", dZ, " dX=", dX, ")\n" );
    fi;
od;
Print( "screen complete\n" );
QUIT;
