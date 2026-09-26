if LoadPackage( "QDistRnd" ) = fail then Error("no QDistRnd"); fi;
rw := function( M ) return List( M, r -> Number( r, y -> not IsZero(y) ) ); end;
cw := function( HX, HZ )
    local a, b, c, d;
    a := rw(HX);  b := rw(HZ);
    c := rw(TransposedMat(HX));  d := rw(TransposedMat(HZ));
    return Maximum( Maximum(a), Maximum(b), Maximum(c+d) );
end;
qbar := function( w )
    if w = 5 then return 6.0; elif w = 6 then return 12.0; elif w = 7 then return 14.0;
    elif w = 8 then return 20.0; elif w = 9 then return 16.0; else return 25.0; fi;
end;
rdN := function( f, N )
    local s, x, i, M;
    s := StringFile(f);
    if s = fail then return fail; fi;
    x := Filtered(s, c -> c='0' or c='1');
    if Length(x) mod N <> 0 then return fail; fi;
    M := [];
    for i in [0..Length(x)/N-1] do
        Add(M, List(x{[i*N+1..i*N+N]}, c -> (IntChar(c)-48)*One(GF(2))));
    od;
    return M;
end;

files := Filtered( SplitString( StringFile("stems.txt"), "\n" ), f -> Length(f) > 0 );;
Print( "stems: ", Length(files), "\n" );
Exec( "mkdir -p pass4" );
PrintTo( "screen4_results.txt", "# stem n k w d Q ratio\n" );
npass := 0;;
for stem in files do
    parts := SplitString( stem, "_" );;
    nstr := parts[Length(parts)];;
    npos := Position( nstr, 'n' );;
    kpos := Position( nstr, 'k' );;
    N := Int( nstr{[npos+1..kpos-1]} );;
    HX := rdN( Concatenation(stem,"_HX.txt"), N );;
    HZ := rdN( Concatenation(stem,"_HZ.txt"), N );;
    if HX = fail or HZ = fail then
        Print( stem, ": unreadable at N=", N, "\n" );
    elif not IsZero( HX*TransposedMat(HZ) ) then
        Print( stem, ": CSS FAILS\n" );
    else
        K := N - RankMat(HX) - RankMat(HZ);;
        w := cw(HX,HZ);;
        dZ := AbsInt( DistRandCSS( HX, HZ, 10^4, 2 : field := GF(2) ) );;
        dX := AbsInt( DistRandCSS( HZ, HX, 10^4, 2 : field := GF(2) ) );;
        d := Minimum(dZ,dX);;
        q := Float( K*d^2/N );;
        AppendTo( "screen4_results.txt", stem, " ", N, " ", K, " ", w, " ", d,
                  " ", q, " ", Float(d/Sqrt(Float(N))), "\n" );
        if q > qbar(w) then
            npass := npass + 1;
            Exec( Concatenation("cp ", stem, "_HX.txt ", stem, "_HZ.txt pass4/") );
            Print( "PASS ", stem, " [[", N, ",", K, ",<=", d, "]] w=", w, " Q=", q, "\n" );
        else
            Print( "drop ", stem, " Q=", q, "\n" );
        fi;
    fi;
od;
Print( "\n", npass, " of ", Length(files), " passed 10^4\n" );
QUIT;
