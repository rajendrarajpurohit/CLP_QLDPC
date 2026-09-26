Read("CLP_core.g");;
Read("CLP_patch.g");;
Read("pairs_final.g");;
hist := [];;
for pr in Filtered( PAIRS, p -> p.l > 0 and p.m = 21 ) do
    P := LoadPair( pr );;  PrepareActionsCapped( P, 3000 );;
    Print( "== ", P.desc, " m=", P.m, " r=", P.r, " n=", 18*P.m, "\n" );
    for t in [ 1 .. 400 ] do
        A := GaugeBase( P.nA, 3, 3 );;  B := GaugeBase( P.nL, 3, 3 );;
        if ForAny( Flat(A), k -> not P.AinN[k] ) then
            code := CLP_Build( P, A, B, [3,3], [3,3] );;
            if not HasCancellation( code, A, B, [3,3], [3,3], P.m ) then
                K := code.N - RankMat(code.HX) - RankMat(code.HZ);;
                if K >= 28 then
                    d := AbsInt( DistRandCSS( code.HX, code.HZ, 3000, 2 : field := GF(2) ) );;
                    Add( hist, [ code.N, K, d, P.desc ] );
                    if Float( K*d^2/code.N ) > 13.0 then
                        Print( "  HIT ", P.desc, " [[", code.N, ",", K, ",<=", d,
                               "]] Q=", Float(K*d^2/code.N), "\n" );
                    fi;
                fi;
            fi;
        fi;
    od;
    if Length(hist) > 0 then
        Print( "  ", Length(hist), " samples, best Q=",
               Maximum( List(hist, x -> Float(x[2]*x[3]^2/x[1])) ), "\n" );
    fi;
od;
QUIT;
