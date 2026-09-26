Read("CLP_core.g");;
Read("CLP_patch.g");;
Read("pairs_W6.g");;
hist := [];;
for pr in Filtered( PAIRS, p -> p.l > 0 and p.m in [20,21,24] ) do
    P := LoadPair( pr );;  PrepareActionsCapped( P, 3000 );;
    Print( "== ", P.desc, " m=", P.m, " n=", 18*P.m, " nA=", P.nA, " nL=", P.nL, "\n" );
    for t in [ 1 .. 300 ] do
        A := GaugeBase( P.nA, 3, 3 );;  B := GaugeBase( P.nL, 3, 3 );;
        if ForAny( Flat(A), k -> not P.AinN[k] ) then
            code := CLP_Build( P, A, B, [3,3], [3,3] );;
            if not HasCancellation( code, A, B, [3,3], [3,3], P.m ) then
                K := code.N - RankMat(code.HX) - RankMat(code.HZ);;
                if K >= 16 then
                    d := AbsInt( DistRandCSS( code.HX, code.HZ, 3000, 2 : field := GF(2) ) );;
                    Add( hist, [ code.N, K, d ] );
                    if Float( K*d^2/code.N ) > 12.0 then
                        Print( "  HIT n=", code.N, " K=", K, " d=", d,
                               " Q=", Float(K*d^2/code.N), "\n" );
                    fi;
                fi;
            fi;
        fi;
    od;
    Print( "  running: ", Length(hist), " samples,  K hist ",
           Collected( List(hist, x->x[2]) ), "\n" );
od;
Print( "\nFINAL K hist: ", Collected( List(hist, x->x[2]) ), "\n" );
Print( "FINAL d hist: ", Collected( List(hist, x->x[3]) ), "\n" );
Print( "best Q: ", Maximum( List(hist, x -> Float(x[2]*x[3]^2/x[1]) ) ), "\n" );
QUIT;
