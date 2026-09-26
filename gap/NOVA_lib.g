if LoadPackage( "QDistRnd" ) = fail then Error("no QDistRnd"); fi;

rw := function( M ) return List( M, r -> Number( r, y -> not IsZero(y) ) ); end;

cw := function( HX, HZ )
    local a, b, c, d;
    a := rw(HX);  b := rw(HZ);
    c := rw(TransposedMat(HX));  d := rw(TransposedMat(HZ));
    return Maximum( Maximum(a), Maximum(b), Maximum(c+d) );
end;

dupcol := function( HX, HZ )
    local J;
    J := List( [1..Length(HX[1])], i -> Concatenation(
             List(HX, r -> r[i]), List(HZ, r -> r[i]) ) );
    return Length( Set(J) ) < Length(J);
end;

thincol := function( HX, HZ )
    local cx, cz, i;
    cx := rw( TransposedMat(HX) );
    cz := rw( TransposedMat(HZ) );
    return ForAny( [1..Length(cx)], i -> cx[i] + cz[i] < 3 );
end;

kfloor := function( w )
    if w = 5 then return 6; elif w = 6 then return 10; elif w = 7 then return 10;
    elif w = 8 then return 10; elif w = 9 then return 10; else return 12; fi;
end;

qbar := function( w )
    if w = 5 then return 5.0; elif w = 6 then return 11.0; elif w = 7 then return 13.0;
    elif w = 8 then return 21.0; elif w = 9 then return 14.0; else return 27.0; fi;
end;

ent := function( nreps, a )
    if a = 1 then return [ Random([1..nreps]) ]; fi;
    return Set( Shuffle( [1..nreps] ){[1..a]} );
end;

entid := function( nreps, a )
    return Concatenation( [1], Difference( ent(nreps,a), [1] ) );
end;

arrA := function( nreps, r, c )
    local M, i, j;
    M := [];
    for i in [1..r] do
        M[i] := [];
        for j in [1..c] do M[i][j] := ent( nreps, Random([2,3,4,5]) ); od;
    od;
    M[1][1] := entid( nreps, Random([3,3,3,3,4,4,5]) );
    return M;
end;

arrB := function( nreps, r, c )
    local M, i, j;
    M := [];
    for i in [1..r] do
        M[i] := [];
        for j in [1..c] do M[i][j] := ent( nreps, Random([1,1,1,2,2,3]) ); od;
    od;
    M[1][1] := entid( nreps, Random([1,1,1,2,2,3]) );
    return M;
end;
