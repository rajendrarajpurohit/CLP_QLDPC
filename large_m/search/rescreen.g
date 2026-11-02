LoadPackage("QDistRnd");;
ReadBin := function( f )
local s, lines;
s := StringFile( f );
lines := Filtered( SplitString( s, "\n" ), l -> Length( l ) > 0 );
return List( lines, l -> List( l, c -> [ 0*Z(2), Z(2)^0 ][ IntChar(c) - 47 ] ) );
end;;
Rescreen := function( id, rounds )
local HX, HZ, n, k, dZ, dX;
HX := ReadBin( Concatenation( id, "_HX.txt" ) );
HZ := ReadBin( Concatenation( id, "_HZ.txt" ) );
if not IsZero( HX * TransposedMat( HZ ) ) then Error( id, ": not CSS" ); fi;
n := Length( HX[1] );
k := n - RankMat( HX ) - RankMat( HZ );
dZ := DistRandCSS( HX, HZ, rounds, 0 : field := GF(2) );
dX := DistRandCSS( HZ, HX, rounds, 0 : field := GF(2) );
AppendTo( "rescreen_results.txt", id, "  [[", n, ",", k, ",<=", Minimum( dZ, dX ), "]]  dZ<=", dZ, "  dX<=", dX, "  rounds=", rounds, "\n" );
end;;
