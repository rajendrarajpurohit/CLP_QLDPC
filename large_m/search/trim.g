Read( "pairs_merge_fix.g" );;
sel := Filtered( ALLPAIRS, p -> p.lift >= 8 );;
Print( "lift>=8: ", Length( sel ), "\n" );
Reset( GlobalMersenneTwister, 424242 );;
if Length( sel ) > 4000 then sel := Shuffle( sel ){ [ 1 .. 4000 ] }; fi;
Sort( sel, function( a, b ) return a.m < b.m; end );;
PrintTo( "pairs_run.g", "ALLPAIRS := ", sel, ";;\n" );
Print( "run list: ", Length( sel ), " pairs, m from ", sel[1].m, " to ", sel[ Length(sel) ].m, "\n" );
QUIT;
