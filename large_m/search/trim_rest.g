Read( "pairs_merge_fix.g" );;  Read( "pairs_run.g" );;  first := ALLPAIRS;;
Read( "pairs_merge_fix.g" );;
rest := Filtered( ALLPAIRS, p -> p.lift >= 8 and not p in first );;
Sort( rest, function( a, b ) return a.m < b.m; end );;
PrintTo( "pairs_rest.g", "ALLPAIRS := ", rest, ";;\n" );
Print( "rest list: ", Length( rest ), " pairs\n" );
QUIT;
