ALL := [ ];;
for f in Filtered( DirectoryContents( "." ), s -> StartsWith( s, "pairs_big_" ) and EndsWith( s, ".g" ) ) do
  Read( f );  Append( ALL, PAIRS );  Print( f, " ", Length( PAIRS ), "\n" );
od;
PrintTo( "pairs_merge_fix.g", "ALLPAIRS := ", ALL, ";;\n" );
Print( "total pairs ", Length( ALL ), "\n" );
QUIT;
