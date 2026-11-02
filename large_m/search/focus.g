Read( "pairs_merge_fix.g" );;
want := [ ];;
for line in Filtered( SplitString( StringFile( "focus_groups.txt" ), "\n" ), s -> s <> "" ) do
  s := SplitString( line, " " );
  li := SplitString( s[1]{ [ 12 .. Length(s[1]) - 1 ] }, "," );
  Add( want, [ Int( li[1] ), Int( li[2] ), Int( s[2] ) ] );
od;
sel := Filtered( ALLPAIRS, p -> [ p.l, p.i, p.m ] in want );;
Sort( sel, function( a, b ) return a.m < b.m; end );;
PrintTo( "pairs_focus.g", "ALLPAIRS := ", sel, ";;\n" );
Print( "focus list: ", Length( sel ), " pairs\n" );
QUIT;
