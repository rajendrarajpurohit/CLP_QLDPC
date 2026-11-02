# usage:  echo 'LLO:=240;; LHI:=500;; MLO:=120;; MHI:=250;; Read("bigpairs.g");' | gap -q
OUT := Concatenation( "pairs_big_", String(LLO), "-", String(LHI), ".g" );;
PrintTo( OUT, "PAIRS := [\n" );
np := 0;;
for l in [ LLO .. LHI ] do
  if NumberSmallGroups( l ) > 2000 then continue; fi;
  if ForAll( [ 2 .. 9 ], h -> not ( l mod h = 0 and l/h >= MLO and l/h <= MHI ) ) then continue; fi;
  for i in [ 1 .. NumberSmallGroups( l ) ] do
    G := SmallGroup( l, i );
    if not IsSolvableGroup( G ) then continue; fi;
    pc := Pcgs( G );
    seen := [ ];
    for cl in ConjugacyClasses( G ) do
      g := Representative( cl );  o := Order( g );
      if o < 2 or o > 9 then continue; fi;
      m := l / o;
      if m < MLO or m > MHI then continue; fi;
      H := Subgroup( G, [ g ] );
      if ForAny( seen, S -> IsConjugate( G, S, H ) ) then continue; fi;
      Add( seen, H );
      if Size( Core( G, H ) ) > 1 then continue; fi;
      N := Normalizer( G, H );  r := Index( G, N );
      if not r in [ 2, 3, 4 ] then continue; fi;
      lift := Index( N, H );
      if lift < 4 then continue; fi;
      np := np + 1;
      AppendTo( OUT, "rec( l := ", l, ", i := ", i, ", pc := true, Hdesc := \"C", o, "\", Hgens := [ ",
                ExponentsOfPcElement( pc, g ), " ], m := ", m, ", r := ", r, ", lift := ", lift,
                ", stab := ", o, ", atb := false ),\n" );
    od;
  od;
  if l mod 50 = 0 then Print( "order ", l, " done, pairs so far ", np, "\n" ); fi;
od;
AppendTo( OUT, "];\n" );
Print( "wrote ", OUT, " with ", np, " pairs\n" );
QUIT;
