Read("pairs_bigm.g");;
sel := Filtered( PAIRS, p -> p.l=270 and p.i=26 and p.m=90 );;
Print("(270,26) m=90 entries: ", Length(sel), "\n");
for p in sel do
  Print("  Hdesc=", p.Hdesc, " r=", p.r, " lift=", p.lift,
        " stab=", p.stab, " core=", p.core, " Hgens=", p.Hgens, "\n");
od;
QUIT;
