Read("gap/CLP_core.g");;

readmat := function( f, N )
    local s, d, i, M;
    s := StringFile(f);  if s = fail then return fail; fi;
    d := Filtered(s, c -> c='0' or c='1');  M := [];
    for i in [0..Length(d)/N-1] do
        Add(M, List(d{[i*N+1..i*N+N]}, c -> (IntChar(c)-48)*One(GF(2))));
    od;
    return M;
end;;

PosMod := function( x, reps, S )
    local i;
    i := Position( reps, x );
    if i <> fail then return i; fi;
    for i in [1..Length(reps)] do
        if x * reps[i]^-1 in S then return i; fi;
    od;
    return fail;
end;;
Recover := function( l, i, mWant, liftWant, stabWant, hdesc, aWords, bWords, hxfile, hzfile, n, outfile )
    local G, f, aEl, bEl, cands, cc, H, m, N, r, lift, C, stab, gens, rec1, P, allB, someA, Ai, Bi, c, HX, HZ, same, k;
    G := SmallGroup( l, i );  f := GeneratorsOfGroup( G );
    aEl := aWords( f );  bEl := bWords( f );
    HX := readmat( hxfile, n );  HZ := readmat( hzfile, n );
    PrintTo( outfile, "== SmallGroup(", l, ",", i, ")  n=", n, "\n" );
    AppendTo( outfile, "G = ", StructureDescription(G), "\n" );
    AppendTo( outfile, "a words  ", aEl, "\n" );
    AppendTo( outfile, "B words  ", bEl, "\n\n" );
    cands := [];
    for cc in ConjugacyClassesSubgroups( G ) do
        H := Representative( cc );
        if Size(H) = 1 or Size(H) = Size(G) or IsNormal(G,H) then continue; fi;
        m := Index( G, H );  if m <> mWant then continue; fi;
        N := Normalizer( G, H );  r := Index( G, N );  if not r in [2,3,4] then continue; fi;
        lift := Index( N, H );  C := Core( G, H );  stab := Index( H, C );
        if StructureDescription(H) <> hdesc then continue; fi;
        gens := List( GeneratorsOfGroup(H), h -> ExponentsOfPcElement( Pcgs(G), h ) );
        rec1 := rec( l := l, i := i, pc := true, Hgens := gens, Hdesc := hdesc, m := m, r := r, lift := lift, stab := stab, core := Size(C), atb := true );
        allB := ForAll( Flat(bEl), x -> x in N );
        someA := ForAny( aEl, x -> not x in N );
        AppendTo( outfile, "candidate ", gens, " lift=", lift, " stab=", stab, " r=", r, " allB=", allB, " someA=", someA, "\n" );
        if allB and someA and lift = liftWant and stab = stabWant then Add( cands, rec1 ); fi;
    od;
    AppendTo( outfile, "\nsurviving: ", Length(cands), "\n\n" );
    for rec1 in cands do
        P := LoadPair( rec1 );  PrepareActions( P );  f := GeneratorsOfGroup( P.G );  aEl := aWords( f );  bEl := bWords( f );
        Ai := List( aEl, x -> PosMod( x, P.Areps, P.core ) );
        Bi := List( bEl, row -> List( row, e -> List( e, x -> PosMod( x, P.Lreps, P.H ) ) ) );
        AppendTo( outfile, "H gens ", rec1.Hgens, "\n  Aidx ", Ai, "  Bidx ", Bi, "\n" );
        if fail in Ai or fail in Flat(Bi) then AppendTo( outfile, "  index matching FAILED\n\n" ); continue; fi;
        c := CLP_Build( P, [ [ Ai ] ], Bi, [1,1], [2,2] );
        same := ( c.HX = HX ) and ( c.HZ = HZ );
        k := n - RankMat(c.HX) - RankMat(c.HZ);
        AppendTo( outfile, "  rebuild identical: ", same, "  k=", k, "\n" );
        if same then
            AppendTo( outfile, "\nPINNED H = ", List( GeneratorsOfGroup(P.H), String ), " ", StructureDescription(P.H), "\n" );
            AppendTo( outfile, "  N_G(H) = ", StructureDescription(P.N), "  Core = ", StructureDescription(P.core), "\n" );
            AppendTo( outfile, "  m=", P.m, " r=", P.r, " lift=", P.lift, " stab=", P.stab, "\n\n" );
        fi;
    od;
end;;
Recover( 360, 99, 60, 30, 2, "C6",
  f -> [ One(f[1]), f[4]^4*f[5]*f[6]^2, f[1]*f[2]*f[4]^4*f[6]^2 ],
  f -> [ [ [ One(f[1]) ], [ One(f[1]), f[4]^2*f[5]*f[6]^2 ] ], [ [ One(f[1]), f[4]^3*f[5]*f[6]^2 ], [ f[4]^4*f[6]^2 ] ] ],
  "codes/240_8_20/D3_4_HX.txt", "codes/240_8_20/D3_4_HZ.txt", 240, "codes/240_8_20/D3_4_elements.txt" );

Recover( 252, 21, 42, 21, 3, "C6",
  f -> [ One(f[1]), f[3]^2*f[5]^2, f[1]*f[3]^6*f[5] ],
  f -> [ [ [ One(f[1]), f[3]^2*f[5]^2 ], [ One(f[1]) ] ], [ [ One(f[1]) ], [ f[5], f[3]^3 ] ] ],
  "codes/168_10_14/D3_6_HX.txt", "codes/168_10_14/D3_6_HZ.txt", 168, "codes/168_10_14/D3_6_elements.txt" );


RecoverGen := function( l, i, mWant, liftWant, stabWant, hdesc, aWords, bWords, sA, sB, hxfile, hzfile, n, outfile )
    local G, f, aEl, bEl, cands, cc, H, m, N, r, lift, C, stab, gens, rec1, P, allB, someA, Ai, Bi, c, HX, HZ, same, k;
    G := SmallGroup( l, i );  f := GeneratorsOfGroup( G );
    aEl := aWords( f );  bEl := bWords( f );
    HX := readmat( hxfile, n );  HZ := readmat( hzfile, n );
    PrintTo( outfile, "== SmallGroup(", l, ",", i, ")  n=", n, "  layout ", sA, "/", sB, "\n" );
    AppendTo( outfile, "G = ", StructureDescription(G), "\na words ", aEl, "\nB words ", bEl, "\n\n" );
    cands := [];
    for cc in ConjugacyClassesSubgroups( G ) do
        H := Representative( cc );
        if Size(H) = 1 or Size(H) = Size(G) or IsNormal(G,H) then continue; fi;
        m := Index( G, H );  if m <> mWant then continue; fi;
        N := Normalizer( G, H );  r := Index( G, N );  if not r in [2,3,4] then continue; fi;
        lift := Index( N, H );  C := Core( G, H );  stab := Index( H, C );
        if StructureDescription(H) <> hdesc then continue; fi;
        gens := List( GeneratorsOfGroup(H), h -> ExponentsOfPcElement( Pcgs(G), h ) );
        rec1 := rec( l := l, i := i, pc := true, Hgens := gens, Hdesc := hdesc, m := m, r := r, lift := lift, stab := stab, core := Size(C), atb := true );
        allB := ForAll( Flat(bEl), x -> x in N );
        someA := ForAny( Flat(aEl), x -> not x in N );
        AppendTo( outfile, "candidate ", gens, " lift=", lift, " stab=", stab, " r=", r, " allB=", allB, " someA=", someA, "\n" );
        if allB and someA and lift = liftWant and stab = stabWant then Add( cands, rec1 ); fi;
    od;
    AppendTo( outfile, "\nsurviving: ", Length(cands), "\n\n" );
    for rec1 in cands do
        P := LoadPair( rec1 );  PrepareActions( P );  f := GeneratorsOfGroup( P.G );  aEl := aWords( f );  bEl := bWords( f );
        Ai := List( aEl, row -> List( row, e -> List( e, x -> PosMod( x, P.Areps, P.core ) ) ) );
        Bi := List( bEl, row -> List( row, e -> List( e, x -> PosMod( x, P.Lreps, P.H ) ) ) );
        AppendTo( outfile, "H gens ", rec1.Hgens, "\n  Aidx ", Ai, "  Bidx ", Bi, "\n" );
        if fail in Flat(Ai) or fail in Flat(Bi) then AppendTo( outfile, "  index matching FAILED\n\n" ); continue; fi;
        c := CLP_Build( P, Ai, Bi, sA, sB );
        same := ( c.HX = HX ) and ( c.HZ = HZ );
        k := n - RankMat(c.HX) - RankMat(c.HZ);
        AppendTo( outfile, "  rebuild identical: ", same, "  k=", k, "\n" );
        if same then
            AppendTo( outfile, "\nPINNED H = ", List( GeneratorsOfGroup(P.H), String ), " ", StructureDescription(P.H), "\n" );
            AppendTo( outfile, "  N_G(H) = ", StructureDescription(P.N), "  Core = ", StructureDescription(P.core), "  m=", P.m, " r=", P.r, " lift=", P.lift, " stab=", P.stab, "\n\n" );
        fi;
    od;
end;;
EXD := "codes/phase1_extras/";;
RecoverGen( 360, 99, 30, 15, 3, "C6 x C2",
  f -> [ [ [ One(f[1]) ], [ f[1]*f[3]^2*f[4]^2*f[6]^2 ] ], [ [ One(f[1]) ], [ f[1]*f[3]^2*f[6]^2 ] ] ],
  f -> [ [ [ One(f[1]), f[6], f[4]^4*f[6]^2 ] ] ], [2,2], [1,1],
  Concatenation(EXD,"D2_15/D2_15_HX.txt"), Concatenation(EXD,"D2_15/D2_15_HZ.txt"), 120, Concatenation(EXD,"D2_15/D2_15_elements.txt") );
RecoverGen( 135, 3, 45, 15, 3, "C3",
  f -> [ [ [ One(f[1]) ], [ f[1]^2*f[2]^2*f[3]^3 ] ], [ [ One(f[1]) ], [ f[1]^2*f[2]^2*f[3]^4 ] ] ],
  f -> [ [ [ One(f[1]), f[4], f[3]^3*f[4]^2 ] ] ], [2,2], [1,1],
  Concatenation(EXD,"D2_23/D2_23_HX.txt"), Concatenation(EXD,"D2_23/D2_23_HZ.txt"), 180, Concatenation(EXD,"D2_23/D2_23_elements.txt") );
RecoverGen( 360, 99, 60, 30, 2, "C6",
  f -> [ [ [ One(f[1]) ], [ f[1]*f[2]*f[4]^3*f[5]*f[6]^2 ] ], [ [ One(f[1]) ], [ f[1]*f[2]*f[4]*f[5]*f[6]^2 ] ] ],
  f -> [ [ [ One(f[1]), f[4]*f[6], f[4]*f[6]^2 ] ] ], [2,2], [1,1],
  Concatenation(EXD,"D2_4/D2_4_HX.txt"), Concatenation(EXD,"D2_4/D2_4_HZ.txt"), 240, Concatenation(EXD,"D2_4/D2_4_elements.txt") );
RecoverGen( 336, 188, 28, 14, 2, "D12",
  f -> [ [ [ One(f[1]), f[3]*f[4]^4, f[3]*f[4]^6*f[5] ] ] ],
  f -> [ [ [ One(f[1]), f[4]^2 ], [ One(f[1]), f[4]^4*f[5] ] ], [ [ One(f[1]), f[4]^3*f[5] ], [ f[4]^3*f[5], f[4]^6*f[5] ] ] ], [1,1], [2,2],
  Concatenation(EXD,"D3_41/D3_41_HX.txt"), Concatenation(EXD,"D3_41/D3_41_HZ.txt"), 112, Concatenation(EXD,"D3_41/D3_41_elements.txt") );
RecoverGen( 336, 188, 42, 14, 2, "D8",
  f -> [ [ [ One(f[1]), f[3]*f[4]^4*f[6]^2, f[3]*f[4]^6*f[6] ] ] ],
  f -> [ [ [ One(f[1]), f[3]*f[4]^6 ], [ One(f[1]), f[4]^6 ] ], [ [ One(f[1]), f[4]^6 ], [ f[4], f[3]*f[4]^3 ] ] ], [1,1], [2,2],
  Concatenation(EXD,"D3_39/D3_39_HX.txt"), Concatenation(EXD,"D3_39/D3_39_HZ.txt"), 168, Concatenation(EXD,"D3_39/D3_39_elements.txt") );
RecoverGen( 336, 188, 42, 14, 2, "D8",
  f -> [ [ [ One(f[1]), f[4]*f[6]^2, f[4]^5*f[6] ] ] ],
  f -> [ [ [ One(f[1]), f[2]*f[4] ], [ One(f[1]), f[2]*f[4]^3 ] ], [ [ One(f[1]), f[2]*f[4]^5 ], [ f[4]^2, f[4]^4 ] ] ], [1,1], [2,2],
  Concatenation(EXD,"D3_56/D3_56_HX.txt"), Concatenation(EXD,"D3_56/D3_56_HZ.txt"), 168, Concatenation(EXD,"D3_56/D3_56_elements.txt") );
RecoverGen( 135, 3, 45, 15, 3, "C3",
  f -> [ [ [ One(f[1]), f[1]*f[2]^2*f[3]^4*f[4]^2, f[1]^2*f[2]*f[3]*f[4] ] ] ],
  f -> [ [ [ One(f[1]), f[3]^2 ], [ One(f[1]), f[3]*f[4]^2 ] ], [ [ One(f[1]), f[3]^3*f[4]^2 ], [ f[3]*f[4]^2, f[3]^2*f[4]^2 ] ] ], [1,1], [2,2],
  Concatenation(EXD,"D3_17/D3_17_HX.txt"), Concatenation(EXD,"D3_17/D3_17_HZ.txt"), 180, Concatenation(EXD,"D3_17/D3_17_elements.txt") );
RecoverGen( 336, 188, 56, 28, 2, "S3",
  f -> [ [ [ One(f[1]), f[3]*f[4], f[3]*f[4]^5*f[5] ] ] ],
  f -> [ [ [ One(f[1]), f[4]*f[5] ], [ One(f[1]), f[4]^6 ] ], [ [ One(f[1]) ], [ f[4], f[4]^2*f[5] ] ] ], [1,1], [2,2],
  Concatenation(EXD,"D3_5/D3_5_HX.txt"), Concatenation(EXD,"D3_5/D3_5_HZ.txt"), 224, Concatenation(EXD,"D3_5/D3_5_elements.txt") );
RecoverGen( 360, 99, 60, 30, 3, "C6",
  f -> [ [ [ One(f[1]), f[1]*f[4]^3*f[6]^2, f[1]*f[2]*f[3]^2*f[4]^4*f[5]*f[6]^2 ] ] ],
  f -> [ [ [ One(f[1]) ], [ One(f[1]), f[2]*f[4]^3*f[6]^2 ] ], [ [ One(f[1]), f[2]*f[6] ], [ f[2]*f[4]^2*f[6]^2, f[4]^4*f[6] ] ] ], [1,1], [2,2],
  Concatenation(EXD,"D3_3/D3_3_HX.txt"), Concatenation(EXD,"D3_3/D3_3_HZ.txt"), 240, Concatenation(EXD,"D3_3/D3_3_elements.txt") );
RecoverGen( 126, 8, 63, 21, 2, "C2",
  f -> [ [ [ One(f[1]), f[3]*f[4]^6, f[3]^2*f[4]^4 ] ] ],
  f -> [ [ [ One(f[1]), f[2]^2*f[4]^4 ], [ One(f[1]), f[2] ] ], [ [ One(f[1]), f[2]*f[4]^3 ], [ f[2]^2, f[2]*f[4]^6 ] ] ], [1,1], [2,2],
  Concatenation(EXD,"D3_8/D3_8_HX.txt"), Concatenation(EXD,"D3_8/D3_8_HZ.txt"), 252, Concatenation(EXD,"D3_8/D3_8_elements.txt") );
Print("extras done\n");
QUIT;
