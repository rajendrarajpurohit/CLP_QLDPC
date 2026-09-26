#############################################################################
##  CLP_pairs_final.g  --  enumerate non-normal (G,H) pairs, write pairs_final.g
##  Run once:   gap -q < CLP_pairs_final.g
##  Then in a driver:   Read( "pairs_final.g" );   # defines PAIRS
#############################################################################

ORDERS := Filtered( [ 24 .. 1000 ], o -> NumberSmallGroups( o ) <= 200
                                     and not ( o in [ 256, 512, 768 ] ) );

ATB_IDS := [ [56,9],[72,30],[80,10],[120,32],[126,8],[128,10],[128,26],[135,3],[136,10],[150,8],
             [162,3],[162,4],[168,33],[180,23],[186,3],[216,12],[216,47],[224,53],[248,9],[252,21],
             [288,489],[288,498],[320,10],[320,17],[336,188],[360,99],[384,512] ];

R_ALLOWED  := [ 2, 3, 4, 6 ];
MIN_M      := 10;    MAX_M    := 40;     # index [G:H]  -> block size m
MIN_LIFT   := 4;     MAX_STAB := 4;
PAIRS_FILE := "pairs_final.g";

IDS := Concatenation( ATB_IDS,
         Concatenation( List( ORDERS, o -> List( [ 1 .. NumberSmallGroups( o ) ],
                                                 i -> [ o, i ] ) ) ) );
IDS := Set( IDS );
Print( "scanning ", Length( IDS ), " groups\n" );

PrintTo( PAIRS_FILE, "PAIRS := [\n" );
count := 0;  scanned := 0;

for id in IDS do
    scanned := scanned + 1;
    G := SmallGroup( id[1], id[2] );
    if IsAbelian( G ) then continue; fi;
    isATB := id in ATB_IDS;
    for cc in ConjugacyClassesSubgroups( G ) do
        H := Representative( cc );
        if Size( H ) = 1 or Size( H ) = Size( G ) then continue; fi;
        if IsNormal( G, H ) then continue; fi;          # non-normal only
        m := Index( G, H );
        if m < MIN_M or m > MAX_M then continue; fi;
        N := Normalizer( G, H );
        r := Index( G, N );
        if not ( r in R_ALLOWED ) then continue; fi;
        lift := Index( N, H );
        if lift < MIN_LIFT then continue; fi;
        C := Core( G, H );
        stab := Index( H, C );
        if stab > MAX_STAB then continue; fi;
        if IsPcGroup( G ) then
            gens := List( GeneratorsOfGroup( H ),
                          h -> ExponentsOfPcElement( Pcgs( G ), h ) );
        else
            gens := List( GeneratorsOfGroup( H ), ListPerm );
        fi;
        AppendTo( PAIRS_FILE,
            "rec( l := ", id[1], ", i := ", id[2], ", pc := ", IsPcGroup( G ),
            ", Hgens := ", gens, ", Hdesc := \"", StructureDescription( H ),
            "\", m := ", m, ", r := ", r, ", lift := ", lift, ", stab := ", stab,
            ", core := ", Size( C ), ", atb := ", isATB, " ),\n" );
        count := count + 1;
    od;
    if scanned mod 200 = 0 then
        Print( "  ", scanned, "/", Length( IDS ), " groups, ", count, " pairs\n" );
    fi;
od;

AppendTo( PAIRS_FILE,
    "rec( l := 0, i := 0, pc := false, Hgens := [ ], Hdesc := \"\", m := 0, r := 0, ",
    "lift := 0, stab := 0, core := 0, atb := false ) ];\n" );
Print( "\nwrote ", count, " pairs to ", PAIRS_FILE, "\n" );