#############################################################################
##  CLP_pairs.g  --  enumerate (G,H) pairs ATB-style, write pairs.g
#############################################################################

ORDERS    := Filtered( [ 48 .. 250 ], o -> NumberSmallGroups( o ) <= 80 );
ATB_IDS   := [ [56,9],[72,30],[80,10],[120,32],[126,8],[128,10],[128,26],[135,3],[136,10],[150,8],
               [162,3],[162,4],[168,33],[180,23],[186,3],[216,12],[216,47],[224,53],[248,9],[252,21],
               [288,489],[288,498],[320,10],[320,17],[336,188],[360,99],[384,512] ];
R_ALLOWED := [ 2, 3, 4 ];
MIN_M := 12;  MAX_M := 130;  MIN_LIFT := 6;  MAX_STAB := 3;
PAIRS_FILE := "pairs.g";

IDS := Concatenation( ATB_IDS, Concatenation( List( ORDERS, o -> List( [ 1 .. NumberSmallGroups( o ) ], i -> [ o, i ] ) ) ) );
IDS := Set( IDS );

PrintTo( PAIRS_FILE, "PAIRS := [\n" );
count := 0;
for id in IDS do
    G := SmallGroup( id[1], id[2] );
    if IsAbelian( G ) then continue; fi;
    isATB := id in ATB_IDS;
    for cc in ConjugacyClassesSubgroups( G ) do
        H := Representative( cc );
        if Size( H ) = 1 or Size( H ) = Size( G ) then continue; fi;
        if IsNormal( G, H ) then continue; fi;
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
            gens := List( GeneratorsOfGroup( H ), h -> ExponentsOfPcElement( Pcgs( G ), h ) );
        else
            gens := List( GeneratorsOfGroup( H ), ListPerm );
        fi;
        AppendTo( PAIRS_FILE, "rec( l := ", id[1], ", i := ", id[2], ", pc := ", IsPcGroup( G ),
                  ", Hgens := ", gens, ", Hdesc := \"", StructureDescription( H ), "\", m := ", m,
                  ", r := ", r, ", lift := ", lift, ", stab := ", stab, ", core := ", Size( C ),
                  ", atb := ", isATB, " ),\n" );
        count := count + 1;
    od;
    Print( "SmallGroup", id, " done, pairs so far: ", count, "\n" );
od;
AppendTo( PAIRS_FILE, "rec( l := 0, i := 0, pc := false, Hgens := [ ], Hdesc := \"\", m := 0, r := 0, lift := 0, stab := 0, core := 0, atb := false ) ];\n" );
Print( "wrote ", count, " pairs to ", PAIRS_FILE, "\n" );