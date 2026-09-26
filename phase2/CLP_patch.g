#############################################################################
##  CLP_patch.g  --  final-run overrides.  Read AFTER CLP_core.g
#############################################################################

# ---- hard frontiers: Q = k d^2 / n that a code must BEAT to be logged ----
Q_FRONTIER := [ ];
Q_FRONTIER[5]  :=  6.5;    # essentially empty class; our own [[320,8,16]] is 6.4
Q_FRONTIER[6]  := 15.3;    # [[340,16,18]] Liang, exact, Q=15.25
Q_FRONTIER[7]  := 18.0;    # [[288,16,18]] Qian-Li, exact, Q=18.00
Q_FRONTIER[8]  := 30.6;    # [[378,32,19]] Qian-Li, exact, Q=30.56
Q_FRONTIER[9]  := 20.3;    # [[288,18,18]] Qian-Li, exact, Q=20.25
Q_FRONTIER[10] := 38.8;    # [[234,28,18]] Qian-Li, exact, Q=38.77

RATE_FLOOR  := 0.07;       # drivers override
MIN_D_ABS   := 8;          # never log d < 8 whatever Q says
CRED        := 1.3;        # d is capped at CRED*sqrt(n) when scoring
D_CAP       := 40;

QCapped := function( N, K, d )
    return Float( K * Minimum( Float(d), CRED * Sqrt( Float(N) ) )^2 / N );
end;

Worthy := function( N, K, d, w )
    local ww;
    ww := Maximum( 5, Minimum( w, 10 ) );
    return d >= MIN_D_ABS
       and Float( K / N ) >= RATE_FLOOR
       and QCapped( N, K, d ) > Q_FRONTIER[ ww ];
end;

RequiredD := function( N, K, w )
    local d;
    if Float( K / N ) < RATE_FLOOR then return D_CAP + 1; fi;
    for d in [ MIN_D_ABS .. D_CAP ] do
        if Worthy( N, K, d, w ) then return d; fi;
    od;
    return D_CAP + 1;
end;

# ---- guaranteed dimension from the shape ----
KFloor := function( m, sA, sB )
    return m * ( sA[2] - sA[1] ) * ( sB[1] - sB[2] );
end;

# ---- LP gauge: first row and first column of the base matrix are the identity.
#      free entries = (r-1)(c-1), each a single element.                    ----
GaugeBase := function( nreps, r, c )
    local M, i, j;
    M := [ ];
    for i in [ 1 .. r ] do
        M[i] := [ ];
        for j in [ 1 .. c ] do
            if i = 1 or j = 1 then M[i][j] := [ 1 ];
            else M[i][j] := [ Random( [ 2 .. nreps ] ) ]; fi;
        od;
    od;
    return M;
end;

# same, but one designated free entry carries alpha elements (raises w by alpha-1)
GaugeBaseAlpha := function( nreps, r, c, alpha )
    local M, i, j;
    M := GaugeBase( nreps, r, c );
    if alpha > 1 and r >= 2 and c >= 2 then
        M[ Random([2..r]) ][ Random([2..c]) ] := RandomSubset( [ 2 .. nreps ], alpha );
    fi;
    return M;
end;

# ---- staged distance: 300 -> 5000 -> 60000 -> FULL_ROUNDS ----
DIST_STAGES := [ 300, 5000, 60000 ];

Evaluate := function( P, A, B, sA, sB, wMax, Kmin )
    local code, wts, K, kfl, dReq, d, w, s;
    STATS.cand := STATS.cand + 1;
    if not ForAny( Flat( A ), k -> not P.AinN[k] ) then
        STATS.inN := STATS.inN + 1; return fail; fi;
    code := CLP_Build( P, A, B, sA, sB );
    if HasCancellation( code, A, B, sA, sB, P.m ) then
        STATS.cancel := STATS.cancel + 1; return fail; fi;
    wts := CodeWeights( code.HX, code.HZ );
    if wts.wmax > wMax or wts.wmin < MIN_CHECK_W then
        STATS.weight := STATS.weight + 1; return fail; fi;
    w := wts.wmax;
    K := code.N - RankMat( code.HX ) - RankMat( code.HZ );
    kfl := KFloor( P.m, sA, sB );
    if K < kfl then
        Print( "!! K=", K, " below floor ", kfl, " -- BUILD BUG\n" );
        Error( "k floor violated" );
    fi;
    if K > STATS.bestK then STATS.bestK := K; fi;
    if K < Kmin then STATS.K := STATS.K + 1; return fail; fi;
    dReq := RequiredD( code.N, K, w );
    if dReq > D_CAP then STATS.dreq := STATS.dreq + 1; return fail; fi;
    for s in DIST_STAGES do
        d := DistUB( code.HX, code.HZ, s, dReq - 1 );
        if AbsInt( d ) > STATS.bestd then STATS.bestd := AbsInt( d ); fi;
        if d < dReq then STATS.dscreen := STATS.dscreen + 1; return fail; fi;
    od;
    d := DistUB( code.HX, code.HZ, FULL_ROUNDS, dReq - 1 );
    if d < dReq then STATS.dfull := STATS.dfull + 1; return fail; fi;
    return rec( code := code, K := K, w := w, d := d, dReq := dReq );
end;

# ---- self-describing filenames; no more D3_n collisions ----
LogCode := function( P, T, A, B, res )
    local id, fom, qc, Ael, Bel, key;
    key := [ P.l, P.i, P.m, res.code.N, res.K, res.d, res.w ];
    if key in SEEN_PARAMS then return; fi;
    AddSet( SEEN_PARAMS, key );
    CODE_COUNTER := CODE_COUNTER + 1;  STATS.logged := STATS.logged + 1;
    fom := Float( res.K * res.d^2 / res.code.N );
    qc  := QCapped( res.code.N, res.K, res.d );
    id := Concatenation( OUT_PREFIX, "n", String(res.code.N), "k", String(res.K),
                         "d", String(res.d), "w", String(res.w), "_", String(CODE_COUNTER) );
    Ael := List( A, row -> List( row, e -> List( e, k -> P.Areps[k] ) ) );
    Bel := List( B, row -> List( row, e -> List( e, k -> P.Lreps[k] ) ) );
    Print( "\n>>> GOLD  [[", res.code.N, ",", res.K, ",<=", res.d, "]]  w=", res.w,
           "  Q<=", fom, "  Qcap=", qc, "  rate=", Float(res.K/res.code.N),
           "   ", P.desc, "  ", T.sA, "/", T.sB, "\n" );
    AppendTo( RECORD_FILE,
        "\n#### ", id, "  [[", res.code.N, ",", res.K, ",<=", res.d, "]]  w=", res.w,
        "  Q_upper=", fom, "  Q_capped=", qc, "  rate=", Float(res.K/res.code.N), "\n",
        "driver := \"", T.label, "\";  group := \"", P.desc, "\";  l := ", P.l, ";  i := ", P.i, ";\n",
        "m := ", P.m, ";  r := ", P.r, ";  lift := ", P.lift, ";  stab := ", P.stab,
        ";  kfloor := ", KFloor( P.m, T.sA, T.sB ), ";\n",
        "shapeA := ", T.sA, ";  shapeB := ", T.sB, ";\n",
        "A := ", Ael, ";\nB := ", Bel, ";\n",
        "Aidx := ", A, ";  Bidx := ", B, ";\n",
        "matrices := \"", id, "_HX.txt / _HZ.txt\";\n" );
    WriteBinaryMatrix( Concatenation( id, "_HX.txt" ), res.code.HX );
    WriteBinaryMatrix( Concatenation( id, "_HZ.txt" ), res.code.HZ );
end;

# ---- cap the A-transversal for large G (memory) ----
PrepareActionsCapped := function( P, capA )
    local cos, hom, one, reps, T, Lmat, Rmatperm;
    cos := RightCosets( P.G, P.H );
    hom := ActionHomomorphism( P.G, cos, OnRight );
    one := One( P.G );
    reps := Filtered( AsList( RightTransversal( P.G, P.core ) ), g -> not ( g in P.core ) );
    if Length( reps ) > capA then reps := RandomSubset( reps, capA ); fi;
    P.Areps := Concatenation( [ one ], reps );  P.nA := Length( P.Areps );
    T := Filtered( AsList( RightTransversal( P.N, P.H ) ), u -> not ( u in P.H ) );
    P.Lreps := Concatenation( [ one ], T );  P.nL := Length( P.Lreps );
    P.AinN := List( P.Areps, g -> g in P.N );
    Lmat := g -> PermutationMat( Image( hom, g^-1 )^-1, P.m, GF(2) );
    P.Lmats := List( P.Areps, Lmat );
    Rmatperm := function( u )
        return PermutationMat( PermList( List( cos, c ->
            Position( cos, RightCoset( P.H, u^-1 * Representative( c ) ) ) ) )^-1, P.m, GF(2) );
    end;
    P.Rmats := List( P.Lreps, Rmatperm );
end;

RunDriver2 := function( T )
    local pr, P, seen, cnt, key, A, B, res, pairs, t0;
    if not IsBound( PAIRS ) then Error( "read pairs.g first" ); fi;
    AppendTo( RECORD_FILE, "# ", T.label, "  target w=", T.wMax,
              "  Q must beat ", Q_FRONTIER[T.wMax], "  rate floor ", RATE_FLOOR, "\n" );
    pairs := Filtered( PAIRS, pr -> pr.l > 0 and T.pairOK( pr ) );
    Sort( pairs, function( a, b ) return a.lift > b.lift; end );
    Print( T.label, ":  ", Length( pairs ), " pairs,  kfloor/m = ",
           ( T.sA[2]-T.sA[1] ) * ( T.sB[1]-T.sB[2] ), ",  n/m = ",
           T.sA[2]*T.sB[1] + T.sA[1]*T.sB[2], "\n" );
    for pr in pairs do
        P := LoadPair( pr );
        PrepareActionsCapped( P, 3000 );
        ResetStats();  t0 := Runtime();
        Print( "\n== ", P.desc, "  m=", P.m, " r=", P.r, " lift=", P.lift,
               " stab=", P.stab, " |G/Core|=", P.nA, "  n=", T.Nfrom( P.m ),
               "  kfloor=", KFloor( P.m, T.sA, T.sB ), "\n" );
        seen := [ ];  cnt := 0;
        while cnt < CAND_PER_PAIR do
            cnt := cnt + 1;
            A := T.genA( P );  B := T.genB( P );
            key := [ A, B ];
            if key in seen then continue; fi;
            AddSet( seen, key );
            res := Evaluate( P, A, B, T.sA, T.sB, T.wMax, T.Kmin );
            if res <> fail then LogCode( P, T, A, B, res ); fi;
            if cnt mod 5000 = 0 then
                Print( "   [", cnt, "] ", STATS, "  ", Int((Runtime()-t0)/1000), "s\n" );
            fi;
        od;
        Print( "   pair done: ", STATS, "\n" );
    od;
    Print( "\n", T.label, " finished; ", CODE_COUNTER, " logged\n" );
end;