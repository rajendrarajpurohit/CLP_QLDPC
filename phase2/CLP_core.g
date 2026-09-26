#############################################################################
##  CLP_core.g  --  shared machinery for the ATB-style CLP search
##  Drivers: Read("CLP_core.g"); Read("pairs.g"); set options; RunDriver(T);
#############################################################################

if LoadPackage( "QDistRnd" ) = fail then Error( "QDistRnd not found" ); fi;

#############################################################################
##  defaults (drivers override after reading this file)
#############################################################################

SCREEN_ROUNDS := 2000;
FULL_ROUNDS   := 20000;
MIN_CHECK_W   := 4;
CAND_PER_PAIR := 20000;
D_CAP         := 40;
OUT_PREFIX    := "CLP_";
RECORD_FILE   := "CLP_results.txt";
MIN_D_LOG := 8;       # never log below this distance
FOM_FLOOR := 4.0;     # never log below this K d^2 / N
SEEN_PARAMS := [ ];   # de-duplicate equivalent hits per pair
VALIDATE_TARGETS := [ ];          # e.g. [ [56,6,8], [36,4,6] ] : accept these (N,K,d) regardless of tables

FOM_COMP := [ ];
FOM_COMP[6] := 12.0;  FOM_COMP[7] := 14.0;  FOM_COMP[8] := 18.0;  FOM_COMP[9] := 20.0;  FOM_COMP[10] := 25.0;

# [N, K, d, w]: ATB Tables II+VI, multi-agent Tables I+II, BB / twisted-torus
KNOWN := [
 [36,4,6,6],[48,4,8,6],[48,8,6,6],[54,8,6,6],[56,6,8,6],[60,16,4,6],[72,4,10,6],[72,8,8,6],[84,6,10,6],
 [90,8,10,6],[96,4,12,6],[96,8,10,6],[108,8,10,6],[112,6,12,6],[112,12,8,6],[120,8,12,6],[126,10,10,6],
 [144,4,16,6],[144,12,12,6],[150,16,8,6],[162,8,14,6],[168,16,10,6],[180,8,16,6],[186,10,14,6],[192,8,16,6],
 [210,10,16,6],[224,12,16,6],[234,8,18,6],[248,10,18,6],[254,14,16,6],[270,8,20,6],
 [48,10,6,8],[72,16,6,8],[80,8,10,8],[80,10,8,8],[84,10,9,8],[84,16,8,8],[90,18,7,8],[96,10,12,8],[96,12,10,8],
 [96,18,8,8],[108,12,10,8],[112,12,12,8],[112,16,10,8],[120,14,12,8],[120,16,11,8],[120,24,7,8],[128,10,14,8],
 [128,16,12,8],[128,18,10,8],[128,30,8,8],[136,8,15,8],[144,10,15,8],[144,18,12,8],[144,24,8,8],[160,8,18,8],
 [160,12,16,8],[168,10,17,8],[168,16,15,8],[224,22,16,8],[256,18,16,8],[288,24,18,8],
 [288,16,18,7],[288,18,18,9],[320,24,16,9],[170,32,14,10],[234,28,18,10],[248,12,18,10] ];

STATS := rec( cand := 0, inN := 0, cancel := 0, weight := 0, K := 0, dreq := 0, dscreen := 0, dfull := 0,
              logged := 0, bestd := 0, bestK := 0 );

ResetStats := function()
    local f;
    for f in RecNames( STATS ) do STATS.( f ) := 0; od;
end;

#############################################################################
##  criteria
#############################################################################

FOMThreshold := function( w )
    local ww;
    ww := Maximum( 6, Minimum( w, 10 ) );
    if IsBound( FOM_COMP[ww] ) then return FOM_COMP[ww]; fi;
    return Float( 10^9 );
end;

# dominated by a known code or by any t-fold direct sum of one (same weight class or lighter)
IsDominated := function( N, K, d, w )
    return ForAny( KNOWN, c -> ForAny( [ 1 .. Int( N / c[1] ) ],
        t -> t * c[2] >= K and c[3] >= d and c[4] <= w ) );
end;

# smallest d that would make [[N,K,d]] at weight w worth logging
RequiredD := function( N, K, w )
    local d, tgt;
    tgt := First( VALIDATE_TARGETS, t -> t[1] = N and t[2] = K );
    if tgt <> fail then return tgt[3]; fi;
    for d in [ MIN_D_LOG .. D_CAP ] do
        if Float( K * d^2 / N ) >= FOM_FLOOR
           and ( ( not IsDominated( N, K, d, w ) ) or Float( K * d^2 / N ) >= FOMThreshold( w ) ) then
            return d;
        fi;
    od;
    return D_CAP + 1;
end;

#############################################################################
##  utilities
#############################################################################

RowWeights := function( M ) return List( M, r -> Number( r, y -> not IsZero( y ) ) ); end;

WriteBinaryMatrix := function( fname, M )
    local f, r;
    f := OutputTextFile( fname, false );
    SetPrintFormattingStatus( f, false );
    for r in M do
        WriteLine( f, JoinStringsWithSeparator( List( r, y -> String( IntFFE( y ) ) ), "" ) );
    od;
    CloseStream( f );
end;

DistUB := function( HX, HZ, num, mindist )
    local dZ, dX;
    dZ := DistRandCSS( HX, HZ, num, mindist : field := GF(2) );
    if dZ <= mindist then return dZ; fi;
    dX := DistRandCSS( HZ, HX, num, mindist : field := GF(2) );
    if dX <= mindist then return dX; fi;
    return Minimum( dZ, dX );
end;

RandomSubset := function( pool, k )
    local s, x;
    s := [ ];
    while Length( s ) < k do
        x := Random( pool );
        if not ( x in s ) then Add( s, x ); fi;
    od;
    return Set( s );
end;

# support of t indices with the identity (index 1) fixed in it  (Lemma IV.1 normalisation)
SupportWithOne := function( nreps, t )
    if t <= 1 then return [ 1 ]; fi;
    return Concatenation( [ 1 ], RandomSubset( [ 2 .. nreps ], t - 1 ) );
end;

RandomEntry := function( nreps, tmin, tmax )
    return RandomSubset( [ 1 .. nreps ], Random( [ tmin .. tmax ] ) );
end;

#############################################################################
##  pairs
#############################################################################

LoadPair := function( pr )
    local G, H, N, C, P;
    G := SmallGroup( pr.l, pr.i );
    if pr.pc then
        H := Subgroup( G, List( pr.Hgens, v -> PcElementByExponents( Pcgs( G ), v ) ) );
    else
        H := Subgroup( G, List( pr.Hgens, PermList ) );
    fi;
    N := Normalizer( G, H );  C := Core( G, H );
    P := rec( G := G, H := H, N := N, core := C, m := Index( G, H ), r := Index( G, N ),
              lift := Index( N, H ), stab := Index( H, C ), l := pr.l, i := pr.i,
              desc := Concatenation( "SmallGroup(", String( pr.l ), ",", String( pr.i ), ") H=", pr.Hdesc ) );
    return P;
end;

# left action of G on G/H (distinct matrices indexed by G/Core reps, identity first)
# right action of N_G(H) on G/H (distinct matrices indexed by N/H reps, identity first)
PrepareActions := function( P )
    local cos, hom, reps, one, T, t, perm, Lmat, Rmatperm;
    cos := RightCosets( P.G, P.H );
    hom := ActionHomomorphism( P.G, cos, OnRight );
    one := One( P.G );
    reps := Filtered( AsList( RightTransversal( P.G, P.core ) ), g -> not ( g in P.core ) );
    P.Areps := Concatenation( [ one ], reps );
    P.nA := Length( P.Areps );
    T := Filtered( AsList( RightTransversal( P.N, P.H ) ), u -> not ( u in P.H ) );
    P.Lreps := Concatenation( [ one ], T );
    P.nL := Length( P.Lreps );
    P.AinN := List( P.Areps, g -> g in P.N );
    Lmat := function( g ) return PermutationMat( Image( hom, g^-1 )^-1, P.m, GF(2) ); end;
    P.Lmats := List( P.Areps, Lmat );
    Rmatperm := function( u )
        local perm;
        perm := PermList( List( cos, c -> Position( cos, RightCoset( P.H, u^-1 * Representative( c ) ) ) ) );
        return PermutationMat( perm^-1, P.m, GF(2) );
    end;
    P.Rmats := List( P.Lreps, Rmatperm );
end;

#############################################################################
##  build / evaluate
#############################################################################

CLP_Build := function( P, A, B, sA, sB )
    local nA0, nA1, nB0, nB1, m, Nq, off, LA, RB, HX, HZ, i, j, k, l, rb, AddBlock, Blk;
    nA0 := sA[1]; nA1 := sA[2]; nB0 := sB[1]; nB1 := sB[2]; m := P.m;
    Blk := function( supp, mats )
        local M, k;
        M := NullMat( m, m, GF(2) );
        for k in supp do M := M + mats[k]; od;
        return M;
    end;
    LA := List( [ 1 .. nA0 ], i -> List( [ 1 .. nA1 ], j -> Blk( A[i][j], P.Lmats ) ) );
    RB := List( [ 1 .. nB0 ], k -> List( [ 1 .. nB1 ], l -> Blk( B[k][l], P.Rmats ) ) );
    Nq := ( nA1 * nB0 + nA0 * nB1 ) * m;  off := nA1 * nB0 * m;
    HX := NullMat( nA0 * nB0 * m, Nq, GF(2) );  HZ := NullMat( nA1 * nB1 * m, Nq, GF(2) );
    AddBlock := function( M, blk, r0, c0 )
        M{ [ r0+1 .. r0+m ] }{ [ c0+1 .. c0+m ] } := M{ [ r0+1 .. r0+m ] }{ [ c0+1 .. c0+m ] } + blk;
    end;
    for i in [ 1 .. nA0 ] do for k in [ 1 .. nB0 ] do
        rb := ( (i-1) * nB0 + (k-1) ) * m;
        for j in [ 1 .. nA1 ] do AddBlock( HX, LA[i][j], rb, ( (j-1) * nB0 + (k-1) ) * m ); od;
        for l in [ 1 .. nB1 ] do AddBlock( HX, RB[k][l], rb, off + ( (i-1) * nB1 + (l-1) ) * m ); od;
    od; od;
    for j in [ 1 .. nA1 ] do for l in [ 1 .. nB1 ] do
        rb := ( (j-1) * nB1 + (l-1) ) * m;
        for k in [ 1 .. nB0 ] do AddBlock( HZ, TransposedMat( RB[k][l] ), rb, ( (j-1) * nB0 + (k-1) ) * m ); od;
        for i in [ 1 .. nA0 ] do AddBlock( HZ, TransposedMat( LA[i][j] ), rb, off + ( (i-1) * nB1 + (l-1) ) * m ); od;
    od; od;
    return rec( HX := HX, HZ := HZ, N := Nq, LA := LA, RB := RB );
end;

HasCancellation := function( code, A, B, sA, sB, m )
    local i, j, k, l;
    for i in [ 1 .. sA[1] ] do for j in [ 1 .. sA[2] ] do
        if Sum( RowWeights( code.LA[i][j] ) ) <> Length( A[i][j] ) * m then return true; fi;
    od; od;
    for k in [ 1 .. sB[1] ] do for l in [ 1 .. sB[2] ] do
        if Sum( RowWeights( code.RB[k][l] ) ) <> Length( B[k][l] ) * m then return true; fi;
    od; od;
    return false;
end;

CodeWeights := function( HX, HZ )
    local rwX, rwZ, dX, dZ;
    rwX := RowWeights( HX );  rwZ := RowWeights( HZ );
    dX := RowWeights( TransposedMat( HX ) );  dZ := RowWeights( TransposedMat( HZ ) );
    return rec( wmax := Maximum( Maximum( rwX ), Maximum( rwZ ), Maximum( dX + dZ ) ),
                wmin := Minimum( Minimum( rwX ), Minimum( rwZ ) ) );
end;

# full pipeline for one candidate; returns fail or a record
Evaluate := function( P, A, B, sA, sB, wMax, Kmin )
    local code, wts, K, dReq, d, w;
    STATS.cand := STATS.cand + 1;
    if not ForAny( Flat( A ), k -> not P.AinN[k] ) then STATS.inN := STATS.inN + 1; return fail; fi;
    code := CLP_Build( P, A, B, sA, sB );
    if HasCancellation( code, A, B, sA, sB, P.m ) then STATS.cancel := STATS.cancel + 1; return fail; fi;
    wts := CodeWeights( code.HX, code.HZ );
    if wts.wmax > wMax or wts.wmin < MIN_CHECK_W then STATS.weight := STATS.weight + 1; return fail; fi;
    w := wts.wmax;
    K := code.N - RankMat( code.HX ) - RankMat( code.HZ );
    if K > STATS.bestK then STATS.bestK := K; fi;
    if K < Kmin then STATS.K := STATS.K + 1; return fail; fi;
    dReq := RequiredD( code.N, K, w );
    if dReq > D_CAP then STATS.dreq := STATS.dreq + 1; return fail; fi;
    d := DistUB( code.HX, code.HZ, SCREEN_ROUNDS, dReq - 1 );
    if AbsInt( d ) > STATS.bestd then STATS.bestd := AbsInt( d ); fi;
    if d < dReq then STATS.dscreen := STATS.dscreen + 1; return fail; fi;
    d := DistUB( code.HX, code.HZ, FULL_ROUNDS, dReq - 1 );
    if d < dReq then STATS.dfull := STATS.dfull + 1; return fail; fi;
    return rec( code := code, K := K, w := w, d := d, dReq := dReq );
end;

CODE_COUNTER := 0;

LogCode := function( P, T, A, B, res )
    local id, fom, Ael, Bel, tag, key;
    key := [ P.l, P.i, P.m, res.code.N, res.K, res.d, res.w ];
    if key in SEEN_PARAMS then return; fi;
    AddSet( SEEN_PARAMS, key );
    CODE_COUNTER := CODE_COUNTER + 1;  STATS.logged := STATS.logged + 1;
    id := Concatenation( OUT_PREFIX, String( CODE_COUNTER ) );
    fom := Float( res.K * res.d^2 / res.code.N );
    Ael := List( A, row -> List( row, e -> List( e, k -> P.Areps[k] ) ) );
    Bel := List( B, row -> List( row, e -> List( e, k -> P.Lreps[k] ) ) );
    if ForAny( VALIDATE_TARGETS, t -> t[1] = res.code.N and t[2] = res.K ) then tag := "VALIDATION";
    elif not IsDominated( res.code.N, res.K, res.d, res.w ) then tag := "NON-DOMINATED";
    else tag := "FOM"; fi;
    Print( "\n>>> ", tag, "  [[", res.code.N, ",", res.K, ",<=", res.d, "]]  w=", res.w, "  FOM<=", fom,
           "   ", P.desc, "  shapes ", T.sA, "/", T.sB, "\n" );
    AppendTo( RECORD_FILE,
        "\n#### ", id, "  [[", res.code.N, ",", res.K, ",<=", res.d, "]]  w=", res.w, "  FOM_upper=", fom, "  ", tag, "\n",
        "driver := \"", T.label, "\";  group := \"", P.desc, "\";  l := ", P.l, ";  i := ", P.i, ";\n",
        "m := ", P.m, ";  r := ", P.r, ";  lift := ", P.lift, ";  stab := ", P.stab, ";\n",
        "shapeA := ", T.sA, ";  shapeB := ", T.sB, ";\n",
        "A := ", Ael, ";\nB := ", Bel, ";\n",
        "Aidx := ", A, ";  Bidx := ", B, ";   # indices into RightTransversal(G,Core) / RightTransversal(N,H), identity first\n",
        "matrices := \"", id, "_HX.txt / _HZ.txt\";\n" );
    WriteBinaryMatrix( Concatenation( id, "_HX.txt" ), res.code.HX );
    WriteBinaryMatrix( Concatenation( id, "_HZ.txt" ), res.code.HZ );
end;

#############################################################################
##  driver loop
#############################################################################

RunDriver := function( T )
    local pr, P, seen, cnt, key, A, B, res, pairs;
    PrintTo( RECORD_FILE, "# ", T.label, "  -- distances are QDistRnd upper bounds; certify before publishing\n" );
    pairs := Filtered( PAIRS, pr -> pr.l > 0 and T.pairOK( pr ) );
    Sort( pairs, function( a, b )
        if a.atb <> b.atb then return a.atb; fi;
        return a.lift > b.lift;
    end );
    Print( T.label, ":  ", Length( pairs ), " pairs\n" );
    for pr in pairs do
        P := LoadPair( pr );  PrepareActions( P );
        ResetStats();
        Print( "\n== ", P.desc, "  m=", P.m, " r=", P.r, " lift=", P.lift, " stab=", P.stab,
               " |G/Core|=", P.nA, "  N=", T.Nfrom( P.m ), "\n" );
        seen := [ ];  cnt := 0;
        while cnt < CAND_PER_PAIR do
            cnt := cnt + 1;
            A := T.genA( P );  B := T.genB( P );
            key := [ A, B ];
            if key in seen then continue; fi;
            AddSet( seen, key );
            res := Evaluate( P, A, B, T.sA, T.sB, T.wMax, T.Kmin );
            if res <> fail then LogCode( P, T, A, B, res ); fi;
            if cnt mod 2000 = 0 then
                Print( "   [", cnt, "] ", STATS, "\n" );
            fi;
        od;
        Print( "   pair done: ", STATS, "\n" );
    od;
    Print( "\n", T.label, " finished; ", CODE_COUNTER, " codes in ", RECORD_FILE, "\n" );
end;