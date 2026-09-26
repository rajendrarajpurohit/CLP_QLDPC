Read("CLP_core.g");;
Read("pairs_bigm.g");;

CLPReadN := function( f, N )
    local s, d, i, M;
    s := StringFile(f);
    if s = fail then return fail; fi;
    d := Filtered(s, c -> c='0' or c='1');  M := [];
    for i in [0..Length(d)/N-1] do
        Add(M, List(d{[i*N+1..i*N+N]}, c -> (IntChar(c)-48)*One(GF(2))));
    od;
    return M;
end;;

# which indices of MATS sum to BLK ?  greedy: a perm matrix is in the sum
# iff every 1 of that perm matrix is a 1 of BLK (no cancellation in our codes)
MatchBlock := function( blk, mats )
    local hits, i, ok, r, c, tot;
    hits := [];
    for i in [ 1 .. Length(mats) ] do
        ok := true;
        for r in [ 1 .. Length(blk) ] do
            c := Position( mats[i][r], One(GF(2)) );
            if c = fail or IsZero( blk[r][c] ) then ok := false; break; fi;
        od;
        if ok then Add( hits, i ); fi;
    od;
    # verify the sum reproduces blk exactly
    tot := Sum( List( hits, i -> mats[i] ) );
    if tot <> blk then return fail; fi;
    return hits;
end;;

Recover := function( fname, N, lid, iid, mm, label )
    local HX, HZ, K, wts, pr, P, m, off, Ai, B11, B12, B21, B22, Ael, Bel;
    HX := CLPReadN( Concatenation(fname,"_HX.txt"), N );
    HZ := CLPReadN( Concatenation(fname,"_HZ.txt"), N );
    if HX = fail then Print("MISSING ", fname, "\n"); return; fi;
    K := N - RankMat(HX) - RankMat(HZ);
    wts := CodeWeights( HX, HZ );
    pr := First( PAIRS, p -> p.l=lid and p.i=iid and p.m=mm );
    if pr = fail then Print("pair not found\n"); return; fi;
    P := LoadPair( pr );  PrepareActions( P );
    m := P.m;  off := 2*m;

    Ai  := MatchBlock( HX{[1..m]}{[1..m]},                   P.Lmats );
    B11 := MatchBlock( HX{[1..m]}{[off+1..off+m]},           P.Rmats );
    B12 := MatchBlock( HX{[1..m]}{[off+m+1..off+2*m]},       P.Rmats );
    B21 := MatchBlock( HX{[m+1..2*m]}{[off+1..off+m]},       P.Rmats );
    B22 := MatchBlock( HX{[m+1..2*m]}{[off+m+1..off+2*m]},   P.Rmats );

    Print("\n========== ", label, " ==========\n");
    Print("file: ", fname, "\n");
    Print("[[", N, ",", K, "]]  w=", wts.wmax, "\n");
    if fail in [Ai,B11,B12,B21,B22] then
        Print("  RECOVERY FAILED (block did not match)\n");
        Print("  Ai=",Ai," B11=",B11," B12=",B12," B21=",B21," B22=",B22,"\n");
        return;
    fi;
    Ael := List( Ai,  i -> P.Areps[i] );
    Bel := [ [ List(B11,i->P.Lreps[i]), List(B12,i->P.Lreps[i]) ],
             [ List(B21,i->P.Lreps[i]), List(B22,i->P.Lreps[i]) ] ];
    Print("G = SmallGroup(",lid,",",iid,") = ", StructureDescription(P.G), "\n");
    Print("H = ", StructureDescription(P.H), ", |H|=", Size(P.H),
          ", normal=", IsNormal(P.G,P.H), "\n");
    Print("N_G(H) = ", StructureDescription(P.N), ", |N|=", Size(P.N), "\n");
    Print("Core = ", Size(P.core), ", m=",P.m," r=",P.r," lift=",P.lift," stab=",P.stab,"\n");
    Print("a = ", Ael, "\n");
    Print("    orders: ", List(Ael,Order), "\n");
    Print("b = ", Bel, "\n");
    Print("    orders: ", List(Bel, row -> List(row, e -> List(e,Order))), "\n");
    Print("Aidx(this session) = ", [[Ai]], "\n");
    Print("Bidx(this session) = ", [[B11,B12],[B21,B22]], "\n");
end;;

Recover( "G1_n360k16d22", 360, 270, 25, 90, "G1 [[360,16,<=22]] w=7" );
Recover( "G6_n360k16d20", 360, 270, 25, 90, "G6 [[360,16,<=20]] w=6" );
Recover( "G2_n360k16d21", 360, 270, 25, 90, "G2 [[360,16,<=21]] w=7" );
Recover( "G3_n360k12d22", 360, 270, 25, 90, "G3 [[360,12,<=22]] w=6" );
Recover( "G4_n360k12d22", 360, 270, 25, 90, "G4 [[360,12,<=22]] w=6" );
Recover( "G5_n360k12d22", 360, 270, 25, 90, "G5 [[360,12,<=22]] w=6" );
Recover( "G7_n360k12d22", 360, 270, 26, 90, "G7 [[360,12,<=22]] w=6" );
Recover( "G8_n360k12d24", 360, 270, 25, 90, "G8 [[360,12,<=24]] w=?" );
QUIT;
