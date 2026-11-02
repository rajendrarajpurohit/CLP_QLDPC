Read("CLP_core.g");;
W := 6;;  NT := 2000;;  QTH := 30.0;;

TAG := Concatenation( "probe6_", String(ORD) );;
RECORD := Concatenation( TAG, "_results.txt" );;  COUNT := 0;;
Reset( GlobalMersenneTwister, 6000 + ORD );;
NeedD := function( K, N ) local dd; for dd in [ 8 .. 60 ] do if Float(K*dd^2/N) > QTH then return dd; fi; od; return 99; end;;
Bsum := function( B )
  return Maximum( Concatenation( List( B, r -> Sum( List( r, Length ) ) ), List( [1,2], j -> Sum( List( B, r -> Length( r[j] ) ) ) ) ) );
end;;
D18 := First( AllSmallGroups( 18 ), g -> StructureDescription( g ) = "(C3 x C3) : C2" );;
id := IdGroup( DirectProduct( CyclicGroup( ORD / 18 ), D18 ) );;
G := SmallGroup( id[1], id[2] );;  pc := Pcgs( G );;
Print( "G = SmallGroup(", id[1], ",", id[2], ") = ", StructureDescription( G ), "\n" );
cls := Filtered( ConjugacyClassesSubgroups( G ), c -> Size( Representative( c ) ) = 3 );;
Hs := Filtered( List( cls, Representative ), H -> Size( Core( G, H ) ) = 1 and Index( G, Normalizer( G, H ) ) in [ 2, 3, 4 ] );;
Print( "core free C3 classes with [G:N] in {2,3,4}: ", Length( Hs ), "\n" );
for H in Hs do
  pr := rec( l := id[1], i := id[2], pc := true, Hdesc := "C3",
             Hgens := List( GeneratorsOfGroup( H ), h -> ExponentsOfPcElement( pc, h ) ) );;
  P := LoadPair( pr );;  PrepareActions( P );;
  if P.lift < 4 then continue; fi;
  nK := 0;;  nD := 0;;  best500 := 0;;  bestK := 0;;
  Print( "== ", P.desc, " Hgens=", pr.Hgens, " m=", P.m, " r=", P.r, " lift=", P.lift, " n=", 4*P.m, "\n" );
  for t in [ 1 .. NT ] do
    wa := 3;;
    A := [ [ SupportWithOne( P.nA, wa ) ] ];;
    B := [ [ SupportWithOne(P.nL,Random([1,2])), SupportWithOne(P.nL,Random([1,2])) ],
           [ SupportWithOne(P.nL,Random([1,2])), RandomEntry(P.nL,1,2) ] ];;
    if Bsum( B ) <= W - wa and ForAny( Flat(A), x -> not P.AinN[x] ) then
      c := CLP_Build( P, A, B, [1,1], [2,2] );;
      if not HasCancellation( c, A, B, [1,1], [2,2], P.m ) then
        K := c.N - RankMat(c.HX) - RankMat(c.HZ);;
        if K > bestK then bestK := K; fi;
        if K >= 16 then
          nK := nK + 1;
          wts := CodeWeights( c.HX, c.HZ );;
          if wts.wmax <= W then
            need := NeedD( K, c.N );;
            if DistUB( c.HX, c.HZ, 50, 4 ) > 0 then
              d := DistUB( c.HX, c.HZ, 500, 1 );;
              if d > best500 then best500 := d; Print( "   t=", t, " K=", K, " d500<=", d, " need ", need, "\n" ); fi;
              if d >= need then
                nD := nD + 1;
                d := DistUB( c.HX, c.HZ, 20000, need-1 );;
                if d >= need then
                  dZ := DistRandCSS( c.HX, c.HZ, 200000, need-1 : field := GF(2) );;
                  dX := DistRandCSS( c.HZ, c.HX, 200000, need-1 : field := GF(2) );;
                  if dZ >= need and dX >= need then
                    COUNT := COUNT + 1;  d := Minimum( dZ, dX );;
                    cid := Concatenation( TAG, "_", String(COUNT), "_n", String(c.N), "k", String(K), "d", String(d), "w", String(wts.wmax) );;
                    Print( "\n>>> [[", c.N, ",", K, ",<=", d, "]] w=", wts.wmax, " dZ<=", dZ, " dX<=", dX, " Q=", Float(K*d^2/c.N), "  ", P.desc, "\n" );
                    AppendTo( RECORD, "\n#### ", cid, "  [[", c.N, ",", K, ",<=", d, "]]  w=", wts.wmax, "  dZ<=", dZ, "  dX<=", dX, "  Q=", Float(K*d^2/c.N), "\n",
                      "group := \"", P.desc, "\";  l := ", P.l, ";  i := ", P.i, ";  m := ", P.m, ";  lift := ", P.lift, ";  stab := ", P.stab, ";\n",
                      "Hgens := ", pr.Hgens, ";\n",
                      "A := ", List(A, r -> List(r, e -> List(e, x -> P.Areps[x]))), ";\n",
                      "B := ", List(B, r -> List(r, e -> List(e, x -> P.Lreps[x]))), ";\n",
                      "Aidx := ", A, ";  Bidx := ", B, ";\n" );
                    WriteBinaryMatrix( Concatenation( cid, "_HX.txt" ), c.HX );
                    WriteBinaryMatrix( Concatenation( cid, "_HZ.txt" ), c.HZ );
                  fi;
                fi;
              fi;
            fi;
          fi;
        fi;
      fi;
    fi;
  od;
  Print( "   pair done: bestK=", bestK, "  K>=16: ", nK, "  best d500 ", best500, "  passed 20000: ", nD, "  logged: ", COUNT, "\n" );
od;
Print( "\n", TAG, " done, ", COUNT, " codes\n" );
QUIT;
