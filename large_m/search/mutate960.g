# usage: echo 'SEED:=1;; STEPS:=400;; Read("mutate960.g");' | gap -q
Read("CLP_core.g");;
Reset( GlobalMersenneTwister, 9600 + SEED );;
TAG := Concatenation( "mut960_s", String(SEED) );;  RECORD := Concatenation( TAG, "_results.txt" );;
DCUT := 36;;  KMIN := 20;;  COUNT := 0;;
pr := rec( l := 480, i := 1005, pc := true, Hdesc := "C2", Hgens := [ [0,0,1,0,0,0,0] ] );;
P := LoadPair( pr );;  PrepareActions( P );;
A := [ [ [ 1, 130, 422 ] ] ];;  B := [ [ [ 1, 85 ], [ 1 ] ], [ [ 1 ], [ 19, 42 ] ] ];;
Fresh := function( n, avoid ) local x; repeat x := Random( [ 2 .. n ] ); until not x in avoid; return x; end;;
Mutate := function( A, B )
  local A2, B2, r, e, j, k, l, heavy;
  A2 := StructuralCopy( A );  B2 := StructuralCopy( B );
  r := Random( [ 1 .. 10 ] );
  if r <= 4 then
    j := Random( [ 2, 3 ] );  A2[1][1][j] := Fresh( P.nA, A2[1][1] );
  elif r <= 9 then
    heavy := Filtered( [ [1,1],[1,2],[2,1],[2,2] ], e -> Length( B2[e[1]][e[2]] ) = 2 );
    e := Random( heavy );  B2[e[1]][e[2]][2] := Fresh( P.nL, B2[e[1]][e[2]] );
  else
    if Length( B2[1][1] ) = 2 then
      B2 := [ [ [1], [1, Fresh(P.nL,[1])] ], [ [1, Fresh(P.nL,[1])], [1] ] ];
    else
      B2 := [ [ [1, Fresh(P.nL,[1])], [1] ], [ [1], [1, Fresh(P.nL,[1])] ] ];
    fi;
  fi;
  return [ A2, B2 ];
end;;
Score := function( A, B )
  local c, K, wts, d;
  if not ForAny( Flat(A), x -> not P.AinN[x] ) then return fail; fi;
  c := CLP_Build( P, A, B, [1,1], [2,2] );
  if HasCancellation( c, A, B, [1,1], [2,2], P.m ) then return fail; fi;
  K := c.N - RankMat(c.HX) - RankMat(c.HZ);
  if K < KMIN then return fail; fi;
  wts := CodeWeights( c.HX, c.HZ );
  if wts.wmax > 6 then return fail; fi;
  d := DistUB( c.HX, c.HZ, 500, DCUT );
  if d < DCUT then return fail; fi;
  d := DistUB( c.HX, c.HZ, 20000, DCUT );
  if d < DCUT then return fail; fi;
  return rec( c := c, K := K, d := d );
end;;
cur := Score( A, B );;
Print( TAG, " center k=", cur.K, " d20k<=", cur.d, "\n" );
curK := cur.K;;
for step in [ 1 .. STEPS ] do
  m := Mutate( A, B );;
  s := Score( m[1], m[2] );;
  if s = fail then continue; fi;
  Print( "  step ", step, " k=", s.K, " d20k<=", s.d, "\n" );
  dZ := DistRandCSS( s.c.HX, s.c.HZ, 200000, DCUT : field := GF(2) );;
  dX := DistRandCSS( s.c.HZ, s.c.HX, 200000, DCUT : field := GF(2) );;
  if dZ >= DCUT and dX >= DCUT then
    COUNT := COUNT + 1;  d := Minimum( dZ, dX );;
    id := Concatenation( TAG, "_", String(COUNT), "_n960k", String(s.K), "d", String(d), "w6" );;
    Print( ">>> [[960,", s.K, ",<=", d, "]] w=6 dZ<=", dZ, " dX<=", dX, " Q=", Float(s.K*d^2/960), "\n" );
    AppendTo( RECORD, "\n#### ", id, "  [[960,", s.K, ",<=", d, "]]  w=6  dZ<=", dZ, "  dX<=", dX, "  Q=", Float(s.K*d^2/960), "\n",
      "group := \"", P.desc, "\";  l := 480;  i := 1005;  m := 240;\nHgens := ", pr.Hgens, ";\n",
      "A := ", List(m[1], r -> List(r, e -> List(e, x -> P.Areps[x]))), ";\n",
      "B := ", List(m[2], r -> List(r, e -> List(e, x -> P.Lreps[x]))), ";\n",
      "Aidx := ", m[1], ";  Bidx := ", m[2], ";\n" );
    WriteBinaryMatrix( Concatenation( id, "_HX.txt" ), s.c.HX );
    WriteBinaryMatrix( Concatenation( id, "_HZ.txt" ), s.c.HZ );
    if s.K >= curK then A := m[1];  B := m[2];  curK := s.K;  Print( "  moved center to k=", curK, "\n" ); fi;
  fi;
od;
Print( TAG, " done, ", COUNT, " codes logged\n" );
QUIT;
