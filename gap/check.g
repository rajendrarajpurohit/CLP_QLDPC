Read("CLP_core.g");;
Read("NOVA_lib.g");;
Read("pairs_bigm.g");;
pr := First( PAIRS, p -> p.l=270 and p.i=26 and p.m=90 );;
P := LoadPair( pr );;  PrepareActions( P );;
Print("pair: ", P.desc, " m=", P.m, " nA=", P.nA, " nL=", P.nL, "\n");
nb := 0;; nd := 0;; nt := 0;; nk := 0;; ok := 0;;
ws := [];;
for t in [1..1000] do
  A := arrA( P.nA, 1, 1 );;
  B := arrB( P.nL, 2, 2 );;
  if ForAny( Flat(A), x -> not P.AinN[x] ) then
    c := CLP_Build( P, A, B, [1,1], [2,2] );;
    nb := nb + 1;
    if dupcol( c.HX, c.HZ ) then nd := nd + 1;
    elif thincol( c.HX, c.HZ ) then nt := nt + 1;
    else
      w := cw( c.HX, c.HZ );;
      Add( ws, w );
      K := c.N - RankMat(c.HX) - RankMat(c.HZ);;
      if K < kfloor(w) then nk := nk + 1; else ok := ok + 1; fi;
    fi;
  fi;
od;
Print("built ", nb, "  dupcol ", nd, "  thincol ", nt, "  belowk ", nk, "  survive ", ok, "\n");
Print("weight histogram: ", Collected( ws ), "\n");
rdN := function( f, N )
    local s, x, i, M;
    s := StringFile(f);  x := Filtered(s, c -> c='0' or c='1');  M := [];
    for i in [0..Length(x)/N-1] do
        Add(M, List(x{[i*N+1..i*N+N]}, c -> (IntChar(c)-48)*One(GF(2))));
    od;
    return M;
end;;
HX := rdN("CONTROL_HX.txt",360);;  HZ := rdN("CONTROL_HZ.txt",360);;
Print("CONTROL: k=", 360-RankMat(HX)-RankMat(HZ), " w=", cw(HX,HZ),
      " CSS=", IsZero(HX*TransposedMat(HZ)),
      " dupcol=", dupcol(HX,HZ), " thincol=", thincol(HX,HZ), "\n");
QUIT;
