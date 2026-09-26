Read("CLP_core.g");;
Read("NOVA_lib.g");;
SEED := 2102;;  Reset( GlobalMersenneTwister, SEED );;
FAM := "G2";;
OUT := "G2_results.txt";;
PrintTo( OUT, "# l i m n k w d Q ratio\n" );;
NCAND := 4000;;  NPAIR := 12;;
NPM := 6;;  A0 := 1;;  A1 := 1;;  B0 := 3;;  B1 := 3;;
MLO := 34;;  MHI := 66;;
COUNT := 0;;

Read("pairs_final.g");;  R1 := Filtered(PAIRS, p->p.l>0);;
Read("pairs_mid.g");;    R2 := Filtered(PAIRS, p->p.l>0);;
Read("pairs_bigm.g");;   R3 := Filtered(PAIRS, p->p.l>0);;
ALL := Concatenation( R1, R2, R3 );;
Print( "pairs: ", Length(ALL), "  shape ", [A0,A1], "/", [B0,B1], "  n=", NPM, "m\n" );

for mm in [ MLO .. MHI ] do
  sel := Filtered( ALL, p -> p.m = mm and mm*NPM >= 60 and mm*NPM <= 400 );
  if Length(sel) = 0 then continue; fi;
  if Length(sel) > NPAIR then sel := Shuffle( ShallowCopy(sel) ){[1..NPAIR]}; fi;
  Print( "== m=", mm, " n=", mm*NPM, " pairs=", Length(sel), "\n" );
  for pr in sel do
    P := LoadPair( pr );;  PrepareActions( P );;
    best := 0.0;;
    for t in [1..NCAND] do
      A := arrA( P.nA, A0, A1 );;
      B := arrB( P.nL, B0, B1 );;
      if ForAny( Flat(A), x -> not P.AinN[x] ) then
        c := CLP_Build( P, A, B, [A0,A1], [B0,B1] );;
        if not dupcol(c.HX,c.HZ) and not thincol(c.HX,c.HZ) then
          w := cw( c.HX, c.HZ );;
          if w >= 5 and w <= 10 then
            K := c.N - RankMat(c.HX) - RankMat(c.HZ);;
            if K >= kfloor(w) then
              d := AbsInt( DistRandCSS( c.HX, c.HZ, 500, 2 : field := GF(2) ) );;
              q := Float( K*d^2/c.N );;
              if q > best then best := q; fi;
              if q > qbar(w) then
                d := AbsInt( DistRandCSS( c.HX, c.HZ, 20000, 2 : field := GF(2) ) );;
                q := Float( K*d^2/c.N );;
                if q > qbar(w) then
                  d := AbsInt( DistRandCSS( c.HX, c.HZ, 200000, 2 : field := GF(2) ) );;
                  q := Float( K*d^2/c.N );;
                  if q > qbar(w) then
                    COUNT := COUNT + 1;
                    AppendTo( OUT, pr.l," ",pr.i," ",P.m," ",c.N," ",K," ",w," ",d," ",q," ",Float(d/Sqrt(Float(c.N))),"\n" );
                    AppendTo( OUT, "#  A := ", List(A, r -> List(r, e -> List(e, x -> P.Areps[x]))), "\n" );
                    AppendTo( OUT, "#  B := ", List(B, r -> List(r, e -> List(e, x -> P.Lreps[x]))), "\n" );
                    AppendTo( OUT, "#  Aidx := ", A, "  Bidx := ", B, "\n" );
                    AppendTo( OUT, "#  Hgens := ", pr.Hgens, "\n" );
                    Print( "  >>> [[", c.N, ",", K, ",<=", d, "]] w=", w, " Q=", q, "  SG(", pr.l, ",", pr.i, ")\n" );
                    WriteBinaryMatrix( Concatenation(FAM,"_",String(COUNT),"_n",String(c.N),"k",String(K),"d",String(d),"w",String(w),"_HX.txt"), c.HX );
                    WriteBinaryMatrix( Concatenation(FAM,"_",String(COUNT),"_n",String(c.N),"k",String(K),"d",String(d),"w",String(w),"_HZ.txt"), c.HZ );
                  fi;
                fi;
              fi;
            fi;
          fi;
        fi;
      fi;
    od;
    Print( "   SG(", pr.l, ",", pr.i, ") bestQ=", best, "\n" );
  od;
od;
Print( "\n", FAM, " done, ", COUNT, " codes\n" );
QUIT;QUIT;
