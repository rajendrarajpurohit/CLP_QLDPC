CLPReadN := function( f, N )
    local s, d, i, M;
    s := StringFile(f);  d := Filtered(s, c -> c='0' or c='1');  M := [];
    for i in [0..Length(d)/N-1] do
        Add(M, List(d{[i*N+1..i*N+N]}, c -> (IntChar(c)-48)*One(GF(2))));
    od;
    return M;
end;;
HX := CLPReadN("G1_n360k16d22_HX.txt", 360);;
HZ := CLPReadN("G1_n360k16d22_HZ.txt", 360);;
PrintTo("G1_support.txt", "# G1 [[360,16,<=22]] w=7\n# HX: row -> qubit indices (1-based)\n");
for i in [1..Length(HX)] do
    AppendTo("G1_support.txt", "HX ", i, ": ", Positions(HX[i], One(GF(2))), "\n");
od;
AppendTo("G1_support.txt", "# HZ\n");
for i in [1..Length(HZ)] do
    AppendTo("G1_support.txt", "HZ ", i, ": ", Positions(HZ[i], One(GF(2))), "\n");
od;
Print("wrote G1_support.txt\n");
QUIT;
