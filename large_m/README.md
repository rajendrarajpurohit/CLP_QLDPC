# CLP codes at 480 <= n <= 1000, check weight 6

Post paper extension of the search in arXiv:2610.12782 to 120 <= m <= 250, run 24 Oct to 2 Nov 2026.
Same construction (gap/CLP_core.g of the main repository), same 1x1 over 2x2 layout, |a| = 3,
B entry weights 2 1 / 1 2 or 1 2 / 2 1, so every code has check weight 6. Pairs (G,H): H cyclic of
order 2 to 9, core free, [G:N_G(H)] in {2,3,4}, [N_G(H):H] >= 4, 25,240 pairs enumerated by
search/bigpairs.g over SmallGroups of order 240 to 2000 (orders with more than 2000 groups skipped).

All distances are two sided QDistRnd upper bounds with at least 10^6 rounds. The exact round count
of each code is in its screen.txt. None is certified exact; at these distances the integer program
certification of the paper is out of reach.

| folder | code | Q | group | H | rounds |
|---|---|---|---|---|---|
| 576_16_30 | [[576,16,<=30]] | 25.0 | SmallGroup(432,692) | C3 | 1000000 |
| 576_20_24 | [[576,20,<=24]] | 20.0 | SmallGroup(288,1021) | C2 | 1000000 |
| 588_18_28 | [[588,18,<=28]] | 24.0 | SmallGroup(1029,15) | C7 | 1000000 |
| 640_16_32 | [[640,16,<=32]] | 25.6 | SmallGroup(800,639) | C5 | 1000000 |
| 640_24_24 | [[640,24,<=24]] | 21.6 | SmallGroup(800,453) | C5 | 1000000 |
| 672_24_26 | [[672,24,<=26]] | 24.1 | SmallGroup(1176,269) | C7 | 1000000 |
| 720_16_36 | [[720,16,<=36]] | 28.8 | SmallGroup(360,116) | C2 | 1000000 |
| 840_16_40 | [[840,16,<=40]] | 30.5 | SmallGroup(1470,44) | C7 | 1000000 |
| 864_16_44 | [[864,16,<=44]] | 35.9 | SmallGroup(432,525) | C2 | 1000000 |
| 864_20_40 | [[864,20,<=40]] | 37.0 | SmallGroup(648,720) | C3 | 1000000 |
| 960_16_44 | [[960,16,<=44]] | 32.3 | SmallGroup(480,1065) | C2 | 10000000 |
| 960_16_46 | [[960,16,<=46]] | 35.3 | SmallGroup(1200,660) | C5 | 10000000 |
| 960_16_46_G480_671 | [[960,16,<=46]] | 35.3 | SmallGroup(480,671) | C2 | 10000000 |
| 960_16_50 | [[960,16,<=50]] | 41.7 | SmallGroup(1200,229) | C5 | 1000000 |
| 960_20_46 | [[960,20,<=46]] | 44.1 | SmallGroup(480,1005) | C2 | 6000000 |
| 960_20_50 | [[960,20,<=50]] | 52.1 | SmallGroup(480,1005) | C2 | 1000000 |
| 1000_16_58 | [[1000,16,<=58]] | 53.8 | SmallGroup(1250,52) | C5 | 1000000 |

Each folder: HX.txt, HZ.txt (0/1 rows), RECORD.txt (the search record with G, H, a, B as words and as
indices), REBUILD.txt (how to rebuild with gap/rebuild.g), screen.txt (every QDistRnd run on the code),
and board_witnesses.json where a submission draft for qecchallenge.com was built (the witness vectors
their tool found, verified by their verifier).

search/ holds every script, pair list and raw results file of the extension, in the order
bigpairs.g -> merge_pairs.g -> trim.g, trim_rest.g -> W6big.g, W6rest.g (sweep, bar Q>30) ->
focus.g, W6focus.g (n 560 to 700, Pareto bar needd6.g) -> mutate960.g (walk around [[960,20]]) ->
W6shape.g (|a| = 2 and 4 shapes, no code found) -> rescreen.g (screens).

Negative result: the |a| = 2 and |a| = 4 shapes give k from 24 to over 100 at these lengths but
over 9,560 pairs no draw paired that k with a distance above the board frontier for its k.
