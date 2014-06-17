.mode csv
.output allSims.csv
select ID, a.teNum, l, N, G, u, x, calc, t, cross, s, r, K, sel, c, b, e, z, v, generation, class, fitness, output.teNum, silence, silencediff, ne from (select * from simulation) as a left join output on ID = simID;

.mode csv
.output lastGen.csv
select ID, a.teNum, l, N, G, u, x, calc, t, cross, s, r, K, sel, c, b, e, z, v, generation, class, fitness, output.teNum, silence, silencediff, ne from (select * from simulation) as a left join output on ID = simID where (generation = 999 or generation = 4999);

.mode csv
.output finalFreqs.csv
select ID, a.teNum, l, N, G, u, x, calc, t, cross, s, r, K, sel, c, b, e, z, v, generation, freq from (select * from simulation) as a left join finalFreqs on ID = simID ;
