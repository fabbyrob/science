#$1 = file for lazyloop
#$2 = prefix
#$3 = output directory for simulations
python ~/bin/lazyLoopGeneral.py $1 > $2.loop.out 2> $2.loop.err
python ~/bin/addToDB.py $2.TEdb $3/*.txt > $2.db.out 2> $2.db.err
sqlite3 $2.TEdb < ~/repos/trunk/TE\ simulation/getTEdata.sql 2> $2.sql.err
Rscript ~/bin/R/te_plotter_fun.R allSims.csv $2 > $2.R.out 2> $2.R.err 