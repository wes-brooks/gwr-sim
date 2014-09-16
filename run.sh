#! /bin/sh
for i in {1..18}
do
    Rscript R/simulation-gaussian.r $i
done
