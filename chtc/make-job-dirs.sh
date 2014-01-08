#! /bin/sh

for ((i=1; i<=$2; i++))
do
  mkdir -p simdata-$1/$i
  mkdir -p simdata-$1/$i/output
  echo "$1\n$i" > simdata-$1/$i/jobid.txt
done
