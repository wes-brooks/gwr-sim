#! /bin/sh

for ((i=1; i<=$1; i++))
do
  mkdir -p simdata-gaussian/$i
  echo "$2\n$i" > simdata-gaussian/$i/jobid.txt
done
