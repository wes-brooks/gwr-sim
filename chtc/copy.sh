#! /bin/sh

for ((i=1; i<=1800; i++))
do
  cp /home/wbrooks2/gwr-sim/chtc/simoutput-gaussian/$i/CoefEstimates* /home/wbrooks2/gwr-sim/output/gaussian/
  cp /home/wbrooks2/gwr-sim/chtc/simoutput-gaussian/$i/Data* /home/wbrooks2/gwr-sim/output/gaussian/
  cp /home/wbrooks2/gwr-sim/chtc/simoutput-gaussian/$i/MiscParams* /home/wbrooks2/gwr-sim/output/gaussian/
done
