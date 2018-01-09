#!/bin/bash

BIN=$@

echo
echo "Initializing New EnKF Assimilation Exepriment"
echo
# get total number of time steps
echo
N=$( sed -n '/total_number_of_time_steps/ s/^[[:space:]]*total_number_of_time_steps[[:space:]]*=[[:space:]]*\(.*\)$/\1/p' EnKFparameters.ini )
echo "Total Number of Time-steps: $N"
# clean a bit
echo "Cleaning Run Directory"
rm -rf proc_*/ out/ 
mkdir out/
echo
# switch restart off
sed -i  '/restart_flag/ s/^\([[:space:]]*restart_flag[[:space:]]*=[[:space:]]\)*.*$/\10/g' EnKFparameters.ini
# run EnKF (first round)
echo
echo "Starting EnKF (first attempt)"
echo $BIN
$BIN
echo

# read current time step
CN=$( cat backup.info )

while [ $CN -lt $N ]; do
  echo
  echo "EnKF assimilation step failed; attempting restart..."
  # set restart flag
  sed -i  '/restart_flag/ s/^\([[:space:]]*restart_flag[[:space:]]*=[[:space:]]\)*.*$/\11/g' EnKFparameters.ini
  # run again
  echo $BIN
  $BIN
  # read time step again
  CN=$( cat backup.info )
  echo
done # while CN < N

# complete
echo
echo "EnKF Assimilation Experiment completed!!!"
echo
