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
rm -rf proc_*/ out/ backup.info enkf_*.log enkf_*.log.gz tmp.log
mkdir out/
echo
# switch restart off
sed -i  '/restart_flag/ s/^\([[:space:]]*restart_flag[[:space:]]*=[[:space:]]\)*.*$/\10/g' EnKFparameters.ini
# run EnKF (first round)
echo
echo "Starting EnKF (first attempt)"
echo $BIN
echo
$BIN > tmp.log 
echo

# check completion
if [ -f backup.info ]; then
  echo "  ..initialization successful!"
else
  echo
  echo "  >>>  initialization failed  ---  aborting!!!"
  mv tmp.log enkf_failed.log
  echo "  (see log file 'enkf_failed.log' for details)" 
  echo
  exit 1
fi # backup.info



# read current time step
CN=$( cat backup.info )
LOG="enkf_${CN}.log"
echo "Completed $CN times steps; compressing log file to '$LOG'"
mv tmp.log $LOG
gzip $LOG
echo

while [ $CN -lt $N ]; do
  echo
  echo "EnKF assimilation did not complete; attempting restart..."
  # set restart flag
  sed -i  '/restart_flag/ s/^\([[:space:]]*restart_flag[[:space:]]*=[[:space:]]\)*.*$/\11/g' EnKFparameters.ini
  # run again
  echo $BIN
  echo
  $BIN > tmp.log
  echo
  # read time step again
  CN=$( cat backup.info )
  LOG="enkf_${CN}.log"
  echo "Completed $CN times steps; compressing log file to '$LOG'"
  mv tmp.log $LOG
  gzip $LOG
  echo
done # while CN < N

# complete
echo
echo "EnKF Assimilation Experiment completed!!!"
echo
