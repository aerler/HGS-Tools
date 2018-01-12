#!/bin/bash

BIN=$@

if [[ -n "$WEXP" ]]; then
  echo
  echo "Waiting for previous experiment to finish):"
  echo "  ${WEXP}"
  while [ ! -f "${WEXP}/COMPLETED" ]; do sleep 100; done
  echo "Experiment finished - starting new experiment!"
  echo
fi # wait...

echo
echo "Initializing New EnKF Assimilation Experiment"
echo
# get total number of time steps
echo
N=$( sed -n '/total_number_of_time_steps/ s/^[[:space:]]*total_number_of_time_steps[[:space:]]*=[[:space:]]*\(.*\)$/\1/p' EnKFparameters.ini )
echo "Total Number of Time-steps: $N"
# clean a bit
if [ -e out/ ]; then
  echo "Cleaning Run Directory (moving existing output to backup)"
  rm -rf out_backup
  mv out out_backup
  [ -f backup.info ] && mv backup.info out_backup/
  [ -f enkf_*.log ] && mv enkf_*.log out_backup/
  [ -f enkf_*.log.gz ] && mv enkf_*.log.gz out_backup/
else
  echo "Cleaning Run Directory"
  rm -f backup.info enkf_*.log enkf_*.log.gz
fi # -e out/
rm -rf proc_*/ tmp.log
mkdir out/
echo
# switch restart off
sed -i  '/restart_flag/ s/^\([[:space:]]*restart_flag[[:space:]]*=[[:space:]]\)*.*$/\10/g' EnKFparameters.ini
# run EnKF (first round)
echo
echo "Starting EnKF (first attempt)"
echo $BIN
echo
$BIN &> tmp.log 
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
  $BIN &> tmp.log
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
touch COMPLETED
echo
echo "EnKF Assimilation Experiment completed!!!"
echo
