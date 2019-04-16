#!/bin/bash
# a mock executable that simulates HGS output for queue submission tests
# this script also supports restart output; the numebr of restats if controlled
# by the parameter C

# read previous restart info file and count one down
if [ -e "$RSTINFO" ]; then
    C=$( cat "$RSTINFO" )
    C=$(( $C - 1 ))
fi # RSTINFO
# wait a little and write a dummy output file
sleep 5
echo $C > hgs_test.$C
# handle restart process
if [ $C -eq 0 ]; then 
    # terminate restart cycle if C is 0
    echo "Normal Exit"
else 
    # if C > 0, write restartfile
    echo $C > "$RSTINFO"
    echo $C > "$RSTFILE"
    echo $C
fi # C
