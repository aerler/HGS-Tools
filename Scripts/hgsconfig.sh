#!/bin/bash
# this is a sample configuration file for HGS synchronization scripts

HGS_DIR="${HGS_DIR:-"HGS"}" # HGS folder name (same on source/dest)
PROJECT="${PROJECT:-""}" # project subfolder in HGS folder 
# remote host, i.e. where the data is
HOST="${HOST:-""}" # can be defined in ssh config
HOST_ROOT="${HOST_ROOT:-""}" # remote DATA_ROOT
# local machine
DATA_ROOT="${DATA_ROOT:-""}" # root folder of data repository
# general setings
SSH="-o BatchMode=yes -o ControlPath=${HOME}/master-%l-%r@%h:%p -o ControlMaster=auto -o ControlPersist=1" # default SSH login options
NICENESS=${NICENESS:-0} # process priority (0 is standard)