#!/bin/bash
# Bash script to pull/update all HG and SVN repositories
# Andre R. Erler, 2013, GPL v3, revised 30/11/2014

# pre-process arguments using getopt
if [ -z $( getopt -T ) ]; then
  TMP=$( getopt -o d:f:tr:gqspuh --long hgs-root:,hgs-tools:,no-hgs-tools,recurse:,no-git,hg,svn,pull-only,update-only,help -n "$0" -- "$@" ) # pre-process arguments
  [ $? != 0 ] && exit 1 # getopt already prints an error message
  eval set -- "$TMP" # reset positional parameters (arguments) to $TMP list
fi # check if GNU getopt ("enhanced")
# set default parameters
HGSROOT="$HGS_ROOT"
if [ -d "$CODE_ROOT/HGS Tools/" ]; then HGSTOOLS="$CODE_ROOT/HGS Tools/"
elif [ -d "$CODE_ROOT/HGS-Tools/" ]; then HGSTOOLS="$CODE_ROOT/HGS-Tools/"
else HGSTOOLS=''; fi
HGSSYNC=1
RECLEV=1
GITSYNC=1
HGSYNC=0
SVNSYNC=0
HGPULL=1
HGUPDATE=1
# parse arguments
#while getopts 'fs' OPTION; do # getopts version... supports only short options
while true; do
  case "$1" in
    -d | --hgs-root     )   HGSROOT=$2; shift 2;;
    -f | --hgs-tools    )   HGSTOOLS=$2; shift 2;;
    -t | --no-hgs-tools )   HGSSYNC=0; shift;;
    -r | --recurse      )   RECLEV=$2; shift 2;;
    -g | --no-git       )   GITSYNC=0; shift;;
    -q | --hg           )   HGSYNC=1; shift;;
    -s | --svn          )   SVNSYNC=1; shift;;
    -p | --pull-only    )   HGUPDATE=0; shift;;
    -u | --update-only  )   HGPULL=0; shift;;
    -h | --help         )   echo -e " \
                            \n\
    -d | --hgs-root      Root folder for HGS Model repositories (default: $HGSROOT) \n\
    -f | --hgs-tools     Root folder for HGS Tools repositories (default: $HGSTOOLS) \n\
    -t | --no-hgs-tools  Do not update the HGS Tools repository \n\
    -r | --recurse       Set a maximum level for recursing into sub-folders (default: 1; max: 3) \n\
    -g | --no-git        Ignore Git repositories \n\
    -q | --hg            Also scan and synchronize HG repositories \n\
    -s | --svn           Also scan and synchronize SVN repositories \n\
    -p | --pull-only     Only run HG/Git pull (no updates) \n\
    -u | --update-only   Only run HG update (no pull) \n\
    -h | --help          print this help \n\
                             "; exit 0;; # \n\ == 'line break, next line'; for syntax highlighting
    -- ) shift; break;; # this terminates the argument list, if GNU getopt is used
    * ) break;;
  esac # case $@
done # while getopts  


## update functions
function UPDATE () {
    local F="$1" # repo folder (with hidden repo subfolder)
    local LEC=0 # local error counter
    local D=${F%/.*} # get parent of .hg folder
    echo "$D" # ${D#"${HGSROOT}/"}
    cd "${D}"
    if [[ "$F" == *'/.git' ]]; then
	    # pull & update repository
	    if [ ${HGPULL} -eq 1 ] && [ ${HGUPDATE} -eq 1 ]; then 
	      git pull
	      [ $? -gt 0 ] && LEC=$(( $LEC + 1 ))
	    elif [ ${HGPULL} -eq 1 ] && [ ${HGUPDATE} -eq 0 ]; then
	      git fetch
	      [ $? -gt 0 ] && LEC=$(( $LEC + 1 ))
	    elif [ ${HGPULL} -eq 0 ] && [ ${HGUPDATE} -eq 1 ]; then
	      git update
        [ $? -gt 0 ] && LEC=$(( $LEC + 1 ))
	    fi # if pull and/or update
    elif [[ "$F" == *'/.hg' ]]; then
      # pull & update repository
      if [ ${HGPULL} -eq 1 ]; then 
        hg pull
        [ $? -gt 0 ] && LEC=$(( $LEC + 1 ))
      fi # if pull
      if [ ${HGUPDATE} -eq 1 ]; then
        hg update
        [ $? -gt 0 ] && LEC=$(( $LEC + 1 ))
      fi # if update    
    elif [[ "$F" == *'/.svn' ]]; then
      # update repository
      svn update
      [ $? -gt 0 ] && LEC=$(( $LEC + 1 ))
    else
      echo "Folder does not containa valid Git/Hg/Svn repository!"
      LEC=$(( $LEC + 1 ))
    fi # $F repo type  
    # evaluate results
    if [ ${LEC} -eq 0 ] 
      then OK=$(( ${OK} + 1 ))
      else ERR=$(( ${ERR} + 1 ))
    fi # if no error
    cd "${HGSROOT}" # back to root folder
    echo
} # UPDATE


ERR=0 # error counter
OK=0 # success counter

## update HGS Tools repository

if [ $HGSSYNC -gt 0 ]; then
  # feedback
  echo
  echo "   ***   Updating HGS Tools Repository   ***  "
  echo
  # find repo and update
  for R in .git .hg .svn; do   
      # check if this repo type exists here...
      if [ -d "$HGSTOOLS/$R" ]; then 
          UPDATE "$HGSTOOLS/$R" 
          break # there is only one...
      fi # check for repo
  done # $R
fi # $HGSSYNC
    

## update HGS model repositories

# set search expression based on recursion level
PATTERN='*/ */*/ */*/*/' # globbing expressions for search
HGSRCX='' 
GITSRCX='' 
SVNSRCX=''
for L in $( seq $RECLEV ); do
  # check which patterns actually apply
  DIR="${HGSROOT}/$( echo "${PATTERN}" | cut -d ' ' -f ${L} )/"
  ls -d ${DIR}/.hg &> /dev/null # check if any HG repositories present
  [ $? -eq 0 ] && [ ${HGSYNC} -eq 1 ] && HGSRCX="${HGSRCX} ${DIR}/.hg" 
  ls -d ${DIR}/.git &> /dev/null # check if any git repositories present
  [ $? -eq 0 ] && [ ${GITSYNC} -eq 1 ] && GITSRCX="${GITSRCX} ${DIR}/.git" 
  ls -d ${DIR}/.svn &> /dev/null # check if any SVN repositories present
  [ $? -eq 0 ] && [ ${SVNSYNC} -eq 1 ] && SVNSRCX="${SVNSRCX} ${DIR}/.svn" 
done

if [[ -n "${GITSRCX}" ]]; then
  # feedback
  echo
  echo "   ***   Updating HGS Model Repositories (GIT)   ***  "
  echo
  # update Git repositories (and pull)
  for RF in ${GITSRCX}; do 
      UPDATE "$RF"
  done
fi # if Git
  
if [[ -n "${HGSRCX}" ]]; then
  # feedback
  echo
  echo "   ***   Updating HGS Model Repositories (HG)   ***  "
  echo
  # update HG repositories (and pull)
  for RF in ${HGSRCX}; do 
      UPDATE "$RF"
  done
fi # if HG
  
if [[ -n "${SVNSRCX}" ]]; then
# feedback
echo
  echo "   ***   Updating HGS Model Repositories (SVN)   ***  "
echo
  # update SVN repositories (and pull)
  for RF in ${SVNSRCX}; do
      UPDATE "$RF"      
  done
fi # if SVN

echo
if [ $ERR == 0 ]
  then
    echo "   <<<   ALL ${OK} UPDATES OK   >>>   "
    echo
    exit 0
  else
    echo "   ###   WARNING: ${ERR} UPDATES FAILED OR INCOMPLETE!   ###   "
    echo "   >>>                 ${OK} UPDATES OK                <<<   "
    echo
    exit ${ERR}
fi
