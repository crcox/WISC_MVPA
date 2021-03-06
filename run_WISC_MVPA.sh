#!/bin/bash
# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes
# the specified command.
#
cleanup() {
  # Remove the Matlab runtime distribution
  if [ -f "r2018b.tar.gz" ]; then
    rm -v "r2018b.tar.gz"
  fi
  echo "all clean"
}
abort() {
  echo >&2 '
*************
** ABORTED **
*************
'
  echo >&2 "Files at time of error"
  echo >&2 "----------------------"
  ls >&2 -l

  cleanup

  echo "An error occured. Exiting ..." >&2
  exit 1
}
terminated() {
  echo >&2 '
****************
** TERMINATED **
****************
'
  echo >&2 "Files at time of interrupt"
  echo >&2 "--------------------------"
  ls >&2 -l

  cleanup

  echo "An error occured. Exiting ..." >&2
  exit 1
}
success() {
  echo '
*************
** SUCCESS **
*************
'
  cleanup

  exit 0
}

# If an exit or interrupt occurs while the script is executing, run the abort
# function.
trap abort EXIT
trap terminated SIGTERM SIGKILL

set -e
set -x

EXECUTABLE=WISC_MVPA
JOB_DIR=$2
PROXY_ROOT=$3
isOSG=$4

## Download all large data files listed in URLS from SQUID
# touch the files to ensure they exist

# Run the Matlab application
if [ $isOSG = "True" ]; then
  # This script needs to run to the end, even if there are errors.
  set +e
  source /cvmfs/oasis.opensciencegrid.org/osg/modules/lmod/current/init/bash
  set -e
  module load matlab/2015b

else
  # CHTC
  echo "------------------------------------------"
  echo "Setting up environment variables"
  tar xzf "r2018b.tar.gz"

  # This is an attempt to fix broken environments by shipping libraries that are
  # missing on some nodes.
  MCR_ROOT="`pwd`/v95"
  mkdir cache && export MCR_CACHE_ROOT="`pwd`/cache"

  echo "MCR_ROOT: ${MCR_ROOT}"
  echo "MCR_CACHE_ROOT: ${MCR_CACHE_ROOT}"

  echo ---
  LD_LIBRARY_PATH=.:${MCR_ROOT}/runtime/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCR_ROOT}/bin/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCR_ROOT}/sys/os/glnxa64;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:"./lib64"
  XAPPLRESDIR=${MCR_ROOT}/X11/app-defaults ;
  export XAPPLRESDIR;
  export LD_LIBRARY_PATH;
  echo LD_LIBRARY_PATH is ${LD_LIBRARY_PATH};

fi

chmod +x ${EXECUTABLE}
eval "./${EXECUTABLE}"

# Exit successfully. Hooray!
trap success EXIT SIGTERM
