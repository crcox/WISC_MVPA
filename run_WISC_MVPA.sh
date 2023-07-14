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
  
  for x in "$@"; do
    rm -fv $(basename "$x")
  done
  
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

  cleanup "$@"

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

  cleanup "$@"

  echo "An error occured. Exiting ..." >&2
  exit 1
}
success() {
  echo '
*************
** SUCCESS **
*************
'
  cleanup "$@"

  exit 0
}

# If an exit or interrupt occurs while the script is executing, run the abort
# function.
trap 'abort "$@"' EXIT
trap 'terminated  "$@"' SIGTERM SIGKILL

set -e
set -x

EXECUTABLE=$1
STAGING_PATH=$2
JOBID=$3
shift 3
for x in "$@"; do
  cp "$x" .
done

# Run the Matlab application
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

chmod +x ${EXECUTABLE}
eval "./${EXECUTABLE}"

# Copy results to STAGING_PATH
mkdir ${STAGING_PATH}/${JOBID}
mv -v results.mat ${STAGING_PATH}/${JOBID}/results.mat

# Exit successfully. Hooray!
trap 'success "$@"' EXIT SIGTERM
