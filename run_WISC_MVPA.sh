#!/bin/bash
# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes
# the specified command.
#
download() {
  url=$1
  maxtries=$2
  try=0
  name=$(basename "$url")
  DOWNLOAD_STATUS=1
  while [ $DOWNLOAD_STATUS -gt 0 ]; do
    rm -f $name
    wget -q "${url}"
    wget -q "${url}.md5"
    md5sum -c "./${name}.md5"
    DOWNLOAD_STATUS=$?
    try=$((try+1))
    if [ $try -gt $maxtries ]; then echo "Download exceeded max tries. Exiting..."; exit; fi
  done
  rm "${name}.md5"
}
cleanup() {
  # Remove the Matlab runtime distribution
  if [ -f "r2015b.tar.gz" ]; then
    rm -v "r2015b.tar.gz"
  fi
  if [ -f "libXmu_libXt.el6.x86_64.tgz" ]; then
    rm -v "libXmu_libXt.el6.x86_64.tgz"
  fi
  # Check the home directory for any transfered files.
  if [ -f ALLURLS ]; then
    while read url; do
      fname=$(basename "$url")
      if [ -f "$fname" ]; then
        rm -v "$fname"
      fi
    done < ALLURLS
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

EXECUTABLE=$1
JOB_DIR=$2
PROXY_ROOT=$3
isOSG=$4

## Download all large data files listed in URLS from SQUID
# touch the files to ensure they exist
touch URLS
touch URLS_SHARED
cat URLS URLS_SHARED > ALLURLS
cat ALLURLS
while read url; do
  download "${PROXY_ROOT}/${url}" 5
done < ALLURLS

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
  echo Setting up environment variables
  ## Download the runtime environment from PROXY_ROOT
  download "${PROXY_ROOT}/crcox/r2015b.tar.gz" 5
  tar xzf "r2015b.tar.gz"

  # This is an attempt to fix broken environments by shipping libraries that are
  # missing on some nodes.
  download "${PROXY_ROOT}/crcox/libXmu_libXt.el6.x86_64.tgz" 5
  tar xzf "./libXmu_libXt.el6.x86_64.tgz"
  MCR_ROOT="`pwd`/v90"
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
