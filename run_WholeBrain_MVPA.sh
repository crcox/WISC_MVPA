#!/bin/bash
# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes
# the specified command.
#
cleanup() {
  # Remove the Matlab runtime distribution
  if [ -f "r2013b.tar.gz" ]; then
    rm -v "r2013b.tar.gz"
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
  echo >&2 "Files at time of error/interrupt"
  echo >&2 "--------------------------------"
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
trap abort EXIT SIGTERM

set -e

## Download the runtime environment from SQUID
wget -q "http://proxy.chtc.wisc.edu/SQUID/r2013b.tar.gz"
#wget -q "http://proxy.chtc.wisc.edu/SQUID/crcox/sha1/r2013b.sha1"
#sha1sum -c "r2013b.sha1"
GET_MATLAB_STATUS=$?
i=0
while [ ! $GET_MATLAB_STATUS -eq 0 ] && [ $i -lt 5 ] ; do
  rm "r2013b.tar.gz"
  wget -q "http://proxy.chtc.wisc.edu/SQUID/r2013b.tar.gz"
  #shasum -c "r2013b.sha1"
  GET_MATLAB_STATUS=$?
  i++
done
tar xzf "r2013b.tar.gz"

## Download all large data files listed in URLS from SQUID
if [ ! -f URLS ]; then
  touch URLS
fi
if [ ! -f URLS_SHARED ]; then
  touch URLS_SHARED
fi
cat URLS URLS_SHARED > ALLURLS
cat ALLURLS
while read url; do
  wget -q "http://proxy.chtc.wisc.edu/SQUID/${url}"
done < ALLURLS

# Run the Matlab application
exe_name=$0
exe_dir=`dirname "$0"`
echo "------------------------------------------"
echo Setting up environment variables
MCRROOT="v82"
echo ---
LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
XAPPLRESDIR=${MCRROOT}/X11/app-defaults ;
export XAPPLRESDIR;
export LD_LIBRARY_PATH;
echo LD_LIBRARY_PATH is ${LD_LIBRARY_PATH};
eval "${exe_dir}/WholeBrain_MVPA"

# Exit successfully. Hooray!
trap success EXIT SIGTERM
