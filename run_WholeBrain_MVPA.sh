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
SQUID="http://proxy.chtc.wisc.edu/SQUID/crcox"
## Download the runtime environment from SQUID
download "${SQUID}/r2013b.tar.gz" 5
tar xzf "r2013b.tar.gz"

# This is an attempt to fix broken environments by shipping libraries that are
# missing on some nodes.
download "${SQUID}/libXmu_libXt.el6.x86_64.tgz" 5
tar xzf "./libXmu_libXt.el6.x86_64.tgz"

## Download all large data files listed in URLS from SQUID
# touch the files to ensure they exist
touch URLS
touch URLS_SHARED
cat URLS URLS_SHARED > ALLURLS
cat ALLURLS
while read url; do
  download "http://proxy.chtc.wisc.edu/SQUID/${url}" 5
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
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:"./lib64"
XAPPLRESDIR=${MCRROOT}/X11/app-defaults ;
export XAPPLRESDIR;
export LD_LIBRARY_PATH;
echo LD_LIBRARY_PATH is ${LD_LIBRARY_PATH};
eval "${exe_dir}/WholeBrain_MVPA"

# Exit successfully. Hooray!
trap success EXIT SIGTERM
