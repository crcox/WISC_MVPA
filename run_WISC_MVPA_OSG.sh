#!/bin/bash

set -e
set -x

EXECUTABLE=$1
# Run the Matlab application
# NB: This script needs to run to the end, even if there are errors.
set +e
source /cvmfs/oasis.opensciencegrid.org/osg/modules/lmod/current/init/bash
set -e
  module load matlab/2015b

chmod +x ${EXECUTABLE}
eval "./${EXECUTABLE}"

