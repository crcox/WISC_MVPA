#!/bin/bash

URL=$(cat DATAURL)
TARFILE=$(basename ${URL})
if [ ! -f $TARFILE ]; then
  wget ${URL}
fi

if [ ! -d data ]; then
  mkdir -v data
fi
tar -C data/ -xzvf FacePlaceObject.tgz
