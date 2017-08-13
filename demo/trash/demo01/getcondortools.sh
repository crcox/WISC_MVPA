#!/bin/bash

if [ ! -d ${HOME}/src ]; then
  mkdir -v ${HOME}/src
fi
git clone https://github.com/crcox/condortools.git ${HOME}/src/condortools

demodir=${pwd}
cd -v ${HOME}/src/condortools
python setup.py install --user --prefix=
cd -v ${demodir}
