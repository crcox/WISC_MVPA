#!/bin/bash

mkdir src && mv source_code.tar.gz src && cd src
tar xzvf source_code.tar.gz

make all
make clean

rm -rf bin/ src/

exit
