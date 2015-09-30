#!/bin/bash

#find source_code/ -type f ! -regex ".*/\..*" ! -regex ".*/ProjectTemplate/.*" ! -regex ".*/bin/.*" ! -regex ".*/log/.*" ! -regex ".*/iterativelasso/runIterativeLasso.m | xargs cp -t build/

mkdir src && mv source_code.tar.gz src && cd src
tar xzvf source_code.tar.gz

cd ../
mkdir build/
mv ./src/src/* build/
mv Makefile build/

cd build/
find ../src/dependencies/jsonlab/ -type f -exec cp '{}' . \;
