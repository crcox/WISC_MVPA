#!/bin/bash -e

INFILE=$1
tmp=${INFILE%.*}
OUTFILE=${tmp##*/}
mnimaster="~/data/FacePlaceObject/CommonBrains/TT_N27_funcres.nii"

3dUndump -master "$mnimaster" -xyz -datum float -prefix $OUTFILE $INFILE 2>tmp
