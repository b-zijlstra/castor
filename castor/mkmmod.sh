#!/bin/bash

# 
# mkmmod.sh
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2014 Inorganic Materials Chemistry
# 
# 

RUNDIR=$(pwd)

H2RATIO=$1
COSITE=$2
H2SITE=$3
COSTICK=$4
H2STICK=$5
CODES=$6
H2DES=$7
PRESSURE=$8
if [ -z "$9" ] ; then
    INPUT="source.mkm"
    DIRNAME="H2ratio"$H2RATIO"_COsite"$COSITE"_H2site"$H2SITE"_COstick"$COSTICK"_H2stick"$H2STICK"_COdes"$CODES"_H2des"$H2DES"_pressure"$PRESSURE
else
    INPUT="$9"
    if [ -z "${10}" ] ; then
    	DIRNAME="H2ratio"$H2RATIO"_COsite"$COSITE"_H2site"$H2SITE"_COstick"$COSTICK"_H2stick"$H2STICK"_COdes"$CODES"_H2des"$H2DES"_pressure"$PRESSURE
    else
        DIRNAME="${10}"
    fi
fi

while [ -d "$DIRNAME" ];
do
	DIRNAME=$DIRNAME"new"
done

mkdir $DIRNAME
cp $INPUT "$DIRNAME/input.mkm"
sed -i 's/@H2ratio@/'$H2RATIO'/g' "$DIRNAME/input.mkm"
sed -i 's/@COsite@/'$COSITE'/g' "$DIRNAME/input.mkm"
sed -i 's/@H2site@/'$H2SITE'/g' "$DIRNAME/input.mkm"
sed -i 's/@COstick@/'$COSTICK'/g' "$DIRNAME/input.mkm"
sed -i 's/@H2stick@/'$H2STICK'/g' "$DIRNAME/input.mkm"
sed -i 's/@COdes@/'$CODES'/g' "$DIRNAME/input.mkm"
sed -i 's/@H2des@/'$H2DES'/g' "$DIRNAME/input.mkm"
sed -i 's/@pressure@/'$PRESSURE'/g' "$DIRNAME/input.mkm"
