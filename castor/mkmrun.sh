#!/bin/bash

# 
# mkmrun.sh
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2014 Inorganic Materials Chemistry
# 
# 

RUNDIR=$(pwd)

if [ -z "$1" ] ; then
	REGEX="^.*/?input.mkm$"
else
	REGEX="$1"
fi

while IFS= read -d $'\0' -r i;
do
	cd "$(dirname "$i")"
	mkmcxx -i "$i"
	REGEX2="^.*.plt$"
	while IFS= read -d $'\0' -r j;
	do
		cd "$(dirname "$j")"
		gnuplot "$(basename "$j")"
		cd "$(dirname "$i")"
	done < <(find "`pwd`" -regextype posix-extended -regex $REGEX2 -print0)
	cd "$RUNDIR"
done < <(find "`pwd`" -regextype posix-extended -regex $REGEX -print0)