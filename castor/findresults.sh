#!/bin/bash

# 
# findresults.sh
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2014 Inorganic Materials Chemistry
# 
# 

RUNDIR=$(pwd)

if [ -z "$1" ] ; then
	REGEX="^.*/?OUTCAR$"
	CHECKNEB="NO"
elif [ "$1" == "NEB" ] ; then
	CHECKNEB="YES"
	if [ -z "$2" ] ; then
		REGEX="^.*/?OUTCAR$"
	else
		REGEX="$2"
	fi
else
	CHECKNEB="NO"
	REGEX="$1"
fi

DATE=$(date -u)
length=`expr length "$DATE"`
echo `printf "%0.s-" $(seq 2 $length)`
echo $DATE

while IFS= read -d $'\0' -r i;
do
	ENERGY=$(grep "y  w" "$i")
	if [[ $ENERGY == "" ]] ; then
		continue
	fi

	if [ $CHECKNEB == "NO" ] ; then
		NEB=$(grep NEB: "$i")
		if [[ $NEB != "" ]] ; then
			continue
		fi
	fi

	FREQ=$(grep "using selective dynamics as specified on POSCAR" "$i")
	DRIFT=$(grep drift "$i")

	echo
	length=`expr length $i`
	echo `printf "%0.s-" $(seq 1 $length)`
	echo "$i"
	echo `printf "%0.s-" $(seq 1 $length)`

	if [[ $FREQ != "" ]] ; then
		echo "$(getfreq.py --less -e H -m 2.0 -o "$i")"
	else
		((count=0))
		for j in $ENERGY
		do
			((count=count+1))
			((mod=count%7))
			if [[ $mod == 0 ]] ; then
			TEMPENERGY=$j
			fi
		done
		((count=0))
		for j in $DRIFT
		do
			((count=count+1))
			((modx=(count-3)%5))
			((mody=(count-4)%5))
			((modz=(count-5)%5))
			if [[ $modx == 0 ]] ; then
			LASTDRIFTX=$j
			fi
			if [[ $mody == 0 ]] ; then
			LASTDRIFTY=$j
			fi
			if [[ $modz == 0 ]] ; then
			LASTDRIFTZ=$j
			fi
		done
		echo $TEMPENERGY
		REACHED=$(grep "reached required accuracy - stopping structural energy minimisation" "$i")
		if [[ $REACHED == "" ]] ; then
			echo "WARNING: DID NOT REACH ENERGY MINIMISATION"
		fi
		x1=$(echo $LASTDRIFTX 0.02 | awk '{print ($1 > $2) ? "true" : "false" }')
		x2=$(echo $LASTDRIFTX -0.02 | awk '{print ($1 < $2) ? "true" : "false" }')
		y1=$(echo $LASTDRIFTY 0.02 | awk '{print ($1 > $2) ? "true" : "false" }')
		y2=$(echo $LASTDRIFTY -0.02 | awk '{print ($1 < $2) ? "true" : "false" }')
		z1=$(echo $LASTDRIFTZ 0.02 | awk '{print ($1 > $2) ? "true" : "false" }')
		z2=$(echo $LASTDRIFTZ -0.02 | awk '{print ($1 < $2) ? "true" : "false" }')

		if [[ "$x1" == "true" || "$x2" == "true" ]] ; then
			echo "WARNING: DRIFT X = " $LASTDRIFTX
		fi
		if [[ "$y1" == "true" || "$y2" == "true" ]] ; then
			echo "WARNING: DRIFT Y = " $LASTDRIFTY
		fi
		if [[ "$z1" == "true" || "$z2" == "true" ]] ; then
			echo "WARNING: DRIFT Z = " $LASTDRIFTZ
		fi
	fi
done < <(find . -regextype posix-extended -regex $REGEX -print0)