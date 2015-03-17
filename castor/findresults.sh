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
	GETFREQ="NO"
	GETD2="NO"
elif [ "$1" == "NEB" ] ; then
	CHECKNEB="YES"
	GETFREQ="NO"
	GETD2="NO"
	if [ -z "$2" ] ; then
		REGEX="^.*/?OUTCAR$"
	else
		REGEX="$2"
	fi
elif [ "$1" == "FREQ" ] ; then
	CHECKNEB="NO"
	GETFREQ="YES"
	GETD2="NO"
	if [ -z "$2" ] ; then
		REGEX="^.*/?OUTCAR$"
	else
		REGEX="$2"
	fi
elif [ "$1" == "D2" ] ; then
	CHECKNEB="NO"
	GETFREQ="YES"
	GETD2="YES"
	if [ -z "$2" ] ; then
		REGEX="^.*/?OUTCAR$"
	else
		REGEX="$2"
	fi

else
	CHECKNEB="NO"
	GETFREQ="NO"
	GETD2="NO"
	REGEX="$1"
fi

DATE=$(date -u)
length=`expr length "$DATE"`
echo `printf "%0.s-" $(seq 2 $length)`
echo $DATE

while IFS= read -d $'\0' -r i;
do
	vasp="$(head -n 1 $i | awk '{print $1;}')"
	mode=0
	if [ $vasp == "vasp.4.6.31" ] ; then
	  mode=4
	fi
	if [ $vasp == "vasp.4.6.38" ] ; then
	  mode=4
	fi
	if [ $vasp == "vasp.5.2" ] ; then
	  mode=5
	fi
	if [ $vasp == "vasp.5.3.3" ] ; then
	  mode=5
	fi
	if [ $vasp == "vasp.5.3.5" ] ; then
	  mode=5
	fi
	if [ $mode == 0 ] ; then
	  echo "WARNING: vasp version (" $vasp ") not recognized!"
	fi

	ENERGY=$(grep "y  w" "$i")
	if [[ $ENERGY == "" ]] ; then
		continue
	fi

	NEB=$(grep NEB: "$i")
	if [[ $CHECKNEB == "NO" && $NEB != "" ]] ; then
		continue
	elif [[ $CHECKNEB == "YES" && $NEB == "" ]] ; then
		continue
	fi

	FREQ=$(grep "using selective dynamics as specified on POSCAR" "$i")
	DRIFT=$(grep drift "$i")

	echo
	length=`expr length $i`
	echo `printf "%0.s-" $(seq 1 $length)`
	echo "$i"
	echo `printf "%0.s-" $(seq 1 $length)`

	if [[ $FREQ != "" ]] ; then
		FREQCOUNT=$(grep meV "$i" | wc -l)
		if [[ $mode == 4 ]] ; then
			IMAGS=$(awk -v a=$FREQCOUNT 'BEGIN{} /THz/{ num++; if($10=="meV") { if(num<=a/2) print num " f/i= "$9" meV"; } } END{}' < "$i")
		fi
		if [[ $mode == 5 ]] ; then
			IMAGS=$(awk 'BEGIN{} /THz/{ num++; if($10=="meV") { print num " f/i= "$9" meV"; } } END{}' < "$i")
		fi
		if [[ $IMAGS != "" ]] ; then
			IMAGCOUNT=$(echo "$IMAGS" | wc -l)
			if [[ $IMAGCOUNT == "1" ]] ; then
				echo "1 imaginary frequency found:"
			else
				echo $IMAGCOUNT" imaginary frequencies found:"
			fi
			echo "$IMAGS"
		else
			echo "No imaginary frequencies found"
		fi

		if [[ $GETFREQ == "YES" ]] ; then
			if [[ $GETD2 == "YES" ]] ; then
				echo "$(getfreq.py --less -e H -m 2.0 -o "$i")"
			else
			echo "$(getfreq.py --less -e NOTHING -m 2.0 -o "$i")"
			fi
		else
			calcfreq "$i"
		fi

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