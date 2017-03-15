#!/bin/bash

# 
# freqresults.sh
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2017 Inorganic Materials Chemistry
# 
# 

RUNDIR=$(pwd)

if [ -z "$1" ] ; then
	REGEX="^.*/?freq/OUTCAR$"
	GETFREQ="NO"
	GETD2="NO"
elif [ "$1" == "FREQ" ] ; then
	GETFREQ="YES"
	GETD2="NO"
	PRINTHIGH="NO"
	if [ -z "$2" ] ; then
		REGEX="^.*/?freq/OUTCAR$"
	else
		REGEX="$2"
	fi
elif [ "$1" == "HIGH" ] ; then
	GETFREQ="NO"
	GETD2="NO"
	PRINTHIGH="YES"
	if [ -z "$2" ] ; then
		REGEX="^.*/?freq/OUTCAR$"
	else
		REGEX="$2"
	fi
elif [ "$1" == "D2" ] ; then
	GETFREQ="YES"
	GETD2="YES"
	PRINTHIGH="NO"
	if [ -z "$2" ] ; then
		REGEX="^.*/?freq/OUTCAR$"
	else
		REGEX="$2"
	fi
elif [ "$1" == "HD2" ] ; then
	GETFREQ="YES"
	GETD2="YES"
	PRINTHIGH="YES"
	if [ -z "$2" ] ; then
		REGEX="^.*/?freq/OUTCAR$"
	else
		REGEX="$2"
	fi
else
	GETFREQ="NO"
	GETD2="NO"
	PRINTHIGH="NO"
	REGEX="$1"
fi

DATE=$(date -u)
length=`expr length "$DATE"`
echo `printf "%0.s-" $(seq 2 $length)`
echo $DATE

while IFS= read -d $'\0' -r i;
do
	FREQ=$(grep "IBRION =      5" "$i")
	if [[ $FREQ == "" ]] ; then
		continue
	fi
	freqdir=$(dirname "${i}")
	optcar=$(dirname "${freqdir}")/OUTCAR
	echo
	length=`expr length $i`
	echo `printf "%0.s-" $(seq 1 $length)`
	echo "$i"
	echo `printf "%0.s-" $(seq 1 $length)`
	if [ ! -r "$optcar" ] ; then
		PRINTOPT="NO"
		echo "No valid OPTCAR found"
	else
		IBRION=$(grep "IBRION =      5" "$optcar")
		if [[ $IBRION != "" ]] ; then
			PRINTOPT="NO"
			echo "No valid OPTCAR found"
		else
			PRINTOPT="YES"
			echo "OPTCAR found: $optcar"
		fi
	fi
######  PRINTOPT  BEGIN  #######
	if [[ $PRINTOPT == "YES" ]] ; then
		vaspfull="$(head -n 1 $optcar | awk '{print $1;}')"
		vaspnum="$(head -n 1 $optcar | tr "." " " | awk '{print $2;}')"
		mode=0
		if [ $vaspnum == "4" ] ; then
		  mode=4
		fi
		if [ $vaspnum == "5" ] ; then
		  mode=5
		fi
		if [ $mode == 0 ] ; then
		  echo "WARNING: vasp version (" $vaspfull ") not recognized!"
		fi
		ENERGY=$(grep "y  w" "$optcar")
		if [[ $ENERGY != "" ]] ; then
			REACHED=$(grep "reached required accuracy - stopping structural energy minimisation" "$optcar")
			if [[ $REACHED == "" ]] ; then
				echo "WARNING: DID NOT REACH ENERGY MINIMISATION"
			fi
			((count=0))
			DRIFT=$(grep drift "$optcar")
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
			((count=0))
			for j in $ENERGY
			do
				((count=count+1))
				((mod=count%7))
				if [[ $mod == 0 ]] ; then
				TEMPENERGY=$j
				fi
			done
			EOPT=$TEMPENERGY
		fi
	fi
######  PRINTOPT  END  #######
	vaspfull="$(head -n 1 $i | awk '{print $1;}')"
	vaspnum="$(head -n 1 $i | tr "." " " | awk '{print $2;}')"
	mode=0
	if [ $vaspnum == "4" ] ; then
	  mode=4
	fi
	if [ $vaspnum == "5" ] ; then
	  mode=5
	fi
	if [ $mode == 0 ] ; then
	  echo "WARNING: vasp version (" $vaspfull ") not recognized!"
	fi
	ZPE="0.0"
	QVIB="0.0"
	ZPED="0.0"
	QVIBD="0.0"
	GFREQ=""
	HFREQ=$(grep "meV" $i | head -n 1 | awk '{print $8}')
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
			GFREQ="$(getfreq.py --least -e H -m 2.0 -o "$i")"
		else
			GFREQ="$(getfreq.py --least -e NOTHING -m 2.0 -o "$i")"
		fi
		ZPE=$(echo $GFREQ | awk '{print $1}')
		QVIB=$(echo $GFREQ | awk '{print $2}')
		ZPED=$(echo $GFREQ | awk '{print $3}')
		QVIBD=$(echo $GFREQ | awk '{print $4}')
		if [[ $(echo $GFREQ | awk '{print $3}') == "" ]] ; then
			ZPED=$ZPE
			QVIBD=$QVIB
		fi
	else
		ZPE=$(calcfreq.sh "$i" -l)
		QVIB=$(calcnu.sh "$i")
	fi
	echo
	if [[ $PRINTOPT == "YES" ]] ; then
		if [[ $GETD2 == "YES" ]] ; then
			if [[ $PRINTHIGH == "YES" ]] ; then
				echo -e "       Energy\t     ZPE\t     Qvib800\t   ZPE_D\t   Qvib800_D\t   HighFreq"
				echo -e $EOPT"\t"$ZPE"\t"$QVIB"\t"$ZPED"\t"$QVIBD"\t"$HFREQ
			else
				echo -e "       Energy\t     ZPE\t     Qvib800\t   ZPE_D\t   Qvib800_D"
				echo -e $EOPT"\t"$ZPE"\t"$QVIB"\t"$ZPED"\t"$QVIBD
			fi
		else
			if [[ $PRINTHIGH == "YES" ]] ; then
				echo -e "       Energy\t     ZPE\t     Qvib800\t   HighFreq"
				echo -e $EOPT"\t"$ZPE"\t"$QVIB"\t"$HFREQ
			else
				echo -e "       Energy\t     ZPE\t     Qvib800"
				echo -e $EOPT"\t"$ZPE"\t"$QVIB
			fi
		fi
	else
		if [[ $GETD2 == "YES" ]] ; then
			if [[ $PRINTHIGH == "YES" ]] ; then
				echo -e "     ZPE\t     Qvib800\t   ZPE_D\t   Qvib800_D\t   HighFreq"
				echo -e $ZPE"\t"$QVIB"\t"$ZPED"\t"$QVIBD"\t"$HFREQ
			else
				echo -e "     ZPE\t     Qvib800\t   ZPE_D\t   Qvib800_D"
				echo -e $ZPE"\t"$QVIB"\t"$ZPED"\t"$QVIBD
			fi
		else
			if [[ $PRINTHIGH == "YES" ]] ; then
				echo -e "     ZPE\t     Qvib800\t   HighFreq"
				echo -e $ZPE"\t"$QVIB"\t"$HFREQ
			else
				echo -e "     ZPE\t     Qvib800"
				echo -e $ZPE"\t"$QVIB
			fi
		fi
	fi
done < <(find . -regextype posix-extended -regex "$REGEX" -print0)
