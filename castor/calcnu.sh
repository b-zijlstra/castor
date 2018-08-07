#!/bin/bash
set -e
if [ -z "$1" ] ; then
  echo 'OUTCAR file name required.'
  exit 1
fi
if [ ! -r "$1" ] ; then
  echo "$1 does not exists or is not readable."
  exit 1
fi

if [ -z "$2" ] ; then
	TEMPERATURE="800"
else
	TEMPERATURE="$2"
fi

vaspfull="$(head -n 1 $1 | awk '{print $1;}')"
vaspnum="$(head -n 1 $1 | tr "." " " | awk '{print $2;}')"
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

# count lines and divide by 2
LINES=`grep "meV" $1 | wc -l`
if [ $mode == 4 ] ; then
  LINES=$(($LINES/2))
fi

# grab all frequency strengths in meV
LIST=`grep "meV" $1 | tail -n $LINES | awk '{print $10}'`

calcnu.py $TEMPERATURE $LIST