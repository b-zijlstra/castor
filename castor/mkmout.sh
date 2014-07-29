#!/bin/bash

# 
# mkmout.sh
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2014 Inorganic Materials Chemistry
# 
# 

RUNDIR=$(pwd)

FIRSTDIR="true"
while IFS= read -d $'\0' -r i;
do
	if [ $FIRSTDIR == "true" ] ; then
		paste "$i" | awk '{print $1}' > dydt.txt
		FIRSTDIR="false"
		DIRSTRING="Directory-> "
	fi
	DIRSTRING=$DIRSTRING" "$i" "$i
	paste "$i" | awk '{print $4}' > temptext.txt
	paste {dydt.txt,temptext.txt} > temptext2.txt
	mv temptext2.txt dydt.txt
done < <(find . -regextype posix-extended -regex "^.*/?dydt.dat$" -print0)
echo $DIRSTRING | sed "s/ /\t/g" > temptext.txt
cat temptext.txt dydt.txt > temptext2.txt
mv temptext2.txt dydt.txt
rm temptext.txt

FIRSTDIR="true"
while IFS= read -d $'\0' -r i;
do
	if [ $FIRSTDIR == "true" ] ; then
		paste "$i" | awk '{print $1}' > cvt.txt
		FIRSTDIR="false"
		DIRSTRING="Directory-> "
	fi
	DIRSTRING=$DIRSTRING" "$i" "$i
	paste "$i" | awk '{print $6, $9}' | sed "s/ /\t/g" > temptext.txt
	paste {cvt.txt,temptext.txt} > temptext2.txt
	mv temptext2.txt cvt.txt
done < <(find . -regextype posix-extended -regex "^.*/?cvt.dat$" -print0)
echo $DIRSTRING | sed "s/ /\t/g" > temptext.txt
cat temptext.txt cvt.txt > temptext2.txt
mv temptext2.txt cvt.txt
rm temptext.txt

FIRSTDIR="true"
while IFS= read -d $'\0' -r i;
do
	if [ $FIRSTDIR == "true" ] ; then
		paste "$i" | awk '{print $1}' > orders.txt
		FIRSTDIR="false"
		DIRSTRING="Directory-> "
	fi
	DIRSTRING=$DIRSTRING" "$i" "$i
	paste "$i" | awk '{print $2, $3}' | sed "s/ /\t/g" > temptext.txt
	paste {orders.txt,temptext.txt} > temptext2.txt
	mv temptext2.txt orders.txt
done < <(find . -regextype posix-extended -regex "^.*/?orders.dat$" -print0)
echo $DIRSTRING | sed "s/ /\t/g" > temptext.txt
cat temptext.txt orders.txt > temptext2.txt
mv temptext2.txt orders.txt
rm temptext.txt

FIRSTDIR="true"
while IFS= read -d $'\0' -r i;
do
	if [ $FIRSTDIR == "true" ] ; then
		paste "$i" | awk '{print $1}' > drc.txt
		FIRSTDIR="false"
		DIRSTRING="Directory-> "
	fi
	DIRSTRING=$DIRSTRING" "$i" "$i" "$i" "$i" "$i" "$i
	paste "$i" | sed "s/ //g" | awk '{print $11, $12, $13, $14, $15, $16}' | sed "s/ /\t/g" > temptext.txt
	paste {drc.txt,temptext.txt} > temptext2.txt
	mv temptext2.txt drc.txt
done < <(find . -regextype posix-extended -regex "^.*/?[^t]drc.dat$" -print0)
echo $DIRSTRING | sed "s/ /\t/g" > temptext.txt
cat temptext.txt drc.txt > temptext2.txt
mv temptext2.txt drc.txt
rm temptext.txt