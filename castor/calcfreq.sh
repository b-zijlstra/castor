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

awk -v mode="$mode" '

BEGIN{
}

/THz/{
  num++;
  if($10!="meV")
  {
  	tot+=$10;
        #print "Frequency #" num": "$10" meV";
  }
  }

END{
if (mode == 4)
{
  reportval=(tot/4)/1000;
}
if (mode == 5)
{
  reportval=(tot/2)/1000;
}
printf("Total contribution of real frequencies in eV: %10.6f", reportval);
print "";
}

' < "$1"
