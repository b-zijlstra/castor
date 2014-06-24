#!/bin/bash

#METAL SETTINGS
metal=Co
lc=2.4933
hc=1.6152
packing=HCP

# metal=Ru
# lc=2.7266
# hc=1.5771
# packing=HCP

#metal=Rh
#lc=3.8034
#hc=1.0000
#packing=FCC

#SLAB SETTINGS (not needed for bulk)
miller_1=0
miller_2=0
miller_3=0
dim_a=3
dim_b=3
layers=5
vacuum=15


#RUN PYTHON SCRIPT FOR BULK
./structuregen.py $metal $lc $hc $packing

#RUN PYTHON SCRIPT FOR SURFACE
#./structuregen.py $metal $lc $hc $packing $miller_1 $miller_2 $miller_3 $dim_a $dim_b $layers $vacuum
