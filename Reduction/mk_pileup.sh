#!/bin/bash

set -x

source "version.sh"

mkdir -p pileup

#dir=/store/group/phys_higgs/Resonant_HH/reduced/${version}
dir=/store/cmst3/user/obondu/H2GGLOBE/Radion/reduced/${version}/mc

cd pileup 
rm *.root
cmsLs $dir | awk '{ print $5}' | grep '/' | grep -v broken | grep -v root > samples.txt 

echo "Hadding and copying"
../parallel --eta --joblog parallel_pileupMerger.log --progress "python ../Macros/pileup/pileupMerger.py --putBack {1}" :::: samples.txt


