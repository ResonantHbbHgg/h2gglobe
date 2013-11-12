#!/bin/bash

#####
##### radion_reduction_v10
#####

### SIGNAL RADION
#./submit_reduction.sh mc_radion_signal RadionToHHTo2G2B\* 5
#./submit_reduction.sh mc_radion_signal GravitonToHH_2Gamma_2b\* 5

### DATA
#./submit_reduction.sh data2012_RERECO \*Run2012B\* 125
#./submit_reduction.sh data2012_RERECO \*Run2012D\* 250
#./submit_reduction.sh data2012_RERECO \*Run2012C\* 250
#./submit_reduction.sh data2012_RERECO \*Run2012A\* 100
## RESUB

### BACKGROUNDS
./submit_reduction.sh mc_radion_SM-HH GluGluToHHTo2B2G\* 10
#./submit_reduction.sh mc_Summer12_RD1 \*HToGG_M-124\* 5
#./submit_reduction.sh mc_Summer12_RD1 \*HToGG_M-125\* 5
#./submit_reduction.sh mc_Summer12_RD1 \*HToGG_M-126\* 5
#./submit_reduction.sh mc_Summer12_RD1 DiPhoton\*sherpa\* 200
#./submit_reduction.sh mc_Summer12_RD1 GJet\*pythia6\* 100
#./submit_reduction.sh mc_Summer12_RD1 QCD_Pt-30\* 20
#./submit_reduction.sh mc_Summer12_RD1 QCD_Pt-40\* 20
#./submit_reduction.sh mc_Summer12_RD1 DY\* 200
#./submit_reduction.sh mc_Summer12_RD1_proxy GJet\*sherpa\* 200
#./submit_reduction.sh mc_Summer12_RD1_proxy Diphoton\*EW4\* 50
## RESUB

wait
