#!/bin/bash

. setup.sh

### RADION SIGNAL
rm mc_radion_signal_v2/*.dat
./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/Summer12_DR53X-PU_RD1_START53_V7N ${storedir}/mc mc_radion_signal_v2.txt

### DATA
rm data2012_RERECO_v2/*.dat
./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V15_00_05/data ${storedir}/data data2012_RERECO_v2.txt

### SM HIGGS AND BACKGROUNDS
rm mc_Summer12_RD1_v2/*.dat
./AnalysisScripts/mk_reduction_dat.py - ${storedir}/mc/Summer12_RD1 mc_Summer12_RD1_v2.txt
rm mc_Summer12_RD1_proxy_v2/*.dat
./AnalysisScripts/mk_reduction_dat.py - ${storedir}/mc/Summer12_RD1 mc_Summer12_RD1_proxy_v2.txt



#### OTHER OLD STUFF
#rm mc_radion_SMhiggs_summer12_rd1/*.dat
#./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/Summer12_RD1 ${storedir} mc_radion_SMhiggs_summer12_rd1.txt

#rm data_2012/*.dat
#./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V15_00_05/data ${storedir}/data data_2012.txt

#rm mc_radion_bkg_diphoton/*.dat
#./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V15_00_05/mc/Summer12_DR53X-PU_S10_START53_V7C ${storedir} mc_radion_bkg_diphoton.txt


#rm mc_dy_summer12_s10/*.dat
#./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/Resonant_HH/processed/V14_00_02/mc/Summer12_S10_8TeV ${storedir}/mc/Summer12_S10_8TeV mc_dy_summer12_s10.txt

#rm mc_radion_bkg_summer12_s10/*.dat
#./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/Resonant_HH/processed/V14_00_03/mc/Summer12_S10_8TeV ${storedir}/mc/Summer12_S10_8TeV mc_radion_bkg_summer12_s10.txt
#rm mc_radion_bkg_summer12_s10_v2/*.dat
#./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/Resonant_HH/processed/V14_00_04/mc/Summer12_S10_8TeV ${storedir}/mc/Summer12_S10_8TeV mc_radion_bkg_summer12_s10_v2.txt


#rm mc_spin2_summer12_s10/*.dat
#./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V15_00_02/mc/Summer12_S10_8TeV ${storedir}/mc/Summer12_S10_8TeV mc_spin2_summer12_s10.txt

#DoubleElectron_Run2012A-22Jan2013-v1_AOD
#DoubleElectron_Run2012C-22Jan2013-v1_AOD
#DoublePhoton_Run2012B-22Jan2013-v1_AOD
#DoublePhoton_Run2012D-22Jan2013-v1_v3
#DoubleElectron_Run2012B-22Jan2013-v1_AOD
#DoubleElectron_Run2012D-22Jan2013-v1_v2
#DoublePhoton_Run2012C-22Jan2013-v2_AOD
#Photon_Run2012A_22Jan2013-v1_AOD

### rm data2012_RERECO/*.dat
### #./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V15_00_05/data /store/group/phys_higgs/cmshgg/reduced/rereco_june2013/data data2012_RERECO.txt
### ./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V15_00_05/data ${storedir}/data data2012_RERECO.txt

## rm data2011_RERECO/*.dat
## ./AnalysisScripts/mk_reduction_dat.py - ${storedir}/data data2011_RERECO.txt

##rm data2012_RERECO/*.dat
##./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V15_00_05/data /store/group/phys_higgs/cmshgg/reduced/rereco_june2013/data data2012_RERECO.txt
##
##
##rm mc2012RD_v2_1/*.dat
##./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V15_00_08/mc ${storedir}/mc mc2012RD_v2_1.txt
##
##rm mc2012RD_v2_2/*.dat
##./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/Summer12_RD1 ${storedir}/mc mc2012RD_v2_2.txt
##
##rm mc2012RD_v2_3/*.dat
##./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/Summer12_DR53X-PU_RD1_START53_V7N ${storedir}/mc mc2012RD_v2_3.txt





wd=$PWD
cd AnalysisScripts
tar cf $wd/${version}.tar *.py $(find common reduction baseline massfac_mva_binned full_mva_binned jetanalysis photonjet -name \*.dat -or -name \*.py) aux common python 
cd -

tar rf ${version}.tar JSON *.sh
gzip -f ${version}.tar

git tag -a ${version} -m "Tag used for reduction ${group}/${version}"

## tar zcf ${version}.tgz  AnalysisScripts/{common,reduction,aux,*.py}
