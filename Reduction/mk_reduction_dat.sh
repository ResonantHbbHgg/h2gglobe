#!/bin/bash

. setup.sh

#rm mc_radion_signal/*.dat
#./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/Resonant_HH/processed/V14_00_08/mc ${storedir} mc_radion_signal.txt
#rm mc_radion_signal_v2/*.dat
#./AnalysisScripts/mk_reduction_dat.py /store/cmst3/user/obondu/H2GGLOBE/Radion/processed/V14_00_08 ${storedir} mc_radion_signal_v2.txt
#rm mc_radion_signal_v3/*.dat
#./AnalysisScripts/mk_reduction_dat.py /store/cmst3/user/obondu/H2GGLOBE/Radion/processed/V14_00_08 ${storedir} mc_radion_signal_v3.txt

#rm mc_radion_SMhiggs_summer12_s10/*.dat
#./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/Resonant_HH/processed/V14_00_03/mc/Summer12_S10_8TeV ${storedir}/mc/Summer12_S10_8TeV mc_radion_SMhiggs_summer12_s10.txt
#rm data_2012/*.dat
#./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/Resonant_HH/processed/V14_00_06/data ${storedir}/data data_2012.txt




#rm mc_dy_summer12_s10/*.dat
#./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/Resonant_HH/processed/V14_00_02/mc/Summer12_S10_8TeV ${storedir}/mc/Summer12_S10_8TeV mc_dy_summer12_s10.txt

#rm mc_radion_bkg_summer12_s10/*.dat
#./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/Resonant_HH/processed/V14_00_03/mc/Summer12_S10_8TeV ${storedir}/mc/Summer12_S10_8TeV mc_radion_bkg_summer12_s10.txt
rm mc_radion_bkg_summer12_s10_v2/*.dat
./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/Resonant_HH/processed/V14_00_04/mc/Summer12_S10_8TeV ${storedir}/mc/Summer12_S10_8TeV mc_radion_bkg_summer12_s10_v2.txt


tar zcf ${version}.tgz  AnalysisScripts/{common,reduction,aux,*.py}
