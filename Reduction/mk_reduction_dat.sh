#!/bin/bash

. setup.sh

#rm mc_radion_signal/*.dat
#./AnalysisScripts/mk_reduction_dat.py /store/cmst3/user/obondu/H2GGLOBE/Radion/processed/V15_00_05 ${storedir} mc_radion_signal.txt

#rm mc_radion_SMhiggs_summer12_rd1/*.dat
#./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V15_00_08/mc/Summer12_RD1 ${storedir} mc_radion_SMhiggs_summer12_rd1.txt

#rm data_2012/*.dat
#./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V15_00_05/data ${storedir}/data data_2012.txt

rm data_2012_v2/*.dat
./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V15_00_05/data ${storedir}/data data_2012_v2.txt

#rm mc_radion_bkg_diphoton/*.dat
#./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V15_00_05/mc/Summer12_DR53X-PU_S10_START53_V7C ${storedir} mc_radion_bkg_diphoton.txt


#rm mc_dy_summer12_s10/*.dat
#./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/Resonant_HH/processed/V14_00_02/mc/Summer12_S10_8TeV ${storedir}/mc/Summer12_S10_8TeV mc_dy_summer12_s10.txt

#rm mc_radion_bkg_summer12_s10/*.dat
#./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/Resonant_HH/processed/V14_00_03/mc/Summer12_S10_8TeV ${storedir}/mc/Summer12_S10_8TeV mc_radion_bkg_summer12_s10.txt
#rm mc_radion_bkg_summer12_s10_v2/*.dat
#./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/Resonant_HH/processed/V14_00_04/mc/Summer12_S10_8TeV ${storedir}/mc/Summer12_S10_8TeV mc_radion_bkg_summer12_s10_v2.txt


tar zcf ${version}.tgz  AnalysisScripts/{common,reduction,aux,*.py}
