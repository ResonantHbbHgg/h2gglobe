#!/bin/bash

#
### DATA
#

## ./submit_reduction.sh    data_2012             Photon_Run2012A-13Jul2012-v1                    10
## ./submit_reduction.sh    data_2012             Photon_Run2012A-recover-06Aug2012-v1            1
## ./submit_reduction.sh    data_2012             DoublePhoton_Run2012B-13Jul2012-v1              80
## ./submit_reduction.sh    data_2012             DoublePhoton_Run2012C-24Aug2012-v2              10
## ./submit_reduction.sh    data_2012             DoublePhoton_Run2012C-EcalRecover_11Dec2012-v1  4
## ./submit_reduction.sh    data_2012             DoublePhoton_Run2012C-PromptReco-v2             150
## ./submit_reduction.sh    data_2012             DoublePhoton_Run2012C-PromptReco-v2_sub2        7
#
##./submit_reduction.sh    data_2012             DoublePhoton_Run2012D-PromptReco-v1             120
## ./submit_reduction.sh    data_2012             DoublePhoton_Run2012D-PromptReco-v1_sub2        10
## ./submit_reduction.sh    data_2012             run2012D_DE                                     80

#
### SIGNAL RADION
#

## bash submit_reduction.sh mc_radion_signal RadionToHHTo2G2B_M-1000_TuneZ2star_8TeV-madgraph
## bash submit_reduction.sh mc_radion_signal RadionToHHTo2G2B_M-1500_TuneZ2star_8TeV-madgraph
## bash submit_reduction.sh mc_radion_signal RadionToHHTo2G2B_M-300_TuneZ2star_8TeV-madgraph
## bash submit_reduction.sh mc_radion_signal RadionToHHTo2G2B_M-500_TuneZ2star_8TeV-madgraph
## bash submit_reduction.sh mc_radion_signal RadionToHHTo2G2B_M-700_TuneZ2star_8TeV-madgraph

bash submit_reduction.sh mc_radion_signal_v2 GravitonToHH_2Gamma_2b_M-1000_TuneZ2star_8TeV-Madgraph_pythia6
bash submit_reduction.sh mc_radion_signal_v2 GravitonToHH_2Gamma_2b_M-1500_TuneZ2star_8TeV-Madgraph_pythia6
bash submit_reduction.sh mc_radion_signal_v2 GravitonToHH_2Gamma_2b_M-300_TuneZ2star_8TeV-Madgraph_pythia6
bash submit_reduction.sh mc_radion_signal_v2 GravitonToHH_2Gamma_2b_M-500_TuneZ2star_8TeV-Madgraph_pythia6
bash submit_reduction.sh mc_radion_signal_v2 GravitonToHH_2Gamma_2b_M-700_TuneZ2star_8TeV-Madgraph_pythia6
bash submit_reduction.sh mc_radion_signal_v2 RadionToHHTo2G2B_M-1000_TuneZ2star_8TeV-nm-madgraph
bash submit_reduction.sh mc_radion_signal_v2 RadionToHHTo2G2B_M-1500_TuneZ2star_8TeV-nm-madgraph
bash submit_reduction.sh mc_radion_signal_v2 RadionToHHTo2G2B_M-300_TuneZ2star_8TeV-nm-madgraph
bash submit_reduction.sh mc_radion_signal_v2 RadionToHHTo2G2B_M-500_TuneZ2star_8TeV-nm-madgraph
bash submit_reduction.sh mc_radion_signal_v2 RadionToHHTo2G2B_M-700_TuneZ2star_8TeV-nm-madgraph

#
### MC HIGGS
#

## bash submit_reduction.sh mc_radion_SMhiggs_summer12_s10 GluGluToHToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1
## bash submit_reduction.sh mc_radion_SMhiggs_summer12_s10 TTH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1
## bash submit_reduction.sh mc_radion_SMhiggs_summer12_s10 VBF_HToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1
## bash submit_reduction.sh mc_radion_SMhiggs_summer12_s10 WH_ZH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1


#
### MC BACKGROUND
#

### bash submit_reduction.sh mc_radion_bkg_summer12_s10 DiPhotonBox_Pt-10To25_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1 10
## bash submit_reduction.sh mc_radion_bkg_summer12_s10 DiPhotonBox_Pt-250ToInf_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1 10
## bash submit_reduction.sh mc_radion_bkg_summer12_s10 DiPhotonBox_Pt-25To250_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1 10
##bash submit_reduction.sh mc_radion_bkg_summer12_s10_v2 DiPhotonJets_8TeV-madgraph-tarball-v2_Summer12_DR53X-PU_S10_START53_V7A-v1 20
##bash submit_reduction.sh mc_radion_bkg_summer12_s10 DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1 50
##bash submit_reduction.sh mc_radion_bkg_summer12_s10 DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_ZMuMuFilter 50
### bash submit_reduction.sh mc_radion_bkg_summer12_s10 GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_ff 50
## bash submit_reduction.sh mc_radion_bkg_summer12_s10 GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_pf 50
## bash submit_reduction.sh mc_radion_bkg_summer12_s10 GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_pp 50
### bash submit_reduction.sh mc_radion_bkg_summer12_s10 GJet_Pt40_v2_ff 50
## bash submit_reduction.sh mc_radion_bkg_summer12_s10 GJet_Pt40_v2_pf 50
## bash submit_reduction.sh mc_radion_bkg_summer12_s10 GJet_Pt40_v2_pp 50
## bash submit_reduction.sh mc_radion_bkg_summer12_s10 QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_ff 50
## bash submit_reduction.sh mc_radion_bkg_summer12_s10 QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_pf 50
### bash submit_reduction.sh mc_radion_bkg_summer12_s10 QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_pp 50
##bash submit_reduction.sh mc_radion_bkg_summer12_s10 QCD_Pt-40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_ff 1
##bash submit_reduction.sh mc_radion_bkg_summer12_s10 QCD_Pt-40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_pf 1
### bash submit_reduction.sh mc_radion_bkg_summer12_s10 QCD_Pt-40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_pp 1


############# EXTRA
#./submit_reduction.sh     mc_bkg_summer12_s10_2 DiPhotonJets_8TeV-madgraph-tarball-v2_Summer12_DR53X-PU_S10_START53_V7A-v1              40
#./submit_reduction.sh     mc_dy_summer12_s10    DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1   50
#
#./submit_reduction.sh     mc_bkg_summer12_s10   DiPhotonBox_Pt-10To25_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1     5
#./submit_reduction.sh     mc_bkg_summer12_s10   DiPhotonBox_Pt-250ToInf_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1   10
#./submit_reduction.sh     mc_bkg_summer12_s10   DiPhotonBox_Pt-25To250_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1    10
#./submit_reduction.sh     mc_bkg_summer12_s10   GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_pp  20
#./submit_reduction.sh     mc_bkg_summer12_s10   GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_pf  30
#./submit_reduction.sh     mc_bkg_summer12_s10   GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_ff  20
#./submit_reduction.sh     mc_bkg_summer12_s10   QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_pp   10
#./submit_reduction.sh     mc_bkg_summer12_s10   QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_pf   20
#./submit_reduction.sh     mc_bkg_summer12_s10   QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_ff   30
#./submit_reduction.sh     mc_bkg_summer12_s10   QCD_Pt-40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_pp       10
#./submit_reduction.sh     mc_bkg_summer12_s10   QCD_Pt-40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_pf       20
#./submit_reduction.sh     mc_bkg_summer12_s10   QCD_Pt-40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_ff       30
#./submit_reduction.sh     mc_bkg_summer12_s10   WGToLNuG_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1   10
#./submit_reduction.sh     mc_bkg_summer12_s10   WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1    10
#./submit_reduction.sh     mc_bkg_summer12_s10   WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1    10
#./submit_reduction.sh     mc_bkg_summer12_s10   WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1   10
#./submit_reduction.sh     mc_bkg_summer12_s10   ZGToLLG_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1    10
#./submit_reduction.sh     mc_bkg_summer12_s10   ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v3    10
#./submit_reduction.sh     mc_bkg_summer12_s10   ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1   10
#./submit_reduction.sh     mc_bkg_summer12_s10   ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1   10
#
#./submit_reduction.sh     mc_bkg_summer12_s10   Wpgg_dR02   8
#./submit_reduction.sh     mc_bkg_summer12_s10   Wmgg_dR02   8
#./submit_reduction.sh     mc_bkg_summer12_s10   Zgg_dR02    8 
#./submit_reduction.sh     mc_bkg_summer12_s10   ttgg_dR02   2
#
#./submit_reduction.sh     mc_bkg_summer12_s10   GJet_Pt40_v2_ff   40
#./submit_reduction.sh     mc_bkg_summer12_s10   GJet_Pt40_v2_pf   80
#./submit_reduction.sh     mc_bkg_summer12_s10   GJet_Pt40_v2_pp   60
#
#./submit_reduction.sh     mc_bkg_summer12_s10   TTJets_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S8_START52_V9-v1    20
#
#./submit_reduction.sh     mc_bkg_summer12_s10   TTJets_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1    20
#./submit_reduction.sh     mc_bkg_summer12_s10   TTJets_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S8_START52_V9-v1    20
#
#
#./submit_reduction.sh     mc_sig_summer12_s10   \*     5
#./submit_reduction.sh     mc_spin2_summer12_s10 \*     5
#


wait
