#include "BTagUtils.h"

//ctor
JetFlavourReader::JetFlavourReader(const std::string name_JetFlavourFile)
{

    name_JetFlavourFile_ = name_JetFlavourFile;

    FILE *f_flav = fopen(name_JetFlavourFile_.c_str(),"r");
    
    int eventId,lumiId,flav,njets,njets_f;
    float pt,eta,phi,energy;
     
    TLorentzVector p4;

    while(fscanf(f_flav,"%d %d %f %f %f %f %d %d %d \n", &lumiId, &eventId, &pt, &eta, &phi, &energy, &flav, &njets, &njets_f) !=EOF ){
            
            p4.SetPtEtaPhiE(pt,eta,phi,energy);

            AODjet_p4_[lumiId][eventId].push_back(p4);  
            AODjet_Flav_[lumiId][eventId].push_back(flav);          
      }

}
//---------------------------------------------------------------------------------------------------------------------------
//dtor
JetFlavourReader::~JetFlavourReader()
{
}
//---------------------------------------------------------------------------------------------------------------------------
int JetFlavourReader::getJetFlavour(const int& lumis, const int& event, const TLorentzVector* jetP4 )
{
    float jet_flavour = -1001.;

    for(unsigned int ii = 0; ii < AODjet_p4_[lumis][event].size(); ii++){
        dR_.push_back(jetP4->DeltaR(AODjet_p4_[lumis][event].at(ii)));
        dR_map_[jetP4->DeltaR(AODjet_p4_[lumis][event].at(ii))] = ii;
    }

    std::sort(dR_.begin(),dR_.end());

    if(dR_.at(0) < 0.3) jet_flavour = AODjet_Flav_[lumis][event].at(dR_map_[dR_.at(0)]);
    else jet_flavour = 0;
    
    dR_.clear();
    dR_map_.clear();

    return jet_flavour;
}
//---------------------------------------------------------------------------------------------------------------------------
//ctor
BtagSFReader::BtagSFReader(const std::string name_btagSFFile)
{

    name_btagSFFile_ = name_btagSFFile;

    TFile* btagSF_File_ = new TFile(name_btagSFFile_.c_str(),"READ");
        
    SFb_CSVL_ = (TF1*)btagSF_File_->Get("SFb_CSVL");
    h1_SFb_CSVL_ = (TH1F*)btagSF_File_->Get("h1_SFb_CSVL");

    SFb_CSVM_ = (TF1*)btagSF_File_->Get("SFb_CSVM");
    h1_SFb_CSVM_ = (TH1F*)btagSF_File_->Get("h1_SFb_CSVM");

    SFb_CSVT_ = (TF1*)btagSF_File_->Get("SFb_CSVT");
    h1_SFb_CSVT_ = (TH1F*)btagSF_File_->Get("h1_SFb_CSVT");

    SFudsg_CSVL_00_05_max_ = (TF1*)btagSF_File_->Get("SFudsg_CSVL_00_05_max");
    SFudsg_CSVL_00_05_mean_ = (TF1*)btagSF_File_->Get("SFudsg_CSVL_00_05_mean");
    SFudsg_CSVL_00_05_min_ = (TF1*)btagSF_File_->Get("SFudsg_CSVL_00_05_min");
    
    SFudsg_CSVL_05_10_max_ = (TF1*)btagSF_File_->Get("SFudsg_CSVL_05_10_max");
    SFudsg_CSVL_05_10_mean_ = (TF1*)btagSF_File_->Get("SFudsg_CSVL_05_10_mean");
    SFudsg_CSVL_05_10_min_ = (TF1*)btagSF_File_->Get("SFudsg_CSVL_05_10_min");

    SFudsg_CSVL_10_15_max_ = (TF1*)btagSF_File_->Get("SFudsg_CSVL_10_15_max");
    SFudsg_CSVL_10_15_mean_ = (TF1*)btagSF_File_->Get("SFudsg_CSVL_10_15_mean");
    SFudsg_CSVL_10_15_min_ = (TF1*)btagSF_File_->Get("SFudsg_CSVL_10_15_min");

    SFudsg_CSVL_15_24_max_ = (TF1*)btagSF_File_->Get("SFudsg_CSVL_15_24_max");
    SFudsg_CSVL_15_24_mean_ = (TF1*)btagSF_File_->Get("SFudsg_CSVL_15_24_mean");
    SFudsg_CSVL_15_24_min_ = (TF1*)btagSF_File_->Get("SFudsg_CSVL_15_24_min");

    SFudsg_CSVM_00_08_max_ = (TF1*)btagSF_File_->Get("SFudsg_CSVM_00_08_max");
    SFudsg_CSVM_00_08_mean_ = (TF1*)btagSF_File_->Get("SFudsg_CSVM_00_08_mean");
    SFudsg_CSVM_00_08_min_ = (TF1*)btagSF_File_->Get("SFudsg_CSVM_00_08_min");

    SFudsg_CSVM_08_16_max_ = (TF1*)btagSF_File_->Get("SFudsg_CSVM_08_16_max");
    SFudsg_CSVM_08_16_mean_ = (TF1*)btagSF_File_->Get("SFudsg_CSVM_08_16_mean");
    SFudsg_CSVM_08_16_min_ = (TF1*)btagSF_File_->Get("SFudsg_CSVM_08_16_min");

    SFudsg_CSVM_16_24_max_ = (TF1*)btagSF_File_->Get("SFudsg_CSVM_16_24_max");
    SFudsg_CSVM_16_24_mean_ = (TF1*)btagSF_File_->Get("SFudsg_CSVM_16_24_mean");
    SFudsg_CSVM_16_24_min_ = (TF1*)btagSF_File_->Get("SFudsg_CSVM_16_24_min");

    SFudsg_CSVT_00_24_max_ = (TF1*)btagSF_File_->Get("SFudsg_CSVT_00_24_max");
    SFudsg_CSVT_00_24_mean_ = (TF1*)btagSF_File_->Get("SFudsg_CSVT_00_24_mean");
    SFudsg_CSVT_00_24_min_ = (TF1*)btagSF_File_->Get("SFudsg_CSVT_00_24_min");
}
//---------------------------------------------------------------------------------------------------------------------------
//dtor
BtagSFReader::~BtagSFReader()
{
}
//---------------------------------------------------------------------------------------------------------------------------
float BtagSFReader::getSF(const TLorentzVector* jetP4,const float& flavour, const float& cvs_Btag)
{
      
    float SF = -1001.;

    if(fabs(flavour) == 5){
       if(cvs_Btag > 0.244){
          SF = SFb_CSVL_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVL_->GetXmax()) SF = SFb_CSVL_->Eval(SFb_CSVL_->GetXmax());
       }
       if(cvs_Btag > 0.679){
          SF = SFb_CSVM_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVM_->GetXmax()) SF = SFb_CSVM_->Eval(SFb_CSVM_->GetXmax());
       }
       if(cvs_Btag > 0.898){
          SF = SFb_CSVT_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVT_->GetXmax()) SF = SFb_CSVT_->Eval(SFb_CSVT_->GetXmax());
       }
    }
    if(fabs(flavour) == 4){
       if(cvs_Btag > 0.244){
          SF = SFb_CSVL_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVL_->GetXmax()) SF = SFb_CSVL_->Eval(SFb_CSVL_->GetXmax());
       }
       if(cvs_Btag > 0.679){
          SF = SFb_CSVM_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVM_->GetXmax()) SF = SFb_CSVM_->Eval(SFb_CSVM_->GetXmax());
       }
       if(cvs_Btag > 0.898){
          SF = SFb_CSVT_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVT_->GetXmax()) SF = SFb_CSVT_->Eval(SFb_CSVT_->GetXmax());
       }
    }
    if(fabs(flavour) != 0 && fabs(flavour) != 5 && fabs(flavour) != 4){
       if(cvs_Btag > 0.244){
          if(jetP4->Eta() > 0. && jetP4->Eta() < 0.5){
             SF = SFudsg_CSVL_00_05_mean_->Eval(jetP4->Pt());
             if(jetP4->Pt() > SFudsg_CSVL_00_05_mean_->GetXmax()) SF = SFudsg_CSVL_00_05_mean_->Eval(SFudsg_CSVL_00_05_mean_->GetXmax());
          }
          if(jetP4->Eta() >= 0.5 && jetP4->Eta() < 1.){
             SF = SFudsg_CSVL_05_10_mean_->Eval(jetP4->Pt());
             if(jetP4->Pt() > SFudsg_CSVL_05_10_mean_->GetXmax()) SF = SFudsg_CSVL_05_10_mean_->Eval(SFudsg_CSVL_05_10_mean_->GetXmax());
          }
          if(jetP4->Eta() >= 1. && jetP4->Eta() < 1.5){
             SF = SFudsg_CSVL_10_15_mean_->Eval(jetP4->Pt());
             if(jetP4->Pt() > SFudsg_CSVL_10_15_mean_->GetXmax()) SF = SFudsg_CSVL_10_15_mean_->Eval(SFudsg_CSVL_10_15_mean_->GetXmax());
          }
          if(jetP4->Eta() >= 1.5){
             SF = SFudsg_CSVL_15_24_mean_->Eval(jetP4->Pt());
             if(jetP4->Pt() > SFudsg_CSVL_15_24_mean_->GetXmax()) SF = SFudsg_CSVL_15_24_mean_->Eval(SFudsg_CSVL_15_24_mean_->GetXmax());
          }
       }
       if(cvs_Btag > 0.679){
          if(jetP4->Eta() > 0. && jetP4->Eta() < 0.8){
             SF = SFudsg_CSVM_00_08_mean_->Eval(jetP4->Pt());
             if(jetP4->Pt() > SFudsg_CSVM_00_08_mean_->GetXmax()) SF = SFudsg_CSVM_00_08_mean_->Eval(SFudsg_CSVM_00_08_mean_->GetXmax());
          }
          if(jetP4->Eta() >= 0.8 && jetP4->Eta() < 1.6){
             SF = SFudsg_CSVM_08_16_mean_->Eval(jetP4->Pt());
             if(jetP4->Pt() > SFudsg_CSVM_08_16_mean_->GetXmax()) SF = SFudsg_CSVM_08_16_mean_->Eval(SFudsg_CSVM_08_16_mean_->GetXmax());
          }
          if(jetP4->Eta() >= 1.6){
             SF = SFudsg_CSVM_16_24_mean_->Eval(jetP4->Pt());
             if(jetP4->Pt() > SFudsg_CSVM_16_24_mean_->GetXmax()) SF = SFudsg_CSVM_16_24_mean_->Eval(SFudsg_CSVM_16_24_mean_->GetXmax());
          }
       }
       if(cvs_Btag > 0.898){
          SF = SFudsg_CSVT_00_24_mean_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFudsg_CSVT_00_24_mean_->GetXmax()) SF = SFudsg_CSVT_00_24_mean_->Eval(SFudsg_CSVT_00_24_mean_->GetXmax());
       }
    }
    
    return SF;
}
//----------------------------------------------------------------------------------------------------------------------------
float BtagSFReader::getSFErrorUp(const TLorentzVector* jetP4,const float& flavour, const float& cvs_Btag)
{
      
    float SFerr = -1001.;
    float SFMax = -1001.;
    float SF = -1001.;

    if(fabs(flavour) == 5){
       if(cvs_Btag > 0.244){
          SF = SFb_CSVL_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVL_->GetXmax()) SF = SFb_CSVL_->Eval(SFb_CSVL_->GetXmax());
          SFerr = h1_SFb_CSVL_-> GetBinError(h1_SFb_CSVL_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVL_-> GetBinCenter(h1_SFb_CSVL_-> GetNbinsX())+0.5*h1_SFb_CSVL_-> GetBinWidth(h1_SFb_CSVL_-> GetNbinsX())){
             SFerr = 2*h1_SFb_CSVL_-> GetBinError(h1_SFb_CSVL_-> FindBin(h1_SFb_CSVL_-> GetBinCenter(h1_SFb_CSVL_-> GetNbinsX())+0.5*h1_SFb_CSVL_-> GetBinWidth(h1_SFb_CSVL_-> GetNbinsX())));
          }
       }
       if(cvs_Btag > 0.679){
          SF = SFb_CSVM_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVM_->GetXmax()) SF = SFb_CSVM_->Eval(SFb_CSVM_->GetXmax());
          SFerr = h1_SFb_CSVM_-> GetBinError(h1_SFb_CSVM_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVM_-> GetBinCenter(h1_SFb_CSVM_-> GetNbinsX())+0.5*h1_SFb_CSVM_-> GetBinWidth(h1_SFb_CSVM_-> GetNbinsX())){
             SFerr = 2*h1_SFb_CSVM_-> GetBinError(h1_SFb_CSVM_-> FindBin(h1_SFb_CSVM_-> GetBinCenter(h1_SFb_CSVM_-> GetNbinsX())+0.5*h1_SFb_CSVM_-> GetBinWidth(h1_SFb_CSVM_-> GetNbinsX())));
          }
       }
       if(cvs_Btag > 0.898){
          SF = SFb_CSVT_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVT_->GetXmax()) SF = SFb_CSVT_->Eval(SFb_CSVT_->GetXmax());
          SFerr = h1_SFb_CSVT_-> GetBinError(h1_SFb_CSVT_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVT_-> GetBinCenter(h1_SFb_CSVT_-> GetNbinsX())+0.5*h1_SFb_CSVT_-> GetBinWidth(h1_SFb_CSVT_-> GetNbinsX())){
             SFerr = 2*h1_SFb_CSVT_-> GetBinError(h1_SFb_CSVT_-> FindBin(h1_SFb_CSVT_-> GetBinCenter(h1_SFb_CSVT_-> GetNbinsX())+0.5*h1_SFb_CSVT_-> GetBinWidth(h1_SFb_CSVT_-> GetNbinsX())));
          }
       }
    }
    if(fabs(flavour) == 4){
       if(cvs_Btag > 0.244){
          SF = SFb_CSVL_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVL_->GetXmax()) SF = SFb_CSVL_->Eval(SFb_CSVL_->GetXmax());
          SFerr = 2*h1_SFb_CSVL_-> GetBinError(h1_SFb_CSVL_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVL_-> GetBinCenter(h1_SFb_CSVL_-> GetNbinsX())+0.5*h1_SFb_CSVL_-> GetBinWidth(h1_SFb_CSVL_-> GetNbinsX())){
             SFerr = 4*h1_SFb_CSVL_-> GetBinError(h1_SFb_CSVL_-> FindBin(h1_SFb_CSVL_-> GetBinCenter(h1_SFb_CSVL_-> GetNbinsX())+0.5*h1_SFb_CSVL_-> GetBinWidth(h1_SFb_CSVL_-> GetNbinsX())));
          }
       }
       if(cvs_Btag > 0.679){
          SF = SFb_CSVM_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVM_->GetXmax()) SF = SFb_CSVM_->Eval(SFb_CSVM_->GetXmax());
          SFerr = 2*h1_SFb_CSVM_-> GetBinError(h1_SFb_CSVM_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVM_-> GetBinCenter(h1_SFb_CSVM_-> GetNbinsX())+0.5*h1_SFb_CSVM_-> GetBinWidth(h1_SFb_CSVM_-> GetNbinsX())){
             SFerr = 4*h1_SFb_CSVM_-> GetBinError(h1_SFb_CSVM_-> FindBin(h1_SFb_CSVM_-> GetBinCenter(h1_SFb_CSVM_-> GetNbinsX())+0.5*h1_SFb_CSVM_-> GetBinWidth(h1_SFb_CSVM_-> GetNbinsX())));
          }
       }
       if(cvs_Btag > 0.898){
          SF = SFb_CSVT_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVT_->GetXmax()) SF = SFb_CSVT_->Eval(SFb_CSVT_->GetXmax());
          SFerr = 2*h1_SFb_CSVT_-> GetBinError(h1_SFb_CSVT_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVT_-> GetBinCenter(h1_SFb_CSVT_-> GetNbinsX())+0.5*h1_SFb_CSVT_-> GetBinWidth(h1_SFb_CSVT_-> GetNbinsX())){
             SFerr = 4*h1_SFb_CSVT_-> GetBinError(h1_SFb_CSVT_-> FindBin(h1_SFb_CSVT_-> GetBinCenter(h1_SFb_CSVT_-> GetNbinsX())+0.5*h1_SFb_CSVT_-> GetBinWidth(h1_SFb_CSVT_-> GetNbinsX())));
          }
       }
    }
    if(fabs(flavour) != 0 && fabs(flavour) != 5 && fabs(flavour) != 4){
       if(cvs_Btag > 0.244){
          if(jetP4->Eta() > 0. && jetP4->Eta() < 0.5){
             SF = SFudsg_CSVL_00_05_mean_->Eval(jetP4->Pt());
             SFMax = SFudsg_CSVL_00_05_max_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFMax);
             if(jetP4->Pt() > SFudsg_CSVL_00_05_mean_->GetXmax()){
                SF = SFudsg_CSVL_00_05_mean_->Eval(SFudsg_CSVL_00_05_mean_->GetXmax());
                SFMax = SFudsg_CSVL_00_05_max_->Eval(SFudsg_CSVL_00_05_max_->GetXmax());
                SFerr = 2*fabs(SF-SFMax);
             } 
          }
          if(jetP4->Eta() > 0.5 && jetP4->Eta() < 1.){
             SF = SFudsg_CSVL_05_10_mean_->Eval(jetP4->Pt());
             SFMax = SFudsg_CSVL_05_10_max_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFMax);
             if(jetP4->Pt() > SFudsg_CSVL_05_10_mean_->GetXmax()){
                SF = SFudsg_CSVL_05_10_mean_->Eval(SFudsg_CSVL_05_10_mean_->GetXmax());
                SFMax = SFudsg_CSVL_05_10_max_->Eval(SFudsg_CSVL_05_10_max_->GetXmax());
                SFerr = 2*fabs(SF-SFMax);
             } 
          }
          if(jetP4->Eta() > 1. && jetP4->Eta() < 1.5){
             SF = SFudsg_CSVL_10_15_mean_->Eval(jetP4->Pt());
             SFMax = SFudsg_CSVL_10_15_max_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFMax);
             if(jetP4->Pt() > SFudsg_CSVL_10_15_mean_->GetXmax()){
                SF = SFudsg_CSVL_10_15_mean_->Eval(SFudsg_CSVL_10_15_mean_->GetXmax());
                SFMax = SFudsg_CSVL_10_15_max_->Eval(SFudsg_CSVL_10_15_max_->GetXmax());
                SFerr = 2*fabs(SF-SFMax);
             } 
          }
          if(jetP4->Eta() > 1.5 && jetP4->Eta() < 2.4){
             SF = SFudsg_CSVL_15_24_mean_->Eval(jetP4->Pt());
             SFMax = SFudsg_CSVL_15_24_max_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFMax);
             if(jetP4->Pt() > SFudsg_CSVL_15_24_mean_->GetXmax()){
                SF = SFudsg_CSVL_15_24_mean_->Eval(SFudsg_CSVL_15_24_mean_->GetXmax());
                SFMax = SFudsg_CSVL_15_24_max_->Eval(SFudsg_CSVL_15_24_max_->GetXmax());
                SFerr = 2*fabs(SF-SFMax);
             } 
          }
          if(jetP4->Eta() >= 2.4){
             SF = SFudsg_CSVL_15_24_mean_->Eval(jetP4->Pt());
             SFMax = SFudsg_CSVL_15_24_max_->Eval(jetP4->Pt());
             SFerr = 2*fabs(SF-SFMax);
             if(jetP4->Pt() > SFudsg_CSVL_15_24_mean_->GetXmax()){
                SF = SFudsg_CSVL_15_24_mean_->Eval(SFudsg_CSVL_15_24_mean_->GetXmax());
                SFMax = SFudsg_CSVL_15_24_max_->Eval(SFudsg_CSVL_15_24_max_->GetXmax());
                SFerr = 4*fabs(SF-SFMax);
             } 
          }  
       }
       if(cvs_Btag > 0.679){
          if(jetP4->Eta() > 0. && jetP4->Eta() < 0.8){
             SF = SFudsg_CSVM_00_08_mean_->Eval(jetP4->Pt());
             SFMax = SFudsg_CSVM_00_08_max_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFMax);
             if(jetP4->Pt() > SFudsg_CSVM_00_08_mean_->GetXmax()){
                SF = SFudsg_CSVM_00_08_mean_->Eval(SFudsg_CSVM_00_08_mean_->GetXmax());
                SFMax = SFudsg_CSVM_00_08_max_->Eval(SFudsg_CSVM_00_08_max_->GetXmax());
                SFerr = 2*fabs(SF-SFMax);
             } 
          }
          if(jetP4->Eta() >= 0.8 && jetP4->Eta() < 1.6){
             SF = SFudsg_CSVM_08_16_mean_->Eval(jetP4->Pt());
             SFMax = SFudsg_CSVM_08_16_max_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFMax);
             if(jetP4->Pt() > SFudsg_CSVM_08_16_mean_->GetXmax()){
                SF = SFudsg_CSVM_08_16_mean_->Eval(SFudsg_CSVM_08_16_mean_->GetXmax());
                SFMax = SFudsg_CSVM_08_16_max_->Eval(SFudsg_CSVM_08_16_max_->GetXmax());
                SFerr = 2*fabs(SF-SFMax);
             } 
          }
          if(jetP4->Eta() >= 1.6 && jetP4->Eta() < 2.4){
             SF = SFudsg_CSVM_16_24_mean_->Eval(jetP4->Pt());
             SFMax = SFudsg_CSVM_16_24_max_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFMax);
             if(jetP4->Pt() > SFudsg_CSVM_16_24_mean_->GetXmax()){
                SF = SFudsg_CSVM_16_24_mean_->Eval(SFudsg_CSVM_16_24_mean_->GetXmax());
                SFMax = SFudsg_CSVM_16_24_max_->Eval(SFudsg_CSVM_16_24_max_->GetXmax());
                SFerr = 2*fabs(SF-SFMax);
             } 
          }
          if(jetP4->Eta() >= 2.4){
             SF = SFudsg_CSVM_16_24_mean_->Eval(jetP4->Pt());
             SFMax = SFudsg_CSVM_16_24_max_->Eval(jetP4->Pt());
             SFerr = 2*fabs(SF-SFMax);
             if(jetP4->Pt() > SFudsg_CSVM_16_24_mean_->GetXmax()){
                SF = SFudsg_CSVM_16_24_mean_->Eval(SFudsg_CSVM_16_24_mean_->GetXmax());
                SFMax = SFudsg_CSVM_16_24_max_->Eval(SFudsg_CSVM_16_24_max_->GetXmax());
                SFerr = 4*fabs(SF-SFMax);
             } 
          }
       }
       if(cvs_Btag > 0.898){
          if(jetP4->Eta() < 2.4){
             SF = SFudsg_CSVT_00_24_mean_->Eval(jetP4->Pt());
             SFMax = SFudsg_CSVT_00_24_max_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFMax);
             if(jetP4->Pt() > SFudsg_CSVT_00_24_mean_->GetXmax()){
                SF = SFudsg_CSVT_00_24_mean_->Eval(SFudsg_CSVT_00_24_mean_->GetXmax());
                SFMax = SFudsg_CSVT_00_24_max_->Eval(SFudsg_CSVT_00_24_max_->GetXmax());
                SFerr = 2*fabs(SF-SFMax);
             } 
          }
          if(jetP4->Eta() >= 2.4){
             SF = SFudsg_CSVT_00_24_mean_->Eval(jetP4->Pt());
             SFMax = SFudsg_CSVT_00_24_max_->Eval(jetP4->Pt());
             SFerr = 2*fabs(SF-SFMax);
             if(jetP4->Pt() > SFudsg_CSVT_00_24_mean_->GetXmax()){
                SF = SFudsg_CSVT_00_24_mean_->Eval(SFudsg_CSVT_00_24_mean_->GetXmax());
                SFMax = SFudsg_CSVT_00_24_max_->Eval(SFudsg_CSVT_00_24_max_->GetXmax());
                SFerr = 4*fabs(SF-SFMax);
             } 
          }
       }
    }

    return SFerr;
}
//----------------------------------------------------------------------------------------------------------------------------
float BtagSFReader::getSFErrorDown(const TLorentzVector* jetP4,const float& flavour, const float& cvs_Btag)
{
      
    float SFerr = -1001.;
    float SFmin = -1001.;
    float SF = -1001.;

    if(fabs(flavour) == 5){
       if(cvs_Btag > 0.244){
          SF = SFb_CSVL_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVL_->GetXmin()) SF = SFb_CSVL_->Eval(SFb_CSVL_->GetXmin());
          SFerr = h1_SFb_CSVL_-> GetBinError(h1_SFb_CSVL_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVL_-> GetBinCenter(h1_SFb_CSVL_-> GetNbinsX())+0.5*h1_SFb_CSVL_-> GetBinWidth(h1_SFb_CSVL_-> GetNbinsX())){
             SFerr = 2*h1_SFb_CSVL_-> GetBinError(h1_SFb_CSVL_-> FindBin(h1_SFb_CSVL_-> GetBinCenter(h1_SFb_CSVL_-> GetNbinsX())+0.5*h1_SFb_CSVL_-> GetBinWidth(h1_SFb_CSVL_-> GetNbinsX())));
          }
       }
       if(cvs_Btag > 0.679){
          SF = SFb_CSVM_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVM_->GetXmin()) SF = SFb_CSVM_->Eval(SFb_CSVM_->GetXmin());
          SFerr = h1_SFb_CSVM_-> GetBinError(h1_SFb_CSVM_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVM_-> GetBinCenter(h1_SFb_CSVM_-> GetNbinsX())+0.5*h1_SFb_CSVM_-> GetBinWidth(h1_SFb_CSVM_-> GetNbinsX())){
             SFerr = 2*h1_SFb_CSVM_-> GetBinError(h1_SFb_CSVM_-> FindBin(h1_SFb_CSVM_-> GetBinCenter(h1_SFb_CSVM_-> GetNbinsX())+0.5*h1_SFb_CSVM_-> GetBinWidth(h1_SFb_CSVM_-> GetNbinsX())));
          }
       }
       if(cvs_Btag > 0.898){
          SF = SFb_CSVT_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVT_->GetXmin()) SF = SFb_CSVT_->Eval(SFb_CSVT_->GetXmin());
          SFerr = h1_SFb_CSVT_-> GetBinError(h1_SFb_CSVT_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVT_-> GetBinCenter(h1_SFb_CSVT_-> GetNbinsX())+0.5*h1_SFb_CSVT_-> GetBinWidth(h1_SFb_CSVT_-> GetNbinsX())){
             SFerr = 2*h1_SFb_CSVT_-> GetBinError(h1_SFb_CSVT_-> FindBin(h1_SFb_CSVT_-> GetBinCenter(h1_SFb_CSVT_-> GetNbinsX())+0.5*h1_SFb_CSVT_-> GetBinWidth(h1_SFb_CSVT_-> GetNbinsX())));
          }
       }
    }
    if(fabs(flavour) == 4){
       if(cvs_Btag > 0.244){
          SF = SFb_CSVL_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVL_->GetXmin()) SF = SFb_CSVL_->Eval(SFb_CSVL_->GetXmin());
          SFerr = 2*h1_SFb_CSVL_-> GetBinError(h1_SFb_CSVL_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVL_-> GetBinCenter(h1_SFb_CSVL_-> GetNbinsX())+0.5*h1_SFb_CSVL_-> GetBinWidth(h1_SFb_CSVL_-> GetNbinsX())){
             SFerr = 4*h1_SFb_CSVL_-> GetBinError(h1_SFb_CSVL_-> FindBin(h1_SFb_CSVL_-> GetBinCenter(h1_SFb_CSVL_-> GetNbinsX())+0.5*h1_SFb_CSVL_-> GetBinWidth(h1_SFb_CSVL_-> GetNbinsX())));
          }
       }
       if(cvs_Btag > 0.679){
          SF = SFb_CSVM_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVM_->GetXmin()) SF = SFb_CSVM_->Eval(SFb_CSVM_->GetXmin());
          SFerr = 2*h1_SFb_CSVM_-> GetBinError(h1_SFb_CSVM_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVM_-> GetBinCenter(h1_SFb_CSVM_-> GetNbinsX())+0.5*h1_SFb_CSVM_-> GetBinWidth(h1_SFb_CSVM_-> GetNbinsX())){
             SFerr = 4*h1_SFb_CSVM_-> GetBinError(h1_SFb_CSVM_-> FindBin(h1_SFb_CSVM_-> GetBinCenter(h1_SFb_CSVM_-> GetNbinsX())+0.5*h1_SFb_CSVM_-> GetBinWidth(h1_SFb_CSVM_-> GetNbinsX())));
          }
       }
       if(cvs_Btag > 0.898){
          SF = SFb_CSVT_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVT_->GetXmin()) SF = SFb_CSVT_->Eval(SFb_CSVT_->GetXmin());
          SFerr = 2*h1_SFb_CSVT_-> GetBinError(h1_SFb_CSVT_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVT_-> GetBinCenter(h1_SFb_CSVT_-> GetNbinsX())+0.5*h1_SFb_CSVT_-> GetBinWidth(h1_SFb_CSVT_-> GetNbinsX())){
             SFerr = 4*h1_SFb_CSVT_-> GetBinError(h1_SFb_CSVT_-> FindBin(h1_SFb_CSVT_-> GetBinCenter(h1_SFb_CSVT_-> GetNbinsX())+0.5*h1_SFb_CSVT_-> GetBinWidth(h1_SFb_CSVT_-> GetNbinsX())));
          }
       }
    }
    if(fabs(flavour) != 0 && fabs(flavour) != 5 && fabs(flavour) != 4){
       if(cvs_Btag > 0.244){
          if(jetP4->Eta() > 0. && jetP4->Eta() < 0.5){
             SF = SFudsg_CSVL_00_05_mean_->Eval(jetP4->Pt());
             SFmin = SFudsg_CSVL_00_05_min_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFmin);
             if(jetP4->Pt() > SFudsg_CSVL_00_05_mean_->GetXmin()){
                SF = SFudsg_CSVL_00_05_mean_->Eval(SFudsg_CSVL_00_05_mean_->GetXmin());
                SFmin = SFudsg_CSVL_00_05_min_->Eval(SFudsg_CSVL_00_05_min_->GetXmin());
                SFerr = 2*fabs(SF-SFmin);
             } 
          }
          if(jetP4->Eta() > 0.5 && jetP4->Eta() < 1.){
             SF = SFudsg_CSVL_05_10_mean_->Eval(jetP4->Pt());
             SFmin = SFudsg_CSVL_05_10_min_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFmin);
             if(jetP4->Pt() > SFudsg_CSVL_05_10_mean_->GetXmin()){
                SF = SFudsg_CSVL_05_10_mean_->Eval(SFudsg_CSVL_05_10_mean_->GetXmin());
                SFmin = SFudsg_CSVL_05_10_min_->Eval(SFudsg_CSVL_05_10_min_->GetXmin());
                SFerr = 2*fabs(SF-SFmin);
             } 
          }
          if(jetP4->Eta() > 1. && jetP4->Eta() < 1.5){
             SF = SFudsg_CSVL_10_15_mean_->Eval(jetP4->Pt());
             SFmin = SFudsg_CSVL_10_15_min_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFmin);
             if(jetP4->Pt() > SFudsg_CSVL_10_15_mean_->GetXmin()){
                SF = SFudsg_CSVL_10_15_mean_->Eval(SFudsg_CSVL_10_15_mean_->GetXmin());
                SFmin = SFudsg_CSVL_10_15_min_->Eval(SFudsg_CSVL_10_15_min_->GetXmin());
                SFerr = 2*fabs(SF-SFmin);
             } 
          }
          if(jetP4->Eta() > 1.5 && jetP4->Eta() < 2.4){
             SF = SFudsg_CSVL_15_24_mean_->Eval(jetP4->Pt());
             SFmin = SFudsg_CSVL_15_24_min_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFmin);
             if(jetP4->Pt() > SFudsg_CSVL_15_24_mean_->GetXmin()){
                SF = SFudsg_CSVL_15_24_mean_->Eval(SFudsg_CSVL_15_24_mean_->GetXmin());
                SFmin = SFudsg_CSVL_15_24_min_->Eval(SFudsg_CSVL_15_24_min_->GetXmin());
                SFerr = 2*fabs(SF-SFmin);
             } 
          }
          if(jetP4->Eta() >= 2.4){
             SF = SFudsg_CSVL_15_24_mean_->Eval(jetP4->Pt());
             SFmin = SFudsg_CSVL_15_24_min_->Eval(jetP4->Pt());
             SFerr = 2*fabs(SF-SFmin);
             if(jetP4->Pt() > SFudsg_CSVL_15_24_mean_->GetXmin()){
                SF = SFudsg_CSVL_15_24_mean_->Eval(SFudsg_CSVL_15_24_mean_->GetXmin());
                SFmin = SFudsg_CSVL_15_24_min_->Eval(SFudsg_CSVL_15_24_min_->GetXmin());
                SFerr = 4*fabs(SF-SFmin);
             } 
          }  
       }
       if(cvs_Btag > 0.679){
          if(jetP4->Eta() > 0. && jetP4->Eta() < 0.8){
             SF = SFudsg_CSVM_00_08_mean_->Eval(jetP4->Pt());
             SFmin = SFudsg_CSVM_00_08_min_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFmin);
             if(jetP4->Pt() > SFudsg_CSVM_00_08_mean_->GetXmin()){
                SF = SFudsg_CSVM_00_08_mean_->Eval(SFudsg_CSVM_00_08_mean_->GetXmin());
                SFmin = SFudsg_CSVM_00_08_min_->Eval(SFudsg_CSVM_00_08_min_->GetXmin());
                SFerr = 2*fabs(SF-SFmin);
             } 
          }
          if(jetP4->Eta() >= 0.8 && jetP4->Eta() < 1.6){
             SF = SFudsg_CSVM_08_16_mean_->Eval(jetP4->Pt());
             SFmin = SFudsg_CSVM_08_16_min_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFmin);
             if(jetP4->Pt() > SFudsg_CSVM_08_16_mean_->GetXmin()){
                SF = SFudsg_CSVM_08_16_mean_->Eval(SFudsg_CSVM_08_16_mean_->GetXmin());
                SFmin = SFudsg_CSVM_08_16_min_->Eval(SFudsg_CSVM_08_16_min_->GetXmin());
                SFerr = 2*fabs(SF-SFmin);
             } 
          }
          if(jetP4->Eta() >= 1.6 && jetP4->Eta() < 2.4){
             SF = SFudsg_CSVM_16_24_mean_->Eval(jetP4->Pt());
             SFmin = SFudsg_CSVM_16_24_min_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFmin);
             if(jetP4->Pt() > SFudsg_CSVM_16_24_mean_->GetXmin()){
                SF = SFudsg_CSVM_16_24_mean_->Eval(SFudsg_CSVM_16_24_mean_->GetXmin());
                SFmin = SFudsg_CSVM_16_24_min_->Eval(SFudsg_CSVM_16_24_min_->GetXmin());
                SFerr = 2*fabs(SF-SFmin);
             } 
          }
          if(jetP4->Eta() >= 2.4){
             SF = SFudsg_CSVM_16_24_mean_->Eval(jetP4->Pt());
             SFmin = SFudsg_CSVM_16_24_min_->Eval(jetP4->Pt());
             SFerr = 2*fabs(SF-SFmin);
             if(jetP4->Pt() > SFudsg_CSVM_16_24_mean_->GetXmin()){
                SF = SFudsg_CSVM_16_24_mean_->Eval(SFudsg_CSVM_16_24_mean_->GetXmin());
                SFmin = SFudsg_CSVM_16_24_min_->Eval(SFudsg_CSVM_16_24_min_->GetXmin());
                SFerr = 4*fabs(SF-SFmin);
             } 
          }
       }
       if(cvs_Btag > 0.898){
          if(jetP4->Eta() < 2.4){
             SF = SFudsg_CSVT_00_24_mean_->Eval(jetP4->Pt());
             SFmin = SFudsg_CSVT_00_24_min_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFmin);
             if(jetP4->Pt() > SFudsg_CSVT_00_24_mean_->GetXmin()){
                SF = SFudsg_CSVT_00_24_mean_->Eval(SFudsg_CSVT_00_24_mean_->GetXmin());
                SFmin = SFudsg_CSVT_00_24_min_->Eval(SFudsg_CSVT_00_24_min_->GetXmin());
                SFerr = 2*fabs(SF-SFmin);
             } 
          }
          if(jetP4->Eta() >= 2.4){
             SF = SFudsg_CSVT_00_24_mean_->Eval(jetP4->Pt());
             SFmin = SFudsg_CSVT_00_24_min_->Eval(jetP4->Pt());
             SFerr = 2*fabs(SF-SFmin);
             if(jetP4->Pt() > SFudsg_CSVT_00_24_mean_->GetXmin()){
                SF = SFudsg_CSVT_00_24_mean_->Eval(SFudsg_CSVT_00_24_mean_->GetXmin());
                SFmin = SFudsg_CSVT_00_24_min_->Eval(SFudsg_CSVT_00_24_min_->GetXmin());
                SFerr = 4*fabs(SF-SFmin);
             } 
          }
       }
    }

    return SFerr;
}
//---------------------------------------------------------------------------------------------------------------------------
//ctor
BtagEfficiencyReader::BtagEfficiencyReader(const std::string name_btagEfficienciesFile)
{

    name_btagEfficienciesFile_ = name_btagEfficienciesFile;
    
    TFile* btagEfficiency_File_ = new TFile(name_btagEfficienciesFile_.c_str(),"READ");

    h2_BTaggingEff_b_L_ = (TH2F*)btagEfficiency_File_->Get("h2_BTaggingEff_b_L");
    h2_BTaggingEff_b_M_ = (TH2F*)btagEfficiency_File_->Get("h2_BTaggingEff_b_M");
    h2_BTaggingEff_b_T_ = (TH2F*)btagEfficiency_File_->Get("h2_BTaggingEff_b_T");
    
    h2_BTaggingEff_c_L_ = (TH2F*)btagEfficiency_File_->Get("h2_BTaggingEff_c_L");
    h2_BTaggingEff_c_M_ = (TH2F*)btagEfficiency_File_->Get("h2_BTaggingEff_c_M");
    h2_BTaggingEff_c_T_ = (TH2F*)btagEfficiency_File_->Get("h2_BTaggingEff_c_T");

    h2_BTaggingEff_udsg_L_ = (TH2F*)btagEfficiency_File_->Get("h2_BTaggingEff_udsg_L");
    h2_BTaggingEff_udsg_M_ = (TH2F*)btagEfficiency_File_->Get("h2_BTaggingEff_udsg_M");
    h2_BTaggingEff_udsg_T_ = (TH2F*)btagEfficiency_File_->Get("h2_BTaggingEff_udsg_T");

    
}
//---------------------------------------------------------------------------------------------------------------------------
//dtor
BtagEfficiencyReader::~BtagEfficiencyReader()
{
}
//---------------------------------------------------------------------------------------------------------------------------
float BtagEfficiencyReader::getBtagEfficiency(const TLorentzVector* jetP4, const float& csv_Btag, const int& jet_flavour)
{

     float eff = -1001.;

     if(fabs(jet_flavour) == 5){
        if(csv_Btag > 0.244) eff = h2_BTaggingEff_b_L_->GetBinContent(h2_BTaggingEff_b_L_->FindBin(jetP4->Pt(),jetP4->Eta()));
        if(csv_Btag > 0.679) eff = h2_BTaggingEff_b_M_->GetBinContent(h2_BTaggingEff_b_M_->FindBin(jetP4->Pt(),jetP4->Eta()));
        if(csv_Btag > 0.898) eff = h2_BTaggingEff_b_T_->GetBinContent(h2_BTaggingEff_b_T_->FindBin(jetP4->Pt(),jetP4->Eta()));
     }
     
     if(fabs(jet_flavour) == 4){
        if(csv_Btag > 0.244) eff = h2_BTaggingEff_c_L_->GetBinContent(h2_BTaggingEff_c_L_->FindBin(jetP4->Pt(),jetP4->Eta()));
        if(csv_Btag > 0.679) eff = h2_BTaggingEff_c_M_->GetBinContent(h2_BTaggingEff_c_M_->FindBin(jetP4->Pt(),jetP4->Eta()));
        if(csv_Btag > 0.898) eff = h2_BTaggingEff_c_T_->GetBinContent(h2_BTaggingEff_c_T_->FindBin(jetP4->Pt(),jetP4->Eta()));
     }
     
     if(fabs(jet_flavour) != 0 && fabs(jet_flavour) != 4 && fabs(jet_flavour) != 5){
        if(csv_Btag > 0.244) eff = h2_BTaggingEff_udsg_L_->GetBinContent(h2_BTaggingEff_udsg_L_->FindBin(jetP4->Pt(),jetP4->Eta()));
        if(csv_Btag > 0.679) eff = h2_BTaggingEff_udsg_M_->GetBinContent(h2_BTaggingEff_udsg_M_->FindBin(jetP4->Pt(),jetP4->Eta()));
        if(csv_Btag > 0.898) eff = h2_BTaggingEff_udsg_T_->GetBinContent(h2_BTaggingEff_udsg_T_->FindBin(jetP4->Pt(),jetP4->Eta()));
     }

     return eff;
}
//---------------------------------------------------------------------------------------------------------------------------
float BtagEfficiencyReader::getBtagEfficiencyError(const TLorentzVector* jetP4, const float& csv_Btag, const int& jet_flavour){

     float eff_err = -1001.;

     if(fabs(jet_flavour) == 5){
        if(csv_Btag > 0.244) eff_err = h2_BTaggingEff_b_L_->GetBinError(h2_BTaggingEff_b_L_->FindBin(jetP4->Pt(),jetP4->Eta()));
        if(csv_Btag > 0.679) eff_err = h2_BTaggingEff_b_M_->GetBinError(h2_BTaggingEff_b_M_->FindBin(jetP4->Pt(),jetP4->Eta()));
        if(csv_Btag > 0.898) eff_err = h2_BTaggingEff_b_T_->GetBinError(h2_BTaggingEff_b_T_->FindBin(jetP4->Pt(),jetP4->Eta()));
     }
     
     if(fabs(jet_flavour) == 4){
        if(csv_Btag > 0.244) eff_err = h2_BTaggingEff_c_L_->GetBinError(h2_BTaggingEff_c_L_->FindBin(jetP4->Pt(),jetP4->Eta()));
        if(csv_Btag > 0.679) eff_err = h2_BTaggingEff_c_M_->GetBinError(h2_BTaggingEff_c_M_->FindBin(jetP4->Pt(),jetP4->Eta()));
        if(csv_Btag > 0.898) eff_err = h2_BTaggingEff_c_T_->GetBinError(h2_BTaggingEff_c_T_->FindBin(jetP4->Pt(),jetP4->Eta()));
     }
     
     if(fabs(jet_flavour) != 0 && fabs(jet_flavour) != 4 && fabs(jet_flavour) != 5){
        if(csv_Btag > 0.244) eff_err = h2_BTaggingEff_udsg_L_->GetBinError(h2_BTaggingEff_udsg_L_->FindBin(jetP4->Pt(),jetP4->Eta()));
        if(csv_Btag > 0.679) eff_err = h2_BTaggingEff_udsg_M_->GetBinError(h2_BTaggingEff_udsg_M_->FindBin(jetP4->Pt(),jetP4->Eta()));
        if(csv_Btag > 0.898) eff_err = h2_BTaggingEff_udsg_T_->GetBinError(h2_BTaggingEff_udsg_T_->FindBin(jetP4->Pt(),jetP4->Eta()));
     }

     return eff_err;
}
//---------------------------------------------------------------------------------------------------------------------------
