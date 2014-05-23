#include "BTagUtils.h"
#define BUDEBUG 0

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
    if(BUDEBUG) std::cerr << "Entering JetFlavourReader::getJetFlavour" << std::endl;
    int jet_flavour = 0;

    for(unsigned int ii = 0; ii < AODjet_p4_[lumis][event].size(); ii++){
        dR_.push_back(jetP4->DeltaR(AODjet_p4_[lumis][event].at(ii)));
        dR_map_[jetP4->DeltaR(AODjet_p4_[lumis][event].at(ii))] = ii;
    }

    std::sort(dR_.begin(),dR_.end());

    if(dR_.at(0) < 0.3) jet_flavour = AODjet_Flav_[lumis][event].at(dR_map_[dR_.at(0)]);
    else jet_flavour = 0;
    
    dR_.clear();
    dR_map_.clear();

    if(BUDEBUG) std::cerr << "Leaving JetFlavourReader::getJetFlavour" << std::endl;
    return jet_flavour;
}
//---------------------------------------------------------------------------------------------------------------------------
//ctor
BtagSFReader::BtagSFReader(const std::string name_btagSFFile)
{

    name_btagSFFile_ = name_btagSFFile;

    TFile* btagSF_File_ = TFile::Open(name_btagSFFile_.c_str());
        
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
float BtagSFReader::getSF(const TLorentzVector* jetP4,const float& flavour, std::string WP)
{
    if(BUDEBUG) std::cerr << "Entering BtagSFReader::getSF" << std::endl;
      
    float SF = -1001.;
    
    if(fabs(flavour) == 5){
       if(WP == "loose"){
          SF = SFb_CSVL_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVL_->GetXmax()) SF = SFb_CSVL_->Eval(SFb_CSVL_->GetXmax());
          if(jetP4->Pt() < SFb_CSVL_->GetXmin()) SF = SFb_CSVL_->Eval(SFb_CSVL_->GetXmin());
       }
       if(WP == "medium"){
          SF = SFb_CSVM_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVM_->GetXmax()) SF = SFb_CSVM_->Eval(SFb_CSVM_->GetXmax());
          if(jetP4->Pt() < SFb_CSVM_->GetXmin()) SF = SFb_CSVM_->Eval(SFb_CSVM_->GetXmin());
       }
       if(WP == "tight"){
          SF = SFb_CSVT_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVT_->GetXmax()) SF = SFb_CSVT_->Eval(SFb_CSVT_->GetXmax());
          if(jetP4->Pt() < SFb_CSVT_->GetXmin()) SF = SFb_CSVT_->Eval(SFb_CSVT_->GetXmin());
       }
    }
    if(fabs(flavour) == 4){
       if(WP == "loose"){
          SF = SFb_CSVL_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVL_->GetXmax()) SF = SFb_CSVL_->Eval(SFb_CSVL_->GetXmax());
          if(jetP4->Pt() < SFb_CSVL_->GetXmin()) SF = SFb_CSVL_->Eval(SFb_CSVL_->GetXmin());
       }
       if(WP == "medium"){
          SF = SFb_CSVM_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVM_->GetXmax()) SF = SFb_CSVM_->Eval(SFb_CSVM_->GetXmax());
          if(jetP4->Pt() < SFb_CSVM_->GetXmin()) SF = SFb_CSVM_->Eval(SFb_CSVM_->GetXmin());
       }
       if(WP == "tight"){
          SF = SFb_CSVT_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVT_->GetXmax()) SF = SFb_CSVT_->Eval(SFb_CSVT_->GetXmax());
          if(jetP4->Pt() < SFb_CSVT_->GetXmin()) SF = SFb_CSVT_->Eval(SFb_CSVT_->GetXmin());
       }
    }
    if(fabs(flavour) != 0 && fabs(flavour) != 5 && fabs(flavour) != 4){
       if(WP == "loose"){
          if(fabs(jetP4->Eta()) > 0. && fabs(jetP4->Eta()) < 0.5){
             SF = SFudsg_CSVL_00_05_mean_->Eval(jetP4->Pt());
             if(jetP4->Pt() > SFudsg_CSVL_00_05_mean_->GetXmax()) SF = SFudsg_CSVL_00_05_mean_->Eval(SFudsg_CSVL_00_05_mean_->GetXmax());
             if(jetP4->Pt() < SFudsg_CSVL_00_05_mean_->GetXmin()) SF = SFudsg_CSVL_00_05_mean_->Eval(SFudsg_CSVL_00_05_mean_->GetXmin());
          }
          if(fabs(jetP4->Eta()) >= 0.5 && fabs(jetP4->Eta()) < 1.){
             SF = SFudsg_CSVL_05_10_mean_->Eval(jetP4->Pt());
             if(jetP4->Pt() > SFudsg_CSVL_05_10_mean_->GetXmax()) SF = SFudsg_CSVL_05_10_mean_->Eval(SFudsg_CSVL_05_10_mean_->GetXmax());
             if(jetP4->Pt() < SFudsg_CSVL_05_10_mean_->GetXmin()) SF = SFudsg_CSVL_05_10_mean_->Eval(SFudsg_CSVL_05_10_mean_->GetXmin());
          }
          if(fabs(jetP4->Eta()) >= 1. && fabs(jetP4->Eta()) < 1.5){
             SF = SFudsg_CSVL_10_15_mean_->Eval(jetP4->Pt());
             if(jetP4->Pt() > SFudsg_CSVL_10_15_mean_->GetXmax()) SF = SFudsg_CSVL_10_15_mean_->Eval(SFudsg_CSVL_10_15_mean_->GetXmax());
             if(jetP4->Pt() < SFudsg_CSVL_10_15_mean_->GetXmin()) SF = SFudsg_CSVL_10_15_mean_->Eval(SFudsg_CSVL_10_15_mean_->GetXmin());
          }
          if(fabs(jetP4->Eta()) >= 1.5){
             SF = SFudsg_CSVL_15_24_mean_->Eval(jetP4->Pt());
             if(jetP4->Pt() > SFudsg_CSVL_15_24_mean_->GetXmax()) SF = SFudsg_CSVL_15_24_mean_->Eval(SFudsg_CSVL_15_24_mean_->GetXmax());
             if(jetP4->Pt() < SFudsg_CSVL_15_24_mean_->GetXmin()) SF = SFudsg_CSVL_15_24_mean_->Eval(SFudsg_CSVL_15_24_mean_->GetXmin());
          }
       }
       if(WP == "medium"){
          if(fabs(jetP4->Eta()) > 0. && fabs(jetP4->Eta()) < 0.8){
             SF = SFudsg_CSVM_00_08_mean_->Eval(jetP4->Pt());
             if(jetP4->Pt() > SFudsg_CSVM_00_08_mean_->GetXmax()) SF = SFudsg_CSVM_00_08_mean_->Eval(SFudsg_CSVM_00_08_mean_->GetXmax());
             if(jetP4->Pt() < SFudsg_CSVM_00_08_mean_->GetXmin()) SF = SFudsg_CSVM_00_08_mean_->Eval(SFudsg_CSVM_00_08_mean_->GetXmin());
          }
          if(fabs(jetP4->Eta()) >= 0.8 && fabs(jetP4->Eta()) < 1.6){
             SF = SFudsg_CSVM_08_16_mean_->Eval(jetP4->Pt());
             if(jetP4->Pt() > SFudsg_CSVM_08_16_mean_->GetXmax()) SF = SFudsg_CSVM_08_16_mean_->Eval(SFudsg_CSVM_08_16_mean_->GetXmax());
             if(jetP4->Pt() < SFudsg_CSVM_08_16_mean_->GetXmin()) SF = SFudsg_CSVM_08_16_mean_->Eval(SFudsg_CSVM_08_16_mean_->GetXmin());
          }
          if(fabs(jetP4->Eta()) >= 1.6){
             SF = SFudsg_CSVM_16_24_mean_->Eval(jetP4->Pt());
             if(jetP4->Pt() > SFudsg_CSVM_16_24_mean_->GetXmax()) SF = SFudsg_CSVM_16_24_mean_->Eval(SFudsg_CSVM_16_24_mean_->GetXmax());
             if(jetP4->Pt() < SFudsg_CSVM_16_24_mean_->GetXmin()) SF = SFudsg_CSVM_16_24_mean_->Eval(SFudsg_CSVM_16_24_mean_->GetXmin());
          }
       }
       if(WP == "tight"){
          SF = SFudsg_CSVT_00_24_mean_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFudsg_CSVT_00_24_mean_->GetXmax()) SF = SFudsg_CSVT_00_24_mean_->Eval(SFudsg_CSVT_00_24_mean_->GetXmax());
          if(jetP4->Pt() < SFudsg_CSVT_00_24_mean_->GetXmin()) SF = SFudsg_CSVT_00_24_mean_->Eval(SFudsg_CSVT_00_24_mean_->GetXmin());
       }
    }

    if(BUDEBUG) std::cerr << "Leaving BtagSFReader::getSF" << std::endl;
    return SF;
}
//----------------------------------------------------------------------------------------------------------------------------
float BtagSFReader::getSFErrorUp(const TLorentzVector* jetP4,const float& flavour, std::string WP)
{
    if(BUDEBUG) std::cerr << "Entering BtagSFReader::getSFErrorUp" << std::endl;
      
    float SFerr = -1001.;
    float SFMax = -1001.;
    float SF = -1001.;

    if(fabs(flavour) == 5){
       if(WP == "loose"){
          SF = SFb_CSVL_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVL_->GetXmax()) SF = SFb_CSVL_->Eval(SFb_CSVL_->GetXmax());
          if(jetP4->Pt() < SFb_CSVL_->GetXmin()) SF = SFb_CSVL_->Eval(SFb_CSVL_->GetXmin());
          SFerr = h1_SFb_CSVL_-> GetBinError(h1_SFb_CSVL_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVL_-> GetBinCenter(h1_SFb_CSVL_-> GetNbinsX())+0.5*h1_SFb_CSVL_-> GetBinWidth(h1_SFb_CSVL_-> GetNbinsX()))
             SFerr = 2*h1_SFb_CSVL_-> GetBinError(h1_SFb_CSVL_-> GetNbinsX());
          if(jetP4->Pt() < h1_SFb_CSVL_-> GetBinCenter(1)-0.5*h1_SFb_CSVL_-> GetBinWidth(1))
             SFerr = 2*h1_SFb_CSVL_-> GetBinError(1);
       }
       if(WP == "medium"){
          SF = SFb_CSVM_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVM_->GetXmax()) SF = SFb_CSVM_->Eval(SFb_CSVM_->GetXmax());
          if(jetP4->Pt() < SFb_CSVM_->GetXmin()) SF = SFb_CSVM_->Eval(SFb_CSVM_->GetXmin());
          SFerr = h1_SFb_CSVM_-> GetBinError(h1_SFb_CSVM_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVM_-> GetBinCenter(h1_SFb_CSVM_-> GetNbinsX())+0.5*h1_SFb_CSVM_-> GetBinWidth(h1_SFb_CSVM_-> GetNbinsX()))
             SFerr = 2*h1_SFb_CSVM_-> GetBinError(h1_SFb_CSVM_-> GetNbinsX());
          if(jetP4->Pt() < h1_SFb_CSVM_-> GetBinCenter(1)-0.5*h1_SFb_CSVM_-> GetBinWidth(1))
             SFerr = 2*h1_SFb_CSVM_-> GetBinError(1);
       }
       if(WP == "tight"){
          SF = SFb_CSVT_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVT_->GetXmax()) SF = SFb_CSVT_->Eval(SFb_CSVT_->GetXmax());
          if(jetP4->Pt() < SFb_CSVT_->GetXmin()) SF = SFb_CSVT_->Eval(SFb_CSVT_->GetXmin());
          SFerr = h1_SFb_CSVT_-> GetBinError(h1_SFb_CSVT_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVT_-> GetBinCenter(h1_SFb_CSVT_-> GetNbinsX())+0.5*h1_SFb_CSVT_-> GetBinWidth(h1_SFb_CSVT_-> GetNbinsX()))
             SFerr = 2*h1_SFb_CSVT_-> GetBinError(h1_SFb_CSVT_-> GetNbinsX());
          if(jetP4->Pt() < h1_SFb_CSVT_-> GetBinCenter(1)-0.5*h1_SFb_CSVT_-> GetBinWidth(1))
             SFerr = 2*h1_SFb_CSVT_-> GetBinError(1);
       }
    }
    if(fabs(flavour) == 4){
       if(WP == "loose"){
          SF = SFb_CSVL_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVL_->GetXmax()) SF = SFb_CSVL_->Eval(SFb_CSVL_->GetXmax());
          if(jetP4->Pt() < SFb_CSVL_->GetXmin()) SF = SFb_CSVL_->Eval(SFb_CSVL_->GetXmin());
          SFerr = 2*h1_SFb_CSVL_-> GetBinError(h1_SFb_CSVL_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVL_-> GetBinCenter(h1_SFb_CSVL_-> GetNbinsX())+0.5*h1_SFb_CSVL_-> GetBinWidth(h1_SFb_CSVL_-> GetNbinsX()))
             SFerr = 4*h1_SFb_CSVL_-> GetBinError(h1_SFb_CSVL_-> GetNbinsX());
          if(jetP4->Pt() < h1_SFb_CSVL_-> GetBinCenter(1)-0.5*h1_SFb_CSVL_-> GetBinWidth(1))
             SFerr = 4*h1_SFb_CSVL_-> GetBinError(1);
       }
       if(WP == "medium"){
          SF = SFb_CSVM_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVM_->GetXmax()) SF = SFb_CSVM_->Eval(SFb_CSVM_->GetXmax());
          if(jetP4->Pt() < SFb_CSVM_->GetXmin()) SF = SFb_CSVM_->Eval(SFb_CSVM_->GetXmin());
          SFerr = 2*h1_SFb_CSVM_-> GetBinError(h1_SFb_CSVM_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVM_-> GetBinCenter(h1_SFb_CSVM_-> GetNbinsX())+0.5*h1_SFb_CSVM_-> GetBinWidth(h1_SFb_CSVM_-> GetNbinsX()))
             SFerr = 4*h1_SFb_CSVM_-> GetBinError(h1_SFb_CSVM_-> GetNbinsX());
          if(jetP4->Pt() < h1_SFb_CSVM_-> GetBinCenter(1)-0.5*h1_SFb_CSVM_-> GetBinWidth(1))
             SFerr = 4*h1_SFb_CSVM_-> GetBinError(1);
       }
       if(WP == "tight"){
          SF = SFb_CSVT_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVT_->GetXmax()) SF = SFb_CSVT_->Eval(SFb_CSVT_->GetXmax());
          if(jetP4->Pt() < SFb_CSVT_->GetXmin()) SF = SFb_CSVT_->Eval(SFb_CSVT_->GetXmin());
          SFerr = 2*h1_SFb_CSVT_-> GetBinError(h1_SFb_CSVT_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVT_-> GetBinCenter(h1_SFb_CSVT_-> GetNbinsX())+0.5*h1_SFb_CSVT_-> GetBinWidth(h1_SFb_CSVT_-> GetNbinsX()))
             SFerr = 4*h1_SFb_CSVT_-> GetBinError(h1_SFb_CSVT_-> GetNbinsX());
          if(jetP4->Pt() < h1_SFb_CSVT_-> GetBinCenter(1)-0.5*h1_SFb_CSVT_-> GetBinWidth(1))
             SFerr = 4*h1_SFb_CSVT_-> GetBinError(1);
       }
    }

    if(fabs(flavour) != 0 && fabs(flavour) != 5 && fabs(flavour) != 4){
       if(WP == "loose"){
          if(fabs(jetP4->Eta()) > 0. && fabs(jetP4->Eta()) < 0.5){
             SF = SFudsg_CSVL_00_05_mean_->Eval(jetP4->Pt());
             SFMax = SFudsg_CSVL_00_05_max_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFMax);
             if(jetP4->Pt() > SFudsg_CSVL_00_05_mean_->GetXmax()){
                SF = SFudsg_CSVL_00_05_mean_->Eval(SFudsg_CSVL_00_05_mean_->GetXmax());
                SFMax = SFudsg_CSVL_00_05_max_->Eval(SFudsg_CSVL_00_05_max_->GetXmax());
                SFerr = 2*fabs(SF-SFMax);
             } 
             if(jetP4->Pt() < SFudsg_CSVL_00_05_mean_->GetXmin()){
                SF = SFudsg_CSVL_00_05_mean_->Eval(SFudsg_CSVL_00_05_mean_->GetXmin());
                SFMax = SFudsg_CSVL_00_05_max_->Eval(SFudsg_CSVL_00_05_max_->GetXmin());
                SFerr = 2*fabs(SF-SFMax);
             } 
          }
          if(fabs(jetP4->Eta()) > 0.5 && fabs(jetP4->Eta()) < 1.){
             SF = SFudsg_CSVL_05_10_mean_->Eval(jetP4->Pt());
             SFMax = SFudsg_CSVL_05_10_max_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFMax);
             if(jetP4->Pt() > SFudsg_CSVL_05_10_mean_->GetXmax()){
                SF = SFudsg_CSVL_05_10_mean_->Eval(SFudsg_CSVL_05_10_mean_->GetXmax());
                SFMax = SFudsg_CSVL_05_10_max_->Eval(SFudsg_CSVL_05_10_max_->GetXmax());
                SFerr = 2*fabs(SF-SFMax);
             } 
             if(jetP4->Pt() < SFudsg_CSVL_05_10_mean_->GetXmin()){
                SF = SFudsg_CSVL_05_10_mean_->Eval(SFudsg_CSVL_05_10_mean_->GetXmin());
                SFMax = SFudsg_CSVL_05_10_max_->Eval(SFudsg_CSVL_05_10_max_->GetXmin());
                SFerr = 2*fabs(SF-SFMax);
             } 
          }
          if(fabs(jetP4->Eta()) > 1. && fabs(jetP4->Eta()) < 1.5){
             SF = SFudsg_CSVL_10_15_mean_->Eval(jetP4->Pt());
             SFMax = SFudsg_CSVL_10_15_max_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFMax);
             if(jetP4->Pt() > SFudsg_CSVL_10_15_mean_->GetXmax()){
                SF = SFudsg_CSVL_10_15_mean_->Eval(SFudsg_CSVL_10_15_mean_->GetXmax());
                SFMax = SFudsg_CSVL_10_15_max_->Eval(SFudsg_CSVL_10_15_max_->GetXmax());
                SFerr = 2*fabs(SF-SFMax);
             } 
             if(jetP4->Pt() < SFudsg_CSVL_10_15_mean_->GetXmin()){
                SF = SFudsg_CSVL_10_15_mean_->Eval(SFudsg_CSVL_10_15_mean_->GetXmin());
                SFMax = SFudsg_CSVL_10_15_max_->Eval(SFudsg_CSVL_10_15_max_->GetXmin());
                SFerr = 2*fabs(SF-SFMax);
             } 
          }
          if(fabs(jetP4->Eta()) > 1.5 && fabs(jetP4->Eta()) < 2.4){
             SF = SFudsg_CSVL_15_24_mean_->Eval(jetP4->Pt());
             SFMax = SFudsg_CSVL_15_24_max_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFMax);
             if(jetP4->Pt() > SFudsg_CSVL_15_24_mean_->GetXmax()){
                SF = SFudsg_CSVL_15_24_mean_->Eval(SFudsg_CSVL_15_24_mean_->GetXmax());
                SFMax = SFudsg_CSVL_15_24_max_->Eval(SFudsg_CSVL_15_24_max_->GetXmax());
                SFerr = 2*fabs(SF-SFMax);
             } 
             if(jetP4->Pt() < SFudsg_CSVL_15_24_mean_->GetXmin()){
                SF = SFudsg_CSVL_15_24_mean_->Eval(SFudsg_CSVL_15_24_mean_->GetXmin());
                SFMax = SFudsg_CSVL_15_24_max_->Eval(SFudsg_CSVL_15_24_max_->GetXmin());
                SFerr = 2*fabs(SF-SFMax);
             } 
          }
          if(fabs(jetP4->Eta()) >= 2.4){
             SF = SFudsg_CSVL_15_24_mean_->Eval(jetP4->Pt());
             SFMax = SFudsg_CSVL_15_24_max_->Eval(jetP4->Pt());
             SFerr = 2*fabs(SF-SFMax);
             if(jetP4->Pt() > SFudsg_CSVL_15_24_mean_->GetXmax()){
                SF = SFudsg_CSVL_15_24_mean_->Eval(SFudsg_CSVL_15_24_mean_->GetXmax());
                SFMax = SFudsg_CSVL_15_24_max_->Eval(SFudsg_CSVL_15_24_max_->GetXmax());
                SFerr = 2*fabs(SF-SFMax);
             } 
             if(jetP4->Pt() < SFudsg_CSVL_15_24_mean_->GetXmin()){
                SF = SFudsg_CSVL_15_24_mean_->Eval(SFudsg_CSVL_15_24_mean_->GetXmin());
                SFMax = SFudsg_CSVL_15_24_max_->Eval(SFudsg_CSVL_15_24_max_->GetXmin());
                SFerr = 2*fabs(SF-SFMax);
             } 
          }  
       }
       if(WP == "medium"){
          if(fabs(jetP4->Eta()) > 0. && fabs(jetP4->Eta()) < 0.8){
             SF = SFudsg_CSVM_00_08_mean_->Eval(jetP4->Pt());
             SFMax = SFudsg_CSVM_00_08_max_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFMax);
             if(jetP4->Pt() > SFudsg_CSVM_00_08_mean_->GetXmax()){
                SF = SFudsg_CSVM_00_08_mean_->Eval(SFudsg_CSVM_00_08_mean_->GetXmax());
                SFMax = SFudsg_CSVM_00_08_max_->Eval(SFudsg_CSVM_00_08_max_->GetXmax());
                SFerr = 2*fabs(SF-SFMax);
             } 
             if(jetP4->Pt() < SFudsg_CSVM_00_08_mean_->GetXmin()){
                SF = SFudsg_CSVM_00_08_mean_->Eval(SFudsg_CSVM_00_08_mean_->GetXmin());
                SFMax = SFudsg_CSVM_00_08_max_->Eval(SFudsg_CSVM_00_08_max_->GetXmin());
                SFerr = 2*fabs(SF-SFMax);
             } 
          }
          if(fabs(jetP4->Eta()) >= 0.8 && fabs(jetP4->Eta()) < 1.6){
             SF = SFudsg_CSVM_08_16_mean_->Eval(jetP4->Pt());
             SFMax = SFudsg_CSVM_08_16_max_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFMax);
             if(jetP4->Pt() > SFudsg_CSVM_08_16_mean_->GetXmax()){
                SF = SFudsg_CSVM_08_16_mean_->Eval(SFudsg_CSVM_08_16_mean_->GetXmax());
                SFMax = SFudsg_CSVM_08_16_max_->Eval(SFudsg_CSVM_08_16_max_->GetXmax());
                SFerr = 2*fabs(SF-SFMax);
             } 
             if(jetP4->Pt() < SFudsg_CSVM_08_16_mean_->GetXmin()){
                SF = SFudsg_CSVM_08_16_mean_->Eval(SFudsg_CSVM_08_16_mean_->GetXmin());
                SFMax = SFudsg_CSVM_08_16_max_->Eval(SFudsg_CSVM_08_16_max_->GetXmin());
                SFerr = 2*fabs(SF-SFMax);
             } 
          }
          if(fabs(jetP4->Eta()) >= 1.6 && fabs(jetP4->Eta()) < 2.4){
             SF = SFudsg_CSVM_16_24_mean_->Eval(jetP4->Pt());
             SFMax = SFudsg_CSVM_16_24_max_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFMax);
             if(jetP4->Pt() > SFudsg_CSVM_16_24_mean_->GetXmax()){
                SF = SFudsg_CSVM_16_24_mean_->Eval(SFudsg_CSVM_16_24_mean_->GetXmax());
                SFMax = SFudsg_CSVM_16_24_max_->Eval(SFudsg_CSVM_16_24_max_->GetXmax());
                SFerr = 2*fabs(SF-SFMax);
             } 
             if(jetP4->Pt() < SFudsg_CSVM_16_24_mean_->GetXmin()){
                SF = SFudsg_CSVM_16_24_mean_->Eval(SFudsg_CSVM_16_24_mean_->GetXmin());
                SFMax = SFudsg_CSVM_16_24_max_->Eval(SFudsg_CSVM_16_24_max_->GetXmin());
                SFerr = 2*fabs(SF-SFMax);
             } 
          }
          if(fabs(jetP4->Eta()) >= 2.4){
             SF = SFudsg_CSVM_16_24_mean_->Eval(jetP4->Pt());
             SFMax = SFudsg_CSVM_16_24_max_->Eval(jetP4->Pt());
             SFerr = 2*fabs(SF-SFMax);
             if(jetP4->Pt() > SFudsg_CSVM_16_24_mean_->GetXmax()){
                SF = SFudsg_CSVM_16_24_mean_->Eval(SFudsg_CSVM_16_24_mean_->GetXmax());
                SFMax = SFudsg_CSVM_16_24_max_->Eval(SFudsg_CSVM_16_24_max_->GetXmax());
                SFerr = 2*fabs(SF-SFMax);
             } 
             if(jetP4->Pt() < SFudsg_CSVM_16_24_mean_->GetXmin()){
                SF = SFudsg_CSVM_16_24_mean_->Eval(SFudsg_CSVM_16_24_mean_->GetXmin());
                SFMax = SFudsg_CSVM_16_24_max_->Eval(SFudsg_CSVM_16_24_max_->GetXmin());
                SFerr = 2*fabs(SF-SFMax);
             } 
          }
       }
       if(WP == "tight"){
          if(fabs(jetP4->Eta()) < 2.4){
             SF = SFudsg_CSVT_00_24_mean_->Eval(jetP4->Pt());
             SFMax = SFudsg_CSVT_00_24_max_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFMax);
             if(jetP4->Pt() > SFudsg_CSVT_00_24_mean_->GetXmax()){
                SF = SFudsg_CSVT_00_24_mean_->Eval(SFudsg_CSVT_00_24_mean_->GetXmax());
                SFMax = SFudsg_CSVT_00_24_max_->Eval(SFudsg_CSVT_00_24_max_->GetXmax());
                SFerr = 2*fabs(SF-SFMax);
             } 
             if(jetP4->Pt() < SFudsg_CSVT_00_24_mean_->GetXmin()){
                SF = SFudsg_CSVT_00_24_mean_->Eval(SFudsg_CSVT_00_24_mean_->GetXmin());
                SFMax = SFudsg_CSVT_00_24_max_->Eval(SFudsg_CSVT_00_24_max_->GetXmin());
                SFerr = 2*fabs(SF-SFMax);
             } 
          }
          if(fabs(jetP4->Eta()) >= 2.4){
             SF = SFudsg_CSVT_00_24_mean_->Eval(jetP4->Pt());
             SFMax = SFudsg_CSVT_00_24_max_->Eval(jetP4->Pt());
             SFerr = 2*fabs(SF-SFMax);
             if(jetP4->Pt() > SFudsg_CSVT_00_24_mean_->GetXmax()){
                SF = SFudsg_CSVT_00_24_mean_->Eval(SFudsg_CSVT_00_24_mean_->GetXmax());
                SFMax = SFudsg_CSVT_00_24_max_->Eval(SFudsg_CSVT_00_24_max_->GetXmax());
                SFerr = 2*fabs(SF-SFMax);
             } 
             if(jetP4->Pt() < SFudsg_CSVT_00_24_mean_->GetXmin()){
                SF = SFudsg_CSVT_00_24_mean_->Eval(SFudsg_CSVT_00_24_mean_->GetXmin());
                SFMax = SFudsg_CSVT_00_24_max_->Eval(SFudsg_CSVT_00_24_max_->GetXmin());
                SFerr = 2*fabs(SF-SFMax);
             } 
          }
       }
    }

    if(BUDEBUG) std::cerr << "Leaving BtagSFReader::getSFErrorUp" << std::endl;
    return SFerr;
}
//----------------------------------------------------------------------------------------------------------------------------
float BtagSFReader::getSFErrorDown(const TLorentzVector* jetP4,const float& flavour, std::string WP)
{
    if(BUDEBUG) std::cerr << "Entering BtagSFReader::getSFErrorDown" << std::endl;
      
    float SFerr = -1001.;
    float SFmin = -1001.;
    float SF = -1001.;

    if(fabs(flavour) == 5){
       if(WP == "loose"){
          SF = SFb_CSVL_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVL_->GetXmax()) SF = SFb_CSVL_->Eval(SFb_CSVL_->GetXmax());
          if(jetP4->Pt() < SFb_CSVL_->GetXmin()) SF = SFb_CSVL_->Eval(SFb_CSVL_->GetXmin());
          SFerr = h1_SFb_CSVL_-> GetBinError(h1_SFb_CSVL_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVL_-> GetBinCenter(h1_SFb_CSVL_-> GetNbinsX())+0.5*h1_SFb_CSVL_-> GetBinWidth(h1_SFb_CSVL_-> GetNbinsX()))
             SFerr = 2*h1_SFb_CSVL_-> GetBinError(h1_SFb_CSVL_-> GetNbinsX());
          if(jetP4->Pt() < h1_SFb_CSVL_-> GetBinCenter(1)-0.5*h1_SFb_CSVL_-> GetBinWidth(1))
             SFerr = 2*h1_SFb_CSVL_-> GetBinError(1);
       }
       if(WP == "medium"){
          SF = SFb_CSVM_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVM_->GetXmax()) SF = SFb_CSVM_->Eval(SFb_CSVM_->GetXmax());
          if(jetP4->Pt() < SFb_CSVM_->GetXmin()) SF = SFb_CSVM_->Eval(SFb_CSVM_->GetXmin());
          SFerr = h1_SFb_CSVM_-> GetBinError(h1_SFb_CSVM_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVM_-> GetBinCenter(h1_SFb_CSVM_-> GetNbinsX())+0.5*h1_SFb_CSVM_-> GetBinWidth(h1_SFb_CSVM_-> GetNbinsX()))
             SFerr = 2*h1_SFb_CSVM_-> GetBinError(h1_SFb_CSVM_-> GetNbinsX());
          if(jetP4->Pt() < h1_SFb_CSVM_-> GetBinCenter(1)-0.5*h1_SFb_CSVM_-> GetBinWidth(1))
             SFerr = 2*h1_SFb_CSVM_-> GetBinError(h1_SFb_CSVM_-> GetNbinsX());
       }
       if(WP == "tight"){
          SF = SFb_CSVT_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVT_->GetXmax()) SF = SFb_CSVT_->Eval(SFb_CSVT_->GetXmax());
          if(jetP4->Pt() < SFb_CSVT_->GetXmin()) SF = SFb_CSVT_->Eval(SFb_CSVT_->GetXmin());
          SFerr = h1_SFb_CSVT_-> GetBinError(h1_SFb_CSVT_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVT_-> GetBinCenter(h1_SFb_CSVT_-> GetNbinsX())+0.5*h1_SFb_CSVT_-> GetBinWidth(h1_SFb_CSVT_-> GetNbinsX()))
             SFerr = 2*h1_SFb_CSVT_-> GetBinError(h1_SFb_CSVT_-> GetNbinsX());
          if(jetP4->Pt() < h1_SFb_CSVT_-> GetBinCenter(1)-0.5*h1_SFb_CSVT_-> GetBinWidth(1))
             SFerr = 2*h1_SFb_CSVT_-> GetBinError(1);
       }
    }
    if(fabs(flavour) == 4){
       if(WP == "loose"){
          SF = SFb_CSVL_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVL_->GetXmax()) SF = SFb_CSVL_->Eval(SFb_CSVL_->GetXmax());
          if(jetP4->Pt() < SFb_CSVL_->GetXmin()) SF = SFb_CSVL_->Eval(SFb_CSVL_->GetXmin());
          SFerr = 2*h1_SFb_CSVL_-> GetBinError(h1_SFb_CSVL_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVL_-> GetBinCenter(h1_SFb_CSVL_-> GetNbinsX())+0.5*h1_SFb_CSVL_-> GetBinWidth(h1_SFb_CSVL_-> GetNbinsX()))
             SFerr = 4*h1_SFb_CSVL_-> GetBinError(h1_SFb_CSVL_-> GetNbinsX());
          if(jetP4->Pt() < h1_SFb_CSVL_-> GetBinCenter(1)-0.5*h1_SFb_CSVL_-> GetBinWidth(1))
             SFerr = 4*h1_SFb_CSVL_-> GetBinError(1);
       }
       if(WP == "medium"){
          SF = SFb_CSVM_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVM_->GetXmax()) SF = SFb_CSVM_->Eval(SFb_CSVM_->GetXmax());
          if(jetP4->Pt() < SFb_CSVM_->GetXmin()) SF = SFb_CSVM_->Eval(SFb_CSVM_->GetXmin());
          SFerr = 2*h1_SFb_CSVM_-> GetBinError(h1_SFb_CSVM_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVM_-> GetBinCenter(h1_SFb_CSVM_-> GetNbinsX())+0.5*h1_SFb_CSVM_-> GetBinWidth(h1_SFb_CSVM_-> GetNbinsX()))
             SFerr = 4*h1_SFb_CSVM_-> GetBinError(h1_SFb_CSVM_-> GetNbinsX());
          if(jetP4->Pt() < h1_SFb_CSVM_-> GetBinCenter(1)-0.5*h1_SFb_CSVM_-> GetBinWidth(1))
             SFerr = 4*h1_SFb_CSVM_-> GetBinError(h1_SFb_CSVM_-> GetNbinsX());
       }
       if(WP == "tight"){
          SF = SFb_CSVT_->Eval(jetP4->Pt());
          if(jetP4->Pt() > SFb_CSVT_->GetXmax()) SF = SFb_CSVT_->Eval(SFb_CSVT_->GetXmax());
          if(jetP4->Pt() < SFb_CSVT_->GetXmin()) SF = SFb_CSVT_->Eval(SFb_CSVT_->GetXmin());
          SFerr = 2*h1_SFb_CSVT_-> GetBinError(h1_SFb_CSVT_-> FindBin(jetP4->Pt()));
          if(jetP4->Pt() > h1_SFb_CSVT_-> GetBinCenter(h1_SFb_CSVT_-> GetNbinsX())+0.5*h1_SFb_CSVT_-> GetBinWidth(h1_SFb_CSVT_-> GetNbinsX()))
             SFerr = 4*h1_SFb_CSVT_-> GetBinError(h1_SFb_CSVT_-> GetNbinsX());
          if(jetP4->Pt() < h1_SFb_CSVT_-> GetBinCenter(1)-0.5*h1_SFb_CSVT_-> GetBinWidth(1))
             SFerr = 4*h1_SFb_CSVT_-> GetBinError(1);
       }
    }

    if(fabs(flavour) != 0 && fabs(flavour) != 5 && fabs(flavour) != 4){
       if(WP == "loose"){
          if(fabs(jetP4->Eta()) > 0. && fabs(jetP4->Eta()) < 0.5){
             SF = SFudsg_CSVL_00_05_mean_->Eval(jetP4->Pt());
             SFmin = SFudsg_CSVL_00_05_min_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFmin);
             if(jetP4->Pt() > SFudsg_CSVL_00_05_mean_->GetXmax()){
                SF = SFudsg_CSVL_00_05_mean_->Eval(SFudsg_CSVL_00_05_mean_->GetXmax());
                SFmin = SFudsg_CSVL_00_05_min_->Eval(SFudsg_CSVL_00_05_min_->GetXmax());
                SFerr = 2*fabs(SF-SFmin);
             } 
          }
          if(fabs(jetP4->Eta()) > 0.5 && fabs(jetP4->Eta()) < 1.){
             SF = SFudsg_CSVL_05_10_mean_->Eval(jetP4->Pt());
             SFmin = SFudsg_CSVL_05_10_min_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFmin);
             if(jetP4->Pt() > SFudsg_CSVL_05_10_mean_->GetXmax()){
                SF = SFudsg_CSVL_05_10_mean_->Eval(SFudsg_CSVL_05_10_mean_->GetXmax());
                SFmin = SFudsg_CSVL_05_10_min_->Eval(SFudsg_CSVL_05_10_min_->GetXmax());
                SFerr = 2*fabs(SF-SFmin);
             } 
          }
          if(fabs(jetP4->Eta()) > 1. && fabs(jetP4->Eta()) < 1.5){
             SF = SFudsg_CSVL_10_15_mean_->Eval(jetP4->Pt());
             SFmin = SFudsg_CSVL_10_15_min_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFmin);
             if(jetP4->Pt() > SFudsg_CSVL_10_15_mean_->GetXmax()){
                SF = SFudsg_CSVL_10_15_mean_->Eval(SFudsg_CSVL_10_15_mean_->GetXmax());
                SFmin = SFudsg_CSVL_10_15_min_->Eval(SFudsg_CSVL_10_15_min_->GetXmax());
                SFerr = 2*fabs(SF-SFmin);
             } 
          }
          if(fabs(jetP4->Eta()) > 1.5 && fabs(jetP4->Eta()) < 2.4){
             SF = SFudsg_CSVL_15_24_mean_->Eval(jetP4->Pt());
             SFmin = SFudsg_CSVL_15_24_min_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFmin);
             if(jetP4->Pt() > SFudsg_CSVL_15_24_mean_->GetXmax()){
                SF = SFudsg_CSVL_15_24_mean_->Eval(SFudsg_CSVL_15_24_mean_->GetXmax());
                SFmin = SFudsg_CSVL_15_24_min_->Eval(SFudsg_CSVL_15_24_min_->GetXmax());
                SFerr = 2*fabs(SF-SFmin);
             } 
          }
          if(fabs(jetP4->Eta()) >= 2.4){
             SF = SFudsg_CSVL_15_24_mean_->Eval(jetP4->Pt());
             SFmin = SFudsg_CSVL_15_24_min_->Eval(jetP4->Pt());
             SFerr = 2*fabs(SF-SFmin);
             if(jetP4->Pt() > SFudsg_CSVL_15_24_mean_->GetXmax()){
                SF = SFudsg_CSVL_15_24_mean_->Eval(SFudsg_CSVL_15_24_mean_->GetXmax());
                SFmin = SFudsg_CSVL_15_24_min_->Eval(SFudsg_CSVL_15_24_min_->GetXmax());
                SFerr = 2*fabs(SF-SFmin);
             } 
          }  
       }
       if(WP == "medium"){
          if(fabs(jetP4->Eta()) > 0. && fabs(jetP4->Eta()) < 0.8){
             SF = SFudsg_CSVM_00_08_mean_->Eval(jetP4->Pt());
             SFmin = SFudsg_CSVM_00_08_min_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFmin);
             if(jetP4->Pt() > SFudsg_CSVM_00_08_mean_->GetXmax()){
                SF = SFudsg_CSVM_00_08_mean_->Eval(SFudsg_CSVM_00_08_mean_->GetXmax());
                SFmin = SFudsg_CSVM_00_08_min_->Eval(SFudsg_CSVM_00_08_min_->GetXmax());
                SFerr = 2*fabs(SF-SFmin);
             } 
          }
          if(fabs(jetP4->Eta()) >= 0.8 && fabs(jetP4->Eta()) < 1.6){
             SF = SFudsg_CSVM_08_16_mean_->Eval(jetP4->Pt());
             SFmin = SFudsg_CSVM_08_16_min_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFmin);
             if(jetP4->Pt() > SFudsg_CSVM_08_16_mean_->GetXmax()){
                SF = SFudsg_CSVM_08_16_mean_->Eval(SFudsg_CSVM_08_16_mean_->GetXmax());
                SFmin = SFudsg_CSVM_08_16_min_->Eval(SFudsg_CSVM_08_16_min_->GetXmax());
                SFerr = 2*fabs(SF-SFmin);
             } 
          }
          if(fabs(jetP4->Eta()) >= 1.6 && fabs(jetP4->Eta()) < 2.4){
             SF = SFudsg_CSVM_16_24_mean_->Eval(jetP4->Pt());
             SFmin = SFudsg_CSVM_16_24_min_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFmin);
             if(jetP4->Pt() > SFudsg_CSVM_16_24_mean_->GetXmax()){
                SF = SFudsg_CSVM_16_24_mean_->Eval(SFudsg_CSVM_16_24_mean_->GetXmax());
                SFmin = SFudsg_CSVM_16_24_min_->Eval(SFudsg_CSVM_16_24_min_->GetXmax());
                SFerr = 2*fabs(SF-SFmin);
             } 
          }
          if(fabs(jetP4->Eta()) >= 2.4){
             SF = SFudsg_CSVM_16_24_mean_->Eval(jetP4->Pt());
             SFmin = SFudsg_CSVM_16_24_min_->Eval(jetP4->Pt());
             SFerr = 2*fabs(SF-SFmin);
             if(jetP4->Pt() > SFudsg_CSVM_16_24_mean_->GetXmax()){
                SF = SFudsg_CSVM_16_24_mean_->Eval(SFudsg_CSVM_16_24_mean_->GetXmax());
                SFmin = SFudsg_CSVM_16_24_min_->Eval(SFudsg_CSVM_16_24_min_->GetXmax());
                SFerr = 2*fabs(SF-SFmin);
             } 
          }
       }
       if(WP == "tight"){
          if(fabs(jetP4->Eta()) < 2.4){
             SF = SFudsg_CSVT_00_24_mean_->Eval(jetP4->Pt());
             SFmin = SFudsg_CSVT_00_24_min_->Eval(jetP4->Pt());
             SFerr = fabs(SF-SFmin);
             if(jetP4->Pt() > SFudsg_CSVT_00_24_mean_->GetXmax()){
                SF = SFudsg_CSVT_00_24_mean_->Eval(SFudsg_CSVT_00_24_mean_->GetXmax());
                SFmin = SFudsg_CSVT_00_24_min_->Eval(SFudsg_CSVT_00_24_min_->GetXmax());
                SFerr = 2*fabs(SF-SFmin);
             } 
          }
          if(fabs(jetP4->Eta()) >= 2.4){
             SF = SFudsg_CSVT_00_24_mean_->Eval(jetP4->Pt());
             SFmin = SFudsg_CSVT_00_24_min_->Eval(jetP4->Pt());
             SFerr = 2*fabs(SF-SFmin);
             if(jetP4->Pt() > SFudsg_CSVT_00_24_mean_->GetXmax()){
                SF = SFudsg_CSVT_00_24_mean_->Eval(SFudsg_CSVT_00_24_mean_->GetXmax());
                SFmin = SFudsg_CSVT_00_24_min_->Eval(SFudsg_CSVT_00_24_min_->GetXmax());
                SFerr = 2*fabs(SF-SFmin);
             } 
          }
       }
    }

    if(BUDEBUG) std::cerr << "Leaving BtagSFReader::getSFErrorDown" << std::endl;
    return SFerr;
}
//---------------------------------------------------------------------------------------------------------------------------
//ctor
BtagEfficiencyReader::BtagEfficiencyReader(){};

void BtagEfficiencyReader::Init(const std::string name_btagEfficienciesFile)
{

    name_btagEfficienciesFile_ = name_btagEfficienciesFile;
    
    TFile* btagEfficiency_File_ = TFile::Open(name_btagEfficienciesFile_.c_str());

    h2_BTaggingEff_b_L_ = (TH2F*)btagEfficiency_File_->Get("h2_BTaggingEff_b_L");
    h2_BTaggingEff_b_M_ = (TH2F*)btagEfficiency_File_->Get("h2_BTaggingEff_b_M");
    h2_BTaggingEff_b_T_ = (TH2F*)btagEfficiency_File_->Get("h2_BTaggingEff_b_T");
    
    h2_BTaggingEff_c_L_ = (TH2F*)btagEfficiency_File_->Get("h2_BTaggingEff_c_L");
    h2_BTaggingEff_c_M_ = (TH2F*)btagEfficiency_File_->Get("h2_BTaggingEff_c_M");
    h2_BTaggingEff_c_T_ = (TH2F*)btagEfficiency_File_->Get("h2_BTaggingEff_c_T");

    h2_BTaggingEff_udsg_L_ = (TH2F*)btagEfficiency_File_->Get("h2_BTaggingEff_udsg_L");
    h2_BTaggingEff_udsg_M_ = (TH2F*)btagEfficiency_File_->Get("h2_BTaggingEff_udsg_M");
    h2_BTaggingEff_udsg_T_ = (TH2F*)btagEfficiency_File_->Get("h2_BTaggingEff_udsg_T");

    return;
}
//---------------------------------------------------------------------------------------------------------------------------
//dtor
BtagEfficiencyReader::~BtagEfficiencyReader()
{
}
//---------------------------------------------------------------------------------------------------------------------------
float BtagEfficiencyReader::getBtagEfficiency(const TLorentzVector* jetP4, std::string WP, const int& jet_flavour)
{
     if(BUDEBUG) std::cerr << "Entering BtagEfficiencyReader::getBtagEfficiency" << std::endl;

     float eff = -1001.;

     if(fabs(jet_flavour) == 5){
        if(WP == "loose") eff = h2_BTaggingEff_b_L_->GetBinContent(h2_BTaggingEff_b_L_->FindBin(jetP4->Pt(),fabs(jetP4->Eta())));
        if(WP == "medium") eff = h2_BTaggingEff_b_M_->GetBinContent(h2_BTaggingEff_b_M_->FindBin(jetP4->Pt(),fabs(jetP4->Eta())));

        if(WP == "tight") eff = h2_BTaggingEff_b_T_->GetBinContent(h2_BTaggingEff_b_T_->FindBin(jetP4->Pt(),fabs(jetP4->Eta())));
     }
     
     if(fabs(jet_flavour) == 4){
        if(WP == "loose") eff = h2_BTaggingEff_c_L_->GetBinContent(h2_BTaggingEff_c_L_->FindBin(jetP4->Pt(),fabs(jetP4->Eta())));
        if(WP == "medium") eff = h2_BTaggingEff_c_M_->GetBinContent(h2_BTaggingEff_c_M_->FindBin(jetP4->Pt(),fabs(jetP4->Eta())));
        if(WP == "tight") eff = h2_BTaggingEff_c_T_->GetBinContent(h2_BTaggingEff_c_T_->FindBin(jetP4->Pt(),fabs(jetP4->Eta())));
     }
     
     if(fabs(jet_flavour) != 0 && fabs(jet_flavour) != 4 && fabs(jet_flavour) != 5){
        if(WP == "loose") eff = h2_BTaggingEff_udsg_L_->GetBinContent(h2_BTaggingEff_udsg_L_->FindBin(jetP4->Pt(),fabs(jetP4->Eta())));
        if(WP == "medium") eff = h2_BTaggingEff_udsg_M_->GetBinContent(h2_BTaggingEff_udsg_M_->FindBin(jetP4->Pt(),fabs(jetP4->Eta())));
        if(WP == "tight") eff = h2_BTaggingEff_udsg_T_->GetBinContent(h2_BTaggingEff_udsg_T_->FindBin(jetP4->Pt(),fabs(jetP4->Eta())));
     }
     
     if(eff == 0) eff = -1001.;
     if(BUDEBUG) std::cerr << "Leaving BtagEfficiencyReader::getBtagEfficiency" << std::endl;
     return eff;
}
//---------------------------------------------------------------------------------------------------------------------------
float BtagEfficiencyReader::getBtagEfficiencyError(const TLorentzVector* jetP4, std::string WP, const int& jet_flavour)
{
     if(BUDEBUG) std::cerr << "Entering BtagEfficiencyReader::getBtagEfficiencyError" << std::endl;

     float eff_err = -1001.;

     if(fabs(jet_flavour) == 5){
        if(WP == "loose") eff_err = h2_BTaggingEff_b_L_->GetBinError(h2_BTaggingEff_b_L_->FindBin(jetP4->Pt(),fabs(jetP4->Eta())));
        if(WP == "medium") eff_err = h2_BTaggingEff_b_M_->GetBinError(h2_BTaggingEff_b_M_->FindBin(jetP4->Pt(),fabs(jetP4->Eta())));
        if(WP == "tight") eff_err = h2_BTaggingEff_b_T_->GetBinError(h2_BTaggingEff_b_T_->FindBin(jetP4->Pt(),fabs(jetP4->Eta())));
     }
     
     if(fabs(jet_flavour) == 4){
        if(WP == "loose") eff_err = h2_BTaggingEff_c_L_->GetBinError(h2_BTaggingEff_c_L_->FindBin(jetP4->Pt(),fabs(jetP4->Eta())));
        if(WP == "medium") eff_err = h2_BTaggingEff_c_M_->GetBinError(h2_BTaggingEff_c_M_->FindBin(jetP4->Pt(),fabs(jetP4->Eta())));
        if(WP == "tight") eff_err = h2_BTaggingEff_c_T_->GetBinError(h2_BTaggingEff_c_T_->FindBin(jetP4->Pt(),fabs(jetP4->Eta())));
     }
     
     if(fabs(jet_flavour) != 0 && fabs(jet_flavour) != 4 && fabs(jet_flavour) != 5){
        if(WP == "loose") eff_err = h2_BTaggingEff_udsg_L_->GetBinError(h2_BTaggingEff_udsg_L_->FindBin(jetP4->Pt(),fabs(jetP4->Eta())));
        if(WP == "medium") eff_err = h2_BTaggingEff_udsg_M_->GetBinError(h2_BTaggingEff_udsg_M_->FindBin(jetP4->Pt(),fabs(jetP4->Eta())));
        if(WP == "tight") eff_err = h2_BTaggingEff_udsg_T_->GetBinError(h2_BTaggingEff_udsg_T_->FindBin(jetP4->Pt(),fabs(jetP4->Eta())));
     }

     if(eff_err == 0) eff_err = -1001.;
     if(BUDEBUG) std::cerr << "Leaving BtagEfficiencyReader::getBtagEfficiencyError" << std::endl;
     return eff_err;
}

//---------------------------------------------------------------------------------------------------------------------------

// Functions to calculate the final event weight and error

float jetWeight(std::string WP, const float jet_SF, const float jet_eff, const float& jet_csvBtag){
      
      float csvWP = 0.;
      
      if(WP == "loose") csvWP = 0.244;
      if(WP == "medium") csvWP = 0.679;
      if(WP == "tight") csvWP = 0.898;
       
      float weight = 1.;
      
      if(jet_csvBtag > csvWP) weight = jet_SF;
      else weight = (1-jet_SF*jet_eff)/(1-jet_eff);

      if(jet_csvBtag <= csvWP && (1-jet_eff) == 0) weight = 0.;
      
      return weight;
}

//---------------------------------------------------------------------------------------------------------------------------

float jetWeight_err(std::string WP, const float jet_SF,const float jet_SF_err, const float jet_eff,const float jet_eff_err, const float& jet_csvBtag){
      
      float csvWP = 0.;
      
      if(WP == "loose") csvWP = 0.244;
      if(WP == "medium") csvWP = 0.679;
      if(WP == "tight") csvWP = 0.898;
      
      float error = 0.;

      if(jet_csvBtag > csvWP) error = jet_SF_err;
      else{
         float first = jet_eff/(1-jet_eff);
         float second = (1-jet_SF)/((1-jet_eff)*(1-jet_eff));
         error = sqrt(first*first*jet_SF_err*jet_SF_err+second*second*jet_eff_err*jet_eff_err);
      }
      
      return error;
}

//---------------------------------------------------------------------------------------------------------------------------

float eventWeight_2jets(std::string WP, const float& j1_SF, const float& j2_SF, const float& j1_eff, const float& j2_eff, const float& j1_csvBtag, const float& j2_csvBtag){

      float csvWP = 0.;
      
      if(WP == "loose") csvWP = 0.244;
      if(WP == "medium") csvWP = 0.679;
      if(WP == "tight") csvWP = 0.898;
      
      float weight1 = 1.;
      float weight2 = 1.;
      
      if(j1_csvBtag > csvWP) weight1 = j1_SF;
      else weight1 = (1-j1_SF*j1_eff)/(1-j1_eff);

      if(j2_csvBtag > csvWP) weight2 = j2_SF;
      else weight2 = (1-j2_SF*j2_eff)/(1-j2_eff);

      if(j1_csvBtag <= csvWP && (1-j1_eff) == 0) weight1 = 0.;
      if(j2_csvBtag <= csvWP && (1-j2_eff) == 0) weight2 = 0.;
      
      return weight1*weight2;
}

//---------------------------------------------------------------------------------------------------------------------------

float eventWeight_error_2jets(std::string WP, const float& j1_SF, const float& j1_SF_error, const float& j2_SF, const float& j2_SF_error, const float& j1_eff, const float& j1_eff_error, const float& j2_eff, const float& j2_eff_error, const float& j1_flavour, const float& j2_flavour, const float& j1_csvBtag, const float& j2_csvBtag){

      float csvWP = 0.;
      
      if(WP == "loose") csvWP = 0.244;
      if(WP == "medium") csvWP = 0.679;
      if(WP == "tight") csvWP = 0.898;

      float weight1 = 0.;
      float weight2 = 0.;
      
      float dW1oSF = 0.;
      float dW2oSF = 0.;
      
      float dW1oEff = 0; 
      float dW2oEff = 0;
      
      if(j1_csvBtag > csvWP) weight1 = j1_SF;
      else weight1 = (1-j1_SF*j1_eff)/(1-j1_eff);
      if(j1_csvBtag > csvWP) dW1oSF = 1.;
      else dW1oSF = -j1_eff/(1-j1_eff);
      if(j1_csvBtag > csvWP) dW1oEff = 0.;
      else dW1oEff = (1-j1_SF)/((1-j1_eff)*(1-j1_eff));

      if(j2_csvBtag > csvWP) weight2 = j2_SF;
      else weight2 = (1-j2_SF*j2_eff)/(1-j2_eff);
      if(j2_csvBtag > csvWP) dW2oSF = 1.;
      else dW2oSF = -j2_eff/(1-j2_eff);
      if(j2_csvBtag > csvWP) dW2oEff = 0.;
      else dW2oEff = (1-j2_SF)/((1-j2_eff)*(1-j2_eff));

      if(j1_csvBtag <= csvWP && (1-j1_eff) == 0){
         weight1 = 0.;
         dW1oSF = 0.;
         dW1oEff = 0.;
      }
      if(j2_csvBtag <= csvWP && (1-j2_eff) == 0){
         weight2 = 0.;
         dW2oSF = 0.;
         dW2oEff = 0.;
      }

      float rho12 = 1.;
      
      if(fabs(j1_flavour) == 5 && fabs(j2_flavour) != 5 && fabs(j2_flavour) != 4) rho12 = 0.;
      if(fabs(j1_flavour) == 4 && fabs(j2_flavour) != 5 && fabs(j2_flavour) != 4) rho12 = 0.;
      if(fabs(j2_flavour) == 5 && fabs(j1_flavour) != 5 && fabs(j1_flavour) != 4) rho12 = 0.;
      if(fabs(j2_flavour) == 4 && fabs(j1_flavour) != 5 && fabs(j1_flavour) != 4) rho12 = 0.;

      float err2 = (dW1oSF*weight2)*(dW1oSF*weight2)*j1_SF_error*j1_SF_error + 
                   (dW1oEff*weight2)*(dW1oEff*weight2)*j1_eff_error*j1_eff_error +
                   (dW2oSF*weight1)*(dW2oSF*weight1)*j2_SF_error*j2_SF_error + 
                   (dW2oEff*weight1)*(dW2oEff*weight1)*j2_eff_error*j2_eff_error +
                   2*(dW1oSF*weight2)*(dW2oSF*weight1)*j1_SF_error*j2_SF_error*rho12;

      return sqrt(err2);
}

//---------------------------------------------------------------------------------------------------------------------------

float eventWeight_3jets(std::string WP, const float& j1_SF, const float& j2_SF, const float& j3_SF, const float& j1_eff, const float& j2_eff, const float& j3_eff, const float& j1_csvBtag, const float& j2_csvBtag, const float& j3_csvBtag){

      float csvWP = 0.;
      
      if(WP == "loose") csvWP = 0.244;
      if(WP == "medium") csvWP = 0.679;
      if(WP == "tight") csvWP = 0.898;
      
      float weight1 = 1.;
      float weight2 = 1.;
      float weight3 = 1.;
     
      if(j1_csvBtag > csvWP) weight1 = j1_SF;
      else weight1 = (1-j1_SF*j1_eff)/(1-j1_eff);

      if(j2_csvBtag > csvWP) weight2 = j2_SF;
      else weight2 = (1-j2_SF*j2_eff)/(1-j2_eff);

      if(j3_csvBtag > csvWP) weight3 = j3_SF;
      else weight3 = (1-j3_SF*j3_eff)/(1-j3_eff);

      if(j1_csvBtag <= csvWP && (1-j1_eff) == 0) weight1 = 0.;
      if(j2_csvBtag <= csvWP && (1-j2_eff) == 0) weight2 = 0.;
      if(j3_csvBtag <= csvWP && (1-j3_eff) == 0) weight3 = 0.;
      
      return weight1*weight2*weight3;
}

//---------------------------------------------------------------------------------------------------------------------------

float eventWeight_error_3jets(std::string WP, const float& j1_SF, const float& j1_SF_error, const float& j2_SF, const float& j2_SF_error, const float& j3_SF, const float& j3_SF_error, const float& j1_eff, const float& j1_eff_error, const float& j2_eff, const float& j2_eff_error, const float& j3_eff, const float& j3_eff_error, const float& j1_flavour, const float& j2_flavour, const float& j3_flavour, const float& j1_csvBtag, const float& j2_csvBtag, const float& j3_csvBtag){

      float csvWP = 0.;
      
      if(WP == "loose") csvWP = 0.244;
      if(WP == "medium") csvWP = 0.679;
      if(WP == "tight") csvWP = 0.898;

      float weight1 = 0.;
      float weight2 = 0.;
      float weight3 = 0.;
      
      float dW1oSF = 0.;
      float dW2oSF = 0.;
      float dW3oSF = 0.;
      
      float dW1oEff = 0; 
      float dW2oEff = 0;
      float dW3oEff = 0;
      
      if(j1_csvBtag > csvWP) weight1 = j1_SF;
      else weight1 = (1-j1_SF*j1_eff)/(1-j1_eff);
      if(j1_csvBtag > csvWP) dW1oSF = 1.;
      else dW1oSF = -j1_eff/(1-j1_eff);
      if(j1_csvBtag > csvWP) dW1oEff = 0.;
      else dW1oEff = (1-j1_SF)/((1-j1_eff)*(1-j1_eff));

      if(j2_csvBtag > csvWP) weight2 = j2_SF;
      else weight2 = (1-j2_SF*j2_eff)/(1-j2_eff);
      if(j2_csvBtag > csvWP) dW2oSF = 1.;
      else dW2oSF = -j2_eff/(1-j2_eff);
      if(j2_csvBtag > csvWP) dW2oEff = 0.;
      else dW2oEff = (1-j2_SF)/((1-j2_eff)*(1-j2_eff));

      if(j3_csvBtag > csvWP) weight3 = j3_SF;
      else weight3 = (1-j3_SF*j3_eff)/(1-j3_eff);
      if(j3_csvBtag > csvWP) dW3oSF = 1.;
      else dW3oSF = -j3_eff/(1-j3_eff);
      if(j3_csvBtag > csvWP) dW3oEff = 0.;
      else dW3oEff = (1-j3_SF)/((1-j3_eff)*(1-j3_eff));

      if(j1_csvBtag <= csvWP && (1-j1_eff) == 0){
         weight1 = 0.;
         dW1oSF = 0.;
         dW1oEff = 0.;
      }
      if(j2_csvBtag <= csvWP && (1-j2_eff) == 0){
         weight2 = 0.;
         dW2oSF = 0.;
         dW2oEff = 0.;
      }
      if(j3_csvBtag <= csvWP && (1-j3_eff) == 0){
         weight3 = 0.;
         dW3oSF = 0.;
         dW3oEff = 0.;
      }

      float rho12 = 1.;
      float rho13 = 1.;
      float rho23 = 1.;
      
      if(fabs(j1_flavour) == 5 && fabs(j2_flavour) != 5 && fabs(j2_flavour) != 4) rho12 = 0.;
      if(fabs(j1_flavour) == 4 && fabs(j2_flavour) != 5 && fabs(j2_flavour) != 4) rho12 = 0.;
      if(fabs(j2_flavour) == 5 && fabs(j1_flavour) != 5 && fabs(j1_flavour) != 4) rho12 = 0.;
      if(fabs(j2_flavour) == 4 && fabs(j1_flavour) != 5 && fabs(j1_flavour) != 4) rho12 = 0.;

      if(fabs(j1_flavour) == 5 && fabs(j3_flavour) != 5 && fabs(j3_flavour) != 4) rho13 = 0.;
      if(fabs(j1_flavour) == 4 && fabs(j3_flavour) != 5 && fabs(j3_flavour) != 4) rho13 = 0.;
      if(fabs(j3_flavour) == 5 && fabs(j1_flavour) != 5 && fabs(j1_flavour) != 4) rho13 = 0.;
      if(fabs(j3_flavour) == 4 && fabs(j1_flavour) != 5 && fabs(j1_flavour) != 4) rho13 = 0.;
       
      if(fabs(j2_flavour) == 5 && fabs(j3_flavour) != 5 && fabs(j3_flavour) != 4) rho23 = 0.;
      if(fabs(j2_flavour) == 4 && fabs(j3_flavour) != 5 && fabs(j3_flavour) != 4) rho23 = 0.;
      if(fabs(j3_flavour) == 5 && fabs(j2_flavour) != 5 && fabs(j2_flavour) != 4) rho23 = 0.;
      if(fabs(j3_flavour) == 4 && fabs(j2_flavour) != 5 && fabs(j2_flavour) != 4) rho23 = 0.;
      
      float err2 = (dW1oSF*weight2*weight3)*(dW1oSF*weight2*weight3)*j1_SF_error*j1_SF_error + 
                   (dW1oEff*weight2*weight3)*(dW1oEff*weight2*weight3)*j1_eff_error*j1_eff_error +
                   (dW2oSF*weight1*weight3)*(dW2oSF*weight1*weight3)*j2_SF_error*j2_SF_error + 
                   (dW2oEff*weight1*weight3)*(dW2oEff*weight1*weight3)*j2_eff_error*j2_eff_error +
                   (dW3oSF*weight1*weight2)*(dW3oSF*weight1*weight2)*j3_SF_error*j3_SF_error + 
                   (dW3oEff*weight1*weight2)*(dW3oEff*weight1*weight2)*j3_eff_error*j3_eff_error +
                   2*(dW1oSF*weight2*weight3)*(dW2oSF*weight1*weight3)*j1_SF_error*j2_SF_error*rho12 +
                   2*(dW1oSF*weight2*weight3)*(dW3oSF*weight1*weight2)*j1_SF_error*j3_SF_error*rho13 +
                   2*(dW2oSF*weight1*weight3)*(dW3oSF*weight1*weight2)*j2_SF_error*j3_SF_error*rho23;

      return sqrt(err2);
}

//---------------------------------------------------------------------------------------------------------------------------

float eventWeight_4jets(std::string WP, const float& j1_SF, const float& j2_SF, const float& j3_SF, const float& j4_SF, const float& j1_eff, const float& j2_eff, const float& j3_eff, const float& j4_eff, const float& j1_csvBtag, const float& j2_csvBtag, const float& j3_csvBtag, const float& j4_csvBtag){

      float csvWP = 0.;
      
      if(WP == "loose") csvWP = 0.244;
      if(WP == "medium") csvWP = 0.679;
      if(WP == "tight") csvWP = 0.898;
      
      float weight1 = 1.;
      float weight2 = 1.;
      float weight3 = 1.;
      float weight4 = 1.;

      if(j1_csvBtag > csvWP) weight1 = j1_SF;
      else weight1 = (1-j1_SF*j1_eff)/(1-j1_eff);

      if(j2_csvBtag > csvWP) weight2 = j2_SF;
      else weight2 = (1-j2_SF*j2_eff)/(1-j2_eff);

      if(j3_csvBtag > csvWP) weight3 = j3_SF;
      else weight3 = (1-j3_SF*j3_eff)/(1-j3_eff);

      if(j4_csvBtag > csvWP) weight4 = j4_SF;
      else weight4 = (1-j4_SF*j4_eff)/(1-j4_eff);

      if(j1_csvBtag <= csvWP && (1-j1_eff) == 0) weight1 = 0.;
      if(j2_csvBtag <= csvWP && (1-j2_eff) == 0) weight2 = 0.;
      if(j3_csvBtag <= csvWP && (1-j3_eff) == 0) weight3 = 0.;
      if(j4_csvBtag <= csvWP && (1-j4_eff) == 0) weight4 = 0.;

      return weight1*weight2*weight3*weight4;
}

//---------------------------------------------------------------------------------------------------------------------------

float eventWeight_error_4jets(std::string WP, const float& j1_SF, const float& j1_SF_error, const float& j2_SF, const float& j2_SF_error, const float& j3_SF, const float& j3_SF_error, const float& j4_SF, const float& j4_SF_error, const float& j1_eff, const float& j1_eff_error, const float& j2_eff, const float& j2_eff_error, const float& j3_eff, const float& j3_eff_error, const float& j4_eff, const float& j4_eff_error, const float& j1_flavour, const float& j2_flavour, const float& j3_flavour, const float& j4_flavour, const float& j1_csvBtag, const float& j2_csvBtag, const float& j3_csvBtag, const float& j4_csvBtag){

      float csvWP = 0.;
      
      if(WP == "loose") csvWP = 0.244;
      if(WP == "medium") csvWP = 0.679;
      if(WP == "tight") csvWP = 0.898;

      float weight1 = 0.;
      float weight2 = 0.;
      float weight3 = 0.;
      float weight4 = 0.;
      
      float dW1oSF = 0.;
      float dW2oSF = 0.;
      float dW3oSF = 0.;
      float dW4oSF = 0.;

      float dW1oEff = 0; 
      float dW2oEff = 0;
      float dW3oEff = 0;
      float dW4oEff = 0; 
      
      if(j1_csvBtag > csvWP) weight1 = j1_SF;
      else weight1 = (1-j1_SF*j1_eff)/(1-j1_eff);
      if(j1_csvBtag > csvWP) dW1oSF = 1.;
      else dW1oSF = -j1_eff/(1-j1_eff);
      if(j1_csvBtag > csvWP) dW1oEff = 0.;
      else dW1oEff = (1-j1_SF)/((1-j1_eff)*(1-j1_eff));

      if(j2_csvBtag > csvWP) weight2 = j2_SF;
      else weight2 = (1-j2_SF*j2_eff)/(1-j2_eff);
      if(j2_csvBtag > csvWP) dW2oSF = 1.;
      else dW2oSF = -j2_eff/(1-j2_eff);
      if(j2_csvBtag > csvWP) dW2oEff = 0.;
      else dW2oEff = (1-j2_SF)/((1-j2_eff)*(1-j2_eff));

      if(j3_csvBtag > csvWP) weight3 = j3_SF;
      else weight3 = (1-j3_SF*j3_eff)/(1-j3_eff);
      if(j3_csvBtag > csvWP) dW3oSF = 1.;
      else dW3oSF = -j3_eff/(1-j3_eff);
      if(j3_csvBtag > csvWP) dW3oEff = 0.;
      else dW3oEff = (1-j3_SF)/((1-j3_eff)*(1-j3_eff));

      if(j4_csvBtag > csvWP) weight4 = j4_SF;
      else weight4 = (1-j4_SF*j4_eff)/(1-j4_eff);
      if(j4_csvBtag > csvWP) dW4oSF = 1.;
      else dW4oSF = -j4_eff/(1-j4_eff);
      if(j4_csvBtag > csvWP) dW4oEff = 0.;
      else dW4oEff = (1-j4_SF)/((1-j4_eff)*(1-j4_eff));

      if(j1_csvBtag <= csvWP && (1-j1_eff) == 0){
         weight1 = 0.;
         dW1oSF = 0.;
         dW1oEff = 0.;
      }
      if(j2_csvBtag <= csvWP && (1-j2_eff) == 0){
         weight2 = 0.;
         dW2oSF = 0.;
         dW2oEff = 0.;
      }
      if(j3_csvBtag <= csvWP && (1-j3_eff) == 0){
         weight3 = 0.;
         dW3oSF = 0.;
         dW3oEff = 0.;
      }
      if(j4_csvBtag <= csvWP && (1-j4_eff) == 0){
         weight4 = 0.;
         dW4oSF = 0.;
         dW4oEff = 0.;
      }
      
      float rho12 = 1.;
      float rho13 = 1.;
      float rho14 = 1.;
      float rho23 = 1.;
      float rho24 = 1.;
      float rho34 = 1.;
       
      if(fabs(j1_flavour) == 5 && fabs(j2_flavour) != 5 && fabs(j2_flavour) != 4) rho12 = 0.;
      if(fabs(j1_flavour) == 4 && fabs(j2_flavour) != 5 && fabs(j2_flavour) != 4) rho12 = 0.;
      if(fabs(j2_flavour) == 5 && fabs(j1_flavour) != 5 && fabs(j1_flavour) != 4) rho12 = 0.;
      if(fabs(j2_flavour) == 4 && fabs(j1_flavour) != 5 && fabs(j1_flavour) != 4) rho12 = 0.;

      if(fabs(j1_flavour) == 5 && fabs(j3_flavour) != 5 && fabs(j3_flavour) != 4) rho13 = 0.;
      if(fabs(j1_flavour) == 4 && fabs(j3_flavour) != 5 && fabs(j3_flavour) != 4) rho13 = 0.;
      if(fabs(j3_flavour) == 5 && fabs(j1_flavour) != 5 && fabs(j1_flavour) != 4) rho13 = 0.;
      if(fabs(j3_flavour) == 4 && fabs(j1_flavour) != 5 && fabs(j1_flavour) != 4) rho13 = 0.;
       
      if(fabs(j1_flavour) == 5 && fabs(j4_flavour) != 5 && fabs(j4_flavour) != 4) rho14 = 0.;
      if(fabs(j1_flavour) == 4 && fabs(j4_flavour) != 5 && fabs(j4_flavour) != 4) rho14 = 0.;
      if(fabs(j4_flavour) == 5 && fabs(j1_flavour) != 5 && fabs(j1_flavour) != 4) rho14 = 0.;
      if(fabs(j4_flavour) == 4 && fabs(j1_flavour) != 5 && fabs(j1_flavour) != 4) rho14 = 0.;
      
      if(fabs(j2_flavour) == 5 && fabs(j3_flavour) != 5 && fabs(j3_flavour) != 4) rho23 = 0.;
      if(fabs(j2_flavour) == 4 && fabs(j3_flavour) != 5 && fabs(j3_flavour) != 4) rho23 = 0.;
      if(fabs(j3_flavour) == 5 && fabs(j2_flavour) != 5 && fabs(j2_flavour) != 4) rho23 = 0.;
      if(fabs(j3_flavour) == 4 && fabs(j2_flavour) != 5 && fabs(j2_flavour) != 4) rho23 = 0.;
      
      if(fabs(j2_flavour) == 5 && fabs(j4_flavour) != 5 && fabs(j4_flavour) != 4) rho24 = 0.;
      if(fabs(j2_flavour) == 4 && fabs(j4_flavour) != 5 && fabs(j4_flavour) != 4) rho24 = 0.;
      if(fabs(j4_flavour) == 5 && fabs(j2_flavour) != 5 && fabs(j2_flavour) != 4) rho24 = 0.;
      if(fabs(j4_flavour) == 4 && fabs(j2_flavour) != 5 && fabs(j2_flavour) != 4) rho24 = 0.;
      
      if(fabs(j3_flavour) == 5 && fabs(j4_flavour) != 5 && fabs(j4_flavour) != 4) rho34 = 0.;
      if(fabs(j3_flavour) == 4 && fabs(j4_flavour) != 5 && fabs(j4_flavour) != 4) rho34 = 0.;
      if(fabs(j4_flavour) == 5 && fabs(j3_flavour) != 5 && fabs(j3_flavour) != 4) rho34 = 0.;
      if(fabs(j4_flavour) == 4 && fabs(j3_flavour) != 5 && fabs(j3_flavour) != 4) rho34 = 0.;
     
      float err2 = (dW1oSF*weight2*weight3*weight4)*(dW1oSF*weight2*weight3*weight4)*j1_SF_error*j1_SF_error + 
                   (dW1oEff*weight2*weight3*weight4)*(dW1oEff*weight2*weight3*weight4)*j1_eff_error*j1_eff_error +
                   (dW2oSF*weight1*weight3*weight4)*(dW2oSF*weight1*weight3*weight4)*j2_SF_error*j2_SF_error + 
                   (dW2oEff*weight1*weight3*weight4)*(dW2oEff*weight1*weight3*weight4)*j2_eff_error*j2_eff_error +
                   (dW3oSF*weight1*weight2*weight4)*(dW3oSF*weight1*weight2*weight4)*j3_SF_error*j3_SF_error + 
                   (dW3oEff*weight1*weight2*weight4)*(dW3oEff*weight1*weight2*weight4)*j3_eff_error*j3_eff_error +
                   (dW4oSF*weight1*weight2*weight3)*(dW4oSF*weight1*weight2*weight3)*j4_SF_error*j4_SF_error + 
                   (dW4oEff*weight1*weight2*weight3)*(dW4oEff*weight1*weight2*weight3)*j4_eff_error*j4_eff_error +
                   2*(dW1oSF*weight2*weight3*weight4)*(dW2oSF*weight1*weight3*weight4)*j1_SF_error*j2_SF_error*rho12 +
                   2*(dW1oSF*weight2*weight3*weight4)*(dW3oSF*weight1*weight2*weight4)*j1_SF_error*j3_SF_error*rho13 + 
                   2*(dW1oSF*weight2*weight3*weight4)*(dW4oSF*weight1*weight2*weight3)*j1_SF_error*j4_SF_error*rho14 +
                   2*(dW2oSF*weight1*weight3*weight4)*(dW3oSF*weight1*weight2*weight4)*j2_SF_error*j3_SF_error*rho23 +
                   2*(dW2oSF*weight1*weight3*weight4)*(dW4oSF*weight1*weight2*weight3)*j2_SF_error*j4_SF_error*rho24 +
                   2*(dW3oSF*weight1*weight2*weight4)*(dW4oSF*weight1*weight2*weight3)*j3_SF_error*j4_SF_error*rho34;

      return sqrt(err2);
}

//---------------------------------------------------------------------------------------------------------------------------
