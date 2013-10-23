#ifndef __BTAGUTILS__
#define __BTAGUTILS__

#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>

class JetFlavourReader
{
public:
	
    JetFlavourReader(const std::string name_JetFlavourFile); 
    virtual ~JetFlavourReader();
    int getJetFlavour(const int& lumis, const int& event, const TLorentzVector* jetP4 );
    
	
protected:
    
    std::string name_JetFlavourFile_;
    std::map<int,std::map<int,std::vector<TLorentzVector> > > AODjet_p4_;
    std::map<int,std::map<int,std::vector<int> > > AODjet_Flav_;
    std::vector<float> dR_;
    std::map<float,int> dR_map_;
    
};


class BtagSFReader
{
public:
	
    BtagSFReader(const std::string name_btagSFFile); 
    virtual ~BtagSFReader();
    float getSF(const TLorentzVector* jetP4,const float& flavour, std::string WP);
    float getSFErrorUp(const TLorentzVector* jetP4,const float& flavour, std::string WP);	
    float getSFErrorDown(const TLorentzVector* jetP4,const float& flavour, std::string WP);
	
protected:
    
    std::string name_btagSFFile_;
    
    TF1*  SFb_CSVL_;
    TH1F* h1_SFb_CSVL_;

    TF1*  SFb_CSVM_;
    TH1F* h1_SFb_CSVM_;

    TF1*  SFb_CSVT_;
    TH1F* h1_SFb_CSVT_;

    TF1* SFudsg_CSVL_00_05_max_;
    TF1* SFudsg_CSVL_00_05_mean_;
    TF1* SFudsg_CSVL_00_05_min_;
    
    TF1* SFudsg_CSVL_05_10_max_;
    TF1* SFudsg_CSVL_05_10_mean_;
    TF1* SFudsg_CSVL_05_10_min_;

    TF1* SFudsg_CSVL_10_15_max_;
    TF1* SFudsg_CSVL_10_15_mean_;
    TF1* SFudsg_CSVL_10_15_min_;

    TF1* SFudsg_CSVL_15_24_max_;
    TF1* SFudsg_CSVL_15_24_mean_;
    TF1* SFudsg_CSVL_15_24_min_;

    TF1* SFudsg_CSVM_00_08_max_;
    TF1* SFudsg_CSVM_00_08_mean_;
    TF1* SFudsg_CSVM_00_08_min_;

    TF1* SFudsg_CSVM_08_16_max_;
    TF1* SFudsg_CSVM_08_16_mean_;
    TF1* SFudsg_CSVM_08_16_min_;

    TF1* SFudsg_CSVM_16_24_max_;
    TF1* SFudsg_CSVM_16_24_mean_;
    TF1* SFudsg_CSVM_16_24_min_;

    TF1* SFudsg_CSVT_00_24_max_;
    TF1* SFudsg_CSVT_00_24_mean_;
    TF1* SFudsg_CSVT_00_24_min_;
    
};

class BtagEfficiencyReader
{
public:
	
    BtagEfficiencyReader(const std::string name_btagEfficienciesFile); 
    virtual ~BtagEfficiencyReader();
    float getBtagEfficiency(const TLorentzVector* jetP4, std::string WP, const int& jet_flavour);
    float getBtagEfficiencyError(const TLorentzVector* jetP4, std::string WP, const int& jet_flavour);
    
	
protected:
    
    std::string name_btagEfficienciesFile_;

    TH2F* h2_BTaggingEff_b_L_;
    TH2F* h2_BTaggingEff_b_M_;
    TH2F* h2_BTaggingEff_b_T_;

    TH2F* h2_BTaggingEff_c_L_;
    TH2F* h2_BTaggingEff_c_M_;
    TH2F* h2_BTaggingEff_c_T_;

    TH2F* h2_BTaggingEff_udsg_L_;
    TH2F* h2_BTaggingEff_udsg_M_;
    TH2F* h2_BTaggingEff_udsg_T_;
    
};

// Functions to calculate the final event weight and error

float jetWeight(std::string WP, const float jet_SF, const float jet_eff, const float& jet_csvBtag);

float jetWeight_err(std::string WP, const float jet_SF,const float jet_SF_err, const float jet_eff,const float jet_eff_err, const float& jet_csvBtag);

float eventWeight(std::string WP, const float j1_SF, const float j2_SF, const float j1_eff, const float j2_eff, const float& j1_csvBtag, const float& j2_csvBtag);

float eventWeight_error(std::string WP, const float j1_SF, const float j1_SF_error, const float j2_SF, const float j2_SF_error , const float j1_eff, const float j1_eff_err, const float j2_eff, const float j2_eff_err , const float j1_flavour, const float j2_flavour, const float& j1_csvBtag, const float& j2_csvBtag);

#endif
