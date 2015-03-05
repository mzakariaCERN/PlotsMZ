////===  author : Pelin Kurt Garberson ======/////
////===  Jet Shape and Jet-track correlation analyses ===///
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include <TF1.h>

#include "assert.h"
#include <fstream>
#include "TMath.h"
#include "TVector3.h"
//#include "hiForest.h"
//#include "commonSetup.h"
#include "class_def/mixing_tree.h"


#include "TH2F.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TCut.h>
#include <TLegend.h>

using namespace std;

//MZ// define is another way to set a constant. 
#define nCBins 4
//#define nCBins 2
#define nPtBins_sub 10
//#define nPtBins 17
//#include "TrackingCorrections2012.h"

//#define nPtBins 4
#define nPtBins 1


//#define nTrkPtBins 1
//#define nTrkPtBins 4
//#define nTrkPtBins 6
#define nTrkPtBins 5
#define nCentralityBins 40

//float trkPtCut=0.5;
float trkPtCut=1.0;
int parti = -999;

bool is_data = false;
bool test_it = false;
bool is_selected = true;
bool is_dijet = false;
bool is_leadingjet = false;
bool is_subleadingjet = false;


enum enum_dataset_types {e_Data2011,e_HydJet30,e_HydJet50,e_HydJet80, e_HydJet100, e_HydJet120,e_HydJet170,e_HydJet200,e_HydJet250, e_HydJet300, e_n_dataset_types};

TString dataset_type_strs[e_n_dataset_types] = { "Data2011","HydJet30","HydJet50","HydJet80", "HydJet100", "HydJet120","HydJet170","HydJet200","HydJet250","HydJet300" };
int dataset_type_code = -999;




enum enum_npart_types {
    e_0_10, e_10_20, e_20_30, e_30_40, e_40_50, e_50_60, e_60_70, e_70_80,e_80_90,e_90_100,
    e_100_110, e_110_120, e_120_130, e_130_140, e_140_150, e_150_160, e_160_170, e_170_180,e_180_190,e_190_200,
    e_200_210, e_210_220, e_220_230, e_230_240, e_240_250, e_250_260, e_260_270, e_270_280,e_280_290,e_290_300,
    e_300_310, e_310_320, e_320_330, e_330_340, e_340_350, e_350_360, e_360_370, e_370_380,e_380_390,e_390_400,
    e_400_410, e_410_420, e_420_430, e_430_440, e_440_450, e_450_460, e_460_70, e_470_480,e_480_490,e_490_500,
    e_500_510, e_510_520, e_520_530, e_530_540, e_540_550, e_550_560, e_560_570, e_570_580,e_580_590,e_590_600,
    e_600_610, e_610_620, e_620_630, e_630_640, e_640_650, e_650_660, e_660_670, e_670_680,e_680_690,e_690_700,
    e_700_710, e_710_720, e_720_730, e_730_740, e_740_750, e_750_760, e_760_770, e_770_780,e_780_790,e_790_800,
    e_800_810, e_810_820, e_820_830, e_830_840, e_840_850, e_850_860, e_860_870,e_870_880, e_880_890, e_890_900,
    e_900_810, e_910_920, e_920_830, e_930_940, e_940_950, e_950_960, e_960_970,e_970_980, e_980_990, e_990_1000,
    e_1000_1010, e_1010_1020, e_1020_1030, e_1030_1040, e_1040_1050, e_1050_1060, e_1060_1070,e_1070_1080, e_1080_1090, e_1090_2000,
    e_2000_2010, e_2010_2020, e_2020_2030, e_2030_2040, e_2040_2050, e_2050_2060, e_2060_2070,e_2070_2080, e_2080_2090, e_2090_2100,
    e_2100_2110, e_2110_2120, e_2120_2130, e_2130_2140, e_2140_2150, e_2150_2160, e_2160_2170,e_2170_2180, e_2180_2190, e_2190_2200, e_AllParts, e_n_npart_types};

TString npart_type_strs[e_n_npart_types] = {
    "part1",  "part2", "part3", "part4", "part5", "part6", "part7", "part8", "part9", "part10",
    "part11", "part12","part13","part14","part15","part16","part17","part18","part19","part20",
    "part21", "part22","part23","part24","part25","part26","part27", "part28", "part29", "part30",
    "part31", "part32", "part33", "part34", "part35", "part36", "part37", "part38", "part39", "part40",
    "part41",  "part42", "part43", "part44", "part45", "part46", "part47", "part48", "part49", "part50",
    "part51",  "part52", "part53", "part54", "part55", "part56", "part57", "part58", "part59", "part60",
    "part61",  "part62", "part63", "part64", "part65", "part66", "part67", "part68", "part69", "part70",
    "part71",  "part72", "part73", "part74", "part75", "part76", "part77", "part78", "part79", "part80",
    "part81",  "part82", "part83", "part84", "part85", "part86", "part87", "part88", "part89", "part90",
    "part91",  "part92", "part93", "part94", "part95", "part96", "part97", "part98", "part99", "part100",
    "part101",  "part102", "part103", "part104", "part105", "part106", "part107", "part108", "part109", "part110",
    "part111",  "part112", "part113", "part114", "part115", "part116", "part117", "part118", "part119", "part120",
    "part121",  "part122", "part123", "part124", "part125", "part126", "part127", "part128", "part129", "part130", "AllParts" };


int npart = -999;



float CBins[nCBins+1] = {0, 20, 60, 100, 200}; //MZ silly naming convention
TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30", "Cent50", "Cent100"};

float CentralityBins[nCentralityBins+1] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200};  //according to binning in PAS


float PtBins[nPtBins+1] = {0, 300};
TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};


float TrkPtBins[nTrkPtBins+1] = {1, 2, 3, 4, 8, 999};
TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt999" };

//float TrkPtBins[nTrkPtBins+1] = {0.5, 999};
//TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt05", "TrkPt999" };







//float TrkPtBins[nTrkPtBins+1] = {1, 999};
//TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt1","TrkPt999"};

//	float PtBins[nPtBins+1] = {100, 110, 120, 130, 150, 200, 300};
//	TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt110", "Pt120","Pt130", "Pt150", "Pt200","Pt300"};


//float PtBins[nPtBins+1] = {100,  120,  150, 200, 300};
//TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt120","Pt150", "Pt200","Pt300"};


//float PtBins[nPtBins+1] = {100, 120, 140, 160, 180, 200,220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 500};
//TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt120", "Pt140","Pt160","Pt180", "Pt200" ,"Pt220" ,"Pt240", "Pt260", "Pt280", "Pt300", "Pt320","Pt340","Pt360","Pt380", "Pt400","Pt500"};


//float PtBins[nPtBins+1] = {100, 110, 120, 130, 140, 150, 160, 200, 300, 500};
//TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt110", "Pt120","Pt130","Pt140" , "Pt150", "Pt160", "Pt200","Pt300", "Pt500"};
Double_t x_Bins_mc[12] = {30,50,70,90,120,160,200,240,280,340,400};



const int nPtBintrk = 27;
Double_t ptBintrk[nPtBintrk+1] = {0, 0.5, 1, 1.203915, 1.449412, 1.74497, 2.100796, 2.529181, 3.04492, 3.665826, 4.413344, 5.313293, 6.396755, 7.701152, 9.271536, 11.16214, 13.43828, 16.17855, 19.47761, 23.44939, 28.23108, 33.98783, 40.91848, 49.26238, 59.30774, 71.40151, 85.96137, 103.4902};





float PtBins_sub[nPtBins_sub+1] = {50, 60, 70, 80, 100, 120, 140, 160, 200, 300, 500};
TString PtBin_sub_strs[nPtBins_sub+1] = {"Pt50", "Pt60", "Pt70", "Pt80",  "Pt100", "Pt120", "Pt140","Pt160", "Pt200" , "Pt300", "Pt500"};


class hist_class {
  public:
    hist_class(TString the_desc, bool is_it_data);
    void Delete();
    void NormalizeByCounts();
    void Write();
    void AddHists(hist_class *more_hists, float wt);
    bool is_data;
    // bool is_selected;
    TString desc;
    int n_evt_raw ;

    TH1F* NumberOfMatches;
    TH1F* NEvents;
    TH1F* NEvents_test;
    TH1F* NEvents_after_noise;
    TH1F* NEvents_after_spike;
    TH1F* NEvents_after_dphi;
    TH1F* NEvents_before_dphi;
    TH1F* NEvents_dijets;
    TH1F* NEvents_after_trigger;


    TH1F* AllJetPt_raw_hist;
    TH1F* AllJetPhi_hist;
    TH1F* AllJetEta_hist;
    TH1F* AllJetPt_hist;
    TH1F* First_AllJetPhi_hist;
    TH1F* First_AllJetEta_hist;
    TH1F* First_AllJetPt_hist;
    TH1F* Sub_AllJetPhi_hist;
    TH1F* Sub_AllJetEta_hist;
    TH1F* Sub_AllJetPt_hist;
    TH1F* JetPt_fraction;
    TH1F* Centrality;
    TH1F* track_vz;
    TH1F* track_vz_weighted;
    TH1F* Centrality_weighted;

    TH1F* jet_pT_hist[nCBins];
    TH1F* jet_phi_hist[nCBins];
    TH1F* jet_eta_hist[nCBins];
    TH1F* jet_corrpT_hist[nCBins];
    TH1F* LeadingJetPt_hist[nCBins];
    TH1F* LeadingJetPhi_hist[nCBins];
    TH1F* LeadingJetEta_hist[nCBins];
    TH1F* SubJetPt_hist[nCBins];
    TH1F* SubJetPhi_hist[nCBins];
    TH1F* SubJetEta_hist[nCBins];

    TH1F* ThirdJetPt_hist[nCBins][nPtBins];
    TH1F* ThirdJetPhi_hist[nCBins][nPtBins];
    TH1F* ThirdJetEta_hist[nCBins][nPtBins];
    TH1F* all_cand_pT_hist[nCBins][nPtBins];
    TH1F* all_cand_phi_hist[nCBins][nPtBins];
    TH1F* all_cand_eta_hist[nCBins][nPtBins];

    TH1D* track_cand_pT_hist[nCBins][nPtBins];
    TH1D* track_cand_phi_hist[nCBins][nPtBins];
    TH1D* track_cand_eta_hist[nCBins][nPtBins];
    TH1F* track_cand_pT_hist_subleadingjet[nCBins][nPtBins];

    TH1F* genparticles_bkg_pT_hist[nCBins][nPtBins];
    TH1F* genparticles_pT_hist[nCBins][nPtBins];
    TH1F* gensube_hist[nCBins][nPtBins];

   //TH1F* dPhi_hist[nCBins][nPtBins];
    TH1F* dPhi_hist[nCBins];

    TH1F* dEta_hist[nCBins][nPtBins];
    TH1D* track_bkg_eta_hist_weighted[nCBins][nPtBins];
    TH1D* track_bkg_eta_hist[nCBins][nPtBins];
    TH1D* track_gen_bkg_eta_hist[nCBins][nPtBins];
    TH1F* alljets_rec;



    //TH1F* dPhi_after_hist[nCBins][nPtBins];
    TH1F* dPhi_after_hist[nCBins];
    TH1D* dPhi_vector[nCBins];

    TH1F* dPhi_leadingjet_hist[nCBins];
    TH1F* dPhi_subleadingjet_hist[nCBins];

    TH1F* all_jets_corrpT[nCBins][nPtBins];
    TH1F* all_jets_phi[nCBins][nPtBins];
    TH1F* all_jets_eta[nCBins][nPtBins];
    TH1F* only_leadingjets_corrpT[nCBins][nPtBins];
    TH1F* only_leadingjets_phi[nCBins][nPtBins];
    TH1F* only_leadingjets_eta[nCBins][nPtBins];
    TH1F* only_subleadingjets_corrpT[nCBins][nPtBins];
    TH1F* only_subleadingjets_phi[nCBins][nPtBins];
    TH1F* only_subleadingjets_eta[nCBins][nPtBins];
    TH1F* neutral_cand_pT_hist[nCBins][nPtBins];
    TH1F* neutral_cand_phi_hist[nCBins][nPtBins];
    TH1F* neutral_cand_eta_hist[nCBins][nPtBins];
    TH1F* photons_cand_pT_hist[nCBins][nPtBins];
    TH1F* photons_cand_phi_hist[nCBins][nPtBins];
    TH1F* photons_cand_eta_hist[nCBins][nPtBins];

    TH1F* NumNeutral[nCBins][nPtBins];
    TH1F* NumPhotons[nCBins][nPtBins];
    TH1F* NumAll[nCBins][nPtBins];
    TH1F* NumCharged[nCBins][nPtBins];

    TH1F* NumAll_bkg[nCBins][nPtBins];
    TH1F* NumNeutral_bkg[nCBins][nPtBins];
    TH1F* NumPhotons_bkg[nCBins][nPtBins];
    TH1F* NumChargedHadrons_bkg[nCBins][nPtBins];
    TH1F* NumChargedParticles_bkg[nCBins][nPtBins];
    TH1F* NumElectrons_bkg[nCBins][nPtBins];
    TH1F* NumMuons_bkg[nCBins][nPtBins];

    TH1F* radius_hist[nCBins][nPtBins];
    TH1F* radius_hist_mine[nCBins][nPtBins];
    TH1F* radius_hist_mine2[nCBins][nPtBins];
    TH1F* radius_hist_mine3[nCBins][nPtBins];

    TH1F* SumJetPtFraction_hist[nCBins][nPtBins];

    //TProfile *JetShapeIntegratedParticles[nCBins][nPtBins];
    // TProfile *JetShapeDiffParticles[nCBins][nPtBins];

    TH1F* JetShapeIntegratedParticles[nCBins][nPtBins];
    TH1F* JetShapeDiffParticles[nCBins][nPtBins];

    TH1F* Multiplicity_vs_radius[nCBins][nPtBins];
    TH1F* Multiplicity_vs_radius_bkg[nCBins][nPtBins];

    TH1F* Multiplicity_vs_radius_gen[nCBins][nPtBins];
    TH1F* Multiplicity_vs_radius_gen_bkg[nCBins][nPtBins];

    TH1F* Multiplicity_vs_radius_gen_all[nCBins][nPtBins];
    TH1F* Multiplicity_vs_radius_gen_all_bkg[nCBins][nPtBins];


    TH1F* Multiplicity_vs_radius_ptwt[nCBins][nPtBins];
    TH1F* Multiplicity_vs_radius_bkg_ptwt[nCBins][nPtBins];

    TH1F* Multiplicity_vs_radius_gen_ptwt[nCBins][nPtBins];
    TH1F* Multiplicity_vs_radius_gen_bkg_ptwt[nCBins][nPtBins];



    TProfile *JetShapeIntegratedParticles_bkgsub[nCBins][nPtBins];
    TProfile *JetShapeDiffParticles_bkgsub[nCBins][nPtBins];
    TProfile *JetShapeIntegratedParticles_bkg[nCBins][nPtBins];
    TProfile *JetShapeDiffParticles_bkg[nCBins][nPtBins];

    TH1F *JetShapeDiffParticles_bkg_1D[nCBins][nPtBins];
    TH1F *JetShapeDiffParticles_1D[nCBins][nPtBins];




    TH1F *JetShapeIntegratedParticles_bkg_1D[nCBins][nPtBins];
    TH1F *JetShapeIntegratedParticles_1D[nCBins][nPtBins];
    TH1F *JetShapeDiffParticlesGen_bkg_1D[nCBins][nPtBins];
    TH1F *JetShapeDiffParticlesGen_1D[nCBins][nPtBins];

    TH1F *JetShapeDiffParticles_1D_check_error[nCBins][nPtBins];


    TH1F* NumNeutral_subleadingjet[nCBins][nPtBins];
    TH1F* NumPhotons_subleadingjet[nCBins][nPtBins];
    TH1F* NumAll_subleadingjet[nCBins][nPtBins];
    TH1F* NumChargedHadrons_subleadingjet[nCBins][nPtBins];



    TH1F* SumJetPtFraction_hist_subleadingjet[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_subleadingjet_photons[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_subleadingjet_neutrals[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_subleadingjet_chargedhadrons[nCBins][nPtBins];

    TH1F* SumJetPtFraction_hist_leadingjet_photons[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_leadingjet_neutrals[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_leadingjet_charged[nCBins][nPtBins];

    TH1F* all_cand_pT_hist_subleadingjet[nCBins][nPtBins];
    TH1F* all_cand_phi_hist_subleadingjet[nCBins][nPtBins];
    TH1F* all_cand_eta_hist_subleadingjet[nCBins][nPtBins];
    TH1F* radius_hist_subleadingjet[nCBins][nPtBins];

    TH1F* Centrality_hist[nCBins][nPtBins];
//    TH1F* Aj[nCBins][nPtBins];
    TH1F* Aj[nCBins];


    TH1F* neutral_cand_pT_hist_subleadingjet[nCBins][nPtBins];
    TH1F* neutral_cand_phi_hist_subleadingjet[nCBins][nPtBins];
    TH1F* neutral_cand_eta_hist_subleadingjet[nCBins][nPtBins];

    TH1F* photons_cand_pT_hist_subleadingjet[nCBins][nPtBins];
    TH1F* photons_cand_phi_hist_subleadingjet[nCBins][nPtBins];
    TH1F* photons_cand_eta_hist_subleadingjet[nCBins][nPtBins];


    TH1F* SumJetPtFraction_hist_leadingjet_chargedhadrons[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_leadingjet_chargedparticles[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_leadingjet_electrons[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_leadingjet_muons[nCBins][nPtBins];

    TH1F* NumChargedHadrons[nCBins][nPtBins];
    TH1F* NumChargedParticles[nCBins][nPtBins];
    TH1F* NumElectrons[nCBins][nPtBins];
    TH1F* NumMuons[nCBins][nPtBins];

    TH1F* chargedhadrons_cand_pT_hist[nCBins][nPtBins];
    TH1F* chargedhadrons_cand_phi_hist[nCBins][nPtBins];
    TH1F* chargedhadrons_cand_eta_hist[nCBins][nPtBins];

    TH1F* chargedparticles_cand_pT_hist[nCBins][nPtBins];
    TH1F* chargedparticles_cand_phi_hist[nCBins][nPtBins];
    TH1F* chargedparticles_cand_eta_hist[nCBins][nPtBins];

    TH1F* electrons_cand_pT_hist[nCBins][nPtBins];
    TH1F* electrons_cand_phi_hist[nCBins][nPtBins];
    TH1F* electrons_cand_eta_hist[nCBins][nPtBins];

    TH1F* muons_cand_pT_hist[nCBins][nPtBins];
    TH1F* muons_cand_phi_hist[nCBins][nPtBins];
    TH1F* muons_cand_eta_hist[nCBins][nPtBins];



    TH1F* chargedhadrons_cand_pT_hist_bkg[nCBins][nPtBins];
    TH1F* chargedparticles_cand_pT_hist_bkg[nCBins][nPtBins];
    TH1F* electrons_cand_pT_hist_bkg[nCBins][nPtBins];
    TH1F* muons_cand_pT_hist_bkg[nCBins][nPtBins];
    TH1F* neutral_cand_pT_hist_bkg[nCBins][nPtBins];
    TH1F* photons_cand_pT_hist_bkg[nCBins][nPtBins];



    TH1F* dN_tracks[nCBins][nPtBins];
    TH1F* dN_chargedhadrons[nCBins][nPtBins];
    TH1F* dN_chargedparticles[nCBins][nPtBins];
    TH1F* dN_electrons[nCBins][nPtBins];
    TH1F* dN_muons[nCBins][nPtBins];
    TH1F* dN_neutrals[nCBins][nPtBins];
    TH1F* dN_photons[nCBins][nPtBins];

    TH1F* SumJetPtFraction_hist_leadingjet_photons_bkg[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_leadingjet_neutrals_bkg[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_leadingjet_chargedhadrons_bkg[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_leadingjet_chargedparticles_bkg[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_leadingjet_electrons_bkg[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_leadingjet_muons_bkg[nCBins][nPtBins];



    /// subleading jet
    TH1F* chargedhadrons_cand_pT_hist_subleadingjet[nCBins][nPtBins];
    TH1F* chargedhadrons_cand_phi_hist_subleadingjet[nCBins][nPtBins];
    TH1F* chargedhadrons_cand_eta_hist_subleadingjet[nCBins][nPtBins];

    TH1F* chargedhadrons_cand_pT_hist_subleadingjet_bkg[nCBins][nPtBins];
    TH1F* neutral_cand_pT_hist_subleadingjet_bkg[nCBins][nPtBins];
    TH1F* photons_cand_pT_hist_subleadingjet_bkg[nCBins][nPtBins];

    TH1F* SumJetPtFraction_hist_subleadingjet_photons_bkg[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_subleadingjet_neutrals_bkg[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_subleadingjet_chargedhadrons_bkg[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_subleadingjet_chargedparticles_bkg[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_subleadingjet_electrons_bkg[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_subleadingjet_muons_bkg[nCBins][nPtBins];

    TH1F* NumAll_subleadingjet_bkg[nCBins][nPtBins];
    TH1F* NumNeutral_subleadingjet_bkg[nCBins][nPtBins];
    TH1F* NumPhotons_subleadingjet_bkg[nCBins][nPtBins];
    TH1F* NumChargedHadrons_subleadingjet_bkg[nCBins][nPtBins];

    TH1F* JetPt_fraction_hist[nCBins][nPtBins];
    TH1F* dPhi_jet_track[nCBins][nPtBins];
    TH1F* dPhi_jet_track_ptweight[nCBins][nPtBins];


    TH1F* JetShapeDiffParticles_bkg_1D_total[nCBins];
    TH1F* JetShapeIntegratedParticles_bkg_1D_total[nCBins];

    TH1D* track_cand_pT_hist_weighted[nCBins][nPtBins];
    TH1D* track_cand_eta_hist_weighted[nCBins][nPtBins];
    TH1D* track_cand_phi_hist_weighted[nCBins][nPtBins];



    TH1F* SumPt_only[nCBins][nPtBins];
    
    TH1D* track_gen_pT_pythia[nCBins][nPtBins];
    TH1D* track_gen_eta_pythia[nCBins][nPtBins];
    TH1D* track_gen_phi_pythia[nCBins][nPtBins];

    TH1D* track_gen_pT_hydjet[nCBins][nPtBins];
    TH1D* track_gen_eta_hydjet[nCBins][nPtBins];
    TH1D* track_gen_phi_hydjet[nCBins][nPtBins];

    TH1D* track_gen_pT_hydjet2[nCBins][nPtBins];
    TH1D* track_gen_pT_hydjet3[nCBins][nPtBins];
    TH1D* track_gen_pT_hydjet4[nCBins][nPtBins];
    TH1D* track_gen_pT_hydjet5[nCBins][nPtBins];
    TH1D* track_gen_pT_hydjet6[nCBins][nPtBins];
    TH1D* track_gen_pT_hydjet7[nCBins][nPtBins];
    TH1D* track_gen_pT_hydjet8[nCBins][nPtBins];
    TH1D* track_gen_pT_hydjet9[nCBins][nPtBins];
    TH1D* track_gen_pT_hydjet10[nCBins][nPtBins];
    TH1D* track_gen_pT_hydjet11[nCBins][nPtBins];

    
    TH1D* track_bkg_pT_hist_weighted[nCBins][nPtBins];
    TH1D* track_bkg_pT_hist[nCBins][nPtBins];
    TH1D* track_gen_bkg_pT_hist[nCBins][nPtBins];

    
    
    TH1F* track_bkg_pT_hist_weighted_oocone[nCBins][nPtBins];
    TH1F* track_bkg_pT_hist_oocone[nCBins][nPtBins];
    TH1F* track_gen_bkg_pT_hist_oocone[nCBins][nPtBins];

    TH1D* Ntracks[nCBins][nPtBins][nTrkPtBins];
    TH1D* Ngen[nCBins][nPtBins][nTrkPtBins];
    TH1D* Ntracks_awayside[nCBins][nPtBins][nTrkPtBins];

    

    TH1D* track_gen_pT_hist[nCBins][nPtBins];
    TH1D* track_gen_eta_hist[nCBins][nPtBins];
    TH1D* track_gen_phi_hist[nCBins][nPtBins];

  /*  TH1F* track_gen_pT_hist[nCBins][nPtBins];
    TH1F* track_gen_eta_hist[nCBins][nPtBins];
    TH1F* track_gen_phi_hist[nCBins][nPtBins];
*/

    TH1F* JetEnergy_resolution[nCBins][nPtBins];
    TH1F* JetEnergy_ratio[nCBins][nPtBins];

    TH1F* JetShapeDiffParticles_subleadingjet_bkg_1DEtaRef[nCBins][nPtBins_sub];
    TH1F* JetShapeIntegratedParticles_subleadingjet_1D[nCBins][nPtBins_sub];
    TH1F* JetShapeDiffParticles_subleadingjet_bkg_1D[nCBins][nPtBins_sub];
    TH1F* JetShapeDiffParticles_subleadingjet_1D[nCBins][nPtBins_sub];

    double sig_count[nCBins][nPtBins];
    double bkg_count[nCBins][nPtBins];
    double sig_count_sub[nCBins][nPtBins_sub];
    double bkg_count_sub[nCBins][nPtBins_sub];

    TH1F* JetEnergy_ratio_new[nCBins][nPtBins];
    TH2F* JetEnergy_gen_vs_rec[nCBins];
    TH1F* fullrecjet_aftercuts[nPtBins];
    TH1F* fullrefjet_aftercuts[nPtBins];

    TH2D* hJetTrackSignalBackground[nCBins][nPtBins][nTrkPtBins];
    TH2D* hJetTrackBackground[nCBins][nPtBins][nTrkPtBins];
    TH1F *JetShapeDiffParticles_1D_EMix[nCBins][nPtBins];
    TH1F* emix_jetpt[nCBins][nPtBins];
    TH1F* emixed_jetpt[nCBins][nPtBins];
    TH1F* emixed_jeteta[nCBins][nPtBins];
    TH1F* emixed_jetphi[nCBins][nPtBins];


    TH2D* hJetGenTrackSignalBackground[nCBins][nPtBins][nTrkPtBins];
    TH2D* hJetGenTrackBackground[nCBins][nPtBins][nTrkPtBins];

    TH2D* hJetTrackSignalBackground_notrkcorr[nCBins][nPtBins][nTrkPtBins];
    TH2D* hJetTrackBackground_notrkcorr[nCBins][nPtBins][nTrkPtBins];

    TH2D* hJetGenTrackSignalBackground_pythia[nCBins][nPtBins][nTrkPtBins];
    TH2D* hJetGenTrackSignalBackground_hydjet[nCBins][nPtBins][nTrkPtBins];

    TH2D* hJetGenTrackSignalBackground_hydro[nCBins][nPtBins][nTrkPtBins];
    TH2D* hJetGenTrackSignalBackground_minijets[nCBins][nPtBins][nTrkPtBins];

    TH2D* hJetGenTrackBackground_pythia[nCBins][nPtBins][nTrkPtBins];
    TH2D* hJetGenTrackBackground_hydjet[nCBins][nPtBins][nTrkPtBins];

    TH1F* TrkPhi[nCBins][nPtBins][nTrkPtBins];
    TH1F* TrkPt[nCBins][nPtBins][nTrkPtBins];
    TH1F* TrkEta[nCBins][nPtBins][nTrkPtBins];

    TH1F* TrkPhi_weighted[nCBins][nPtBins][nTrkPtBins];
    TH1F* TrkPt_weighted[nCBins][nPtBins][nTrkPtBins];
    TH1F* TrkEta_weighted[nCBins][nPtBins][nTrkPtBins];


    TH1F* GenTrkEta[nCBins][nPtBins][nTrkPtBins];
    TH1F* GenTrkPhi[nCBins][nPtBins][nTrkPtBins];
    TH1F* GenTrkPt[nCBins][nPtBins][nTrkPtBins];


    TH1F* GenTrkEta_bkg[nCBins][nPtBins][nTrkPtBins];
    TH1F* GenTrkPhi_bkg[nCBins][nPtBins][nTrkPtBins];
    TH1F* TrkPhi_bkg[nCBins][nPtBins][nTrkPtBins];
    TH1F* TrkEta_bkg[nCBins][nPtBins][nTrkPtBins];


    TH1F* rr_gen_leadingjet[nCBins][nTrkPtBins];
    TH1F* rr_reco_leadingjet[nCBins][nTrkPtBins];
    TH1F* rr_reco_leadingjet_weighted[nCBins][nTrkPtBins];

    TH1F* EtaRef_bkg_pt[nCBins][nPtBins];
    TH1F* EtaRef_bkg_pt_weighted[nCBins][nPtBins];



};

hist_class::hist_class(TString the_desc, bool is_it_data) {
  n_evt_raw = 0;
  desc = the_desc;
  is_data = is_it_data;

  //  JetEnergy_gen_vs_rec = new TH2F((TString) (desc + "_JetEnergy_gen_vs_rec"), "", 500, 0, 1000, 500, 0, 1000);     JetEnergy_gen_vs_rec->Sumw2();
        alljets_rec = new TH1F((TString) (desc + "_alljets_rec"), "", 10, x_Bins_mc);     alljets_rec->Sumw2();

  NumberOfMatches = new TH1F((TString) (desc + "_NumberOfMatches"), "", 100, 0., 100.);     NumberOfMatches->Sumw2();

  NEvents = new TH1F((TString) (desc + "_Nevents"), "This counts the number of events", 100, 0., 100.);   
  NEvents->Sumw2(); 
  NEvents_test = new TH1F((TString) (desc + "_Nevents_test"), "", 100, 0., 100.);  
  NEvents_test->Sumw2();
  NEvents_after_noise = new TH1F((TString) (desc + "_Nevents_after_noise"), "", 100, 0., 100.);     NEvents_after_noise->Sumw2();
  NEvents_after_spike = new TH1F((TString) (desc + "_Nevents_after_spike"), "", 100, 0., 100.);     NEvents_after_spike->Sumw2();

  NEvents_dijets = new TH1F((TString) (desc + "_Nevents_dijets"), "", 100, 0., 100.);     NEvents_dijets->Sumw2();
  NEvents_after_dphi = new TH1F((TString) (desc + "_Nevents_after_dphi"), "", 100, 0., 100.);     NEvents_after_dphi->Sumw2();
  NEvents_before_dphi = new TH1F((TString) (desc + "_Nevents_before_dphi"), "", 100, 0., 100.);     NEvents_before_dphi->Sumw2();
  NEvents_after_trigger = new TH1F((TString) (desc + "_Nevents_after_trigger"), "", 100, 0., 100.);     NEvents_after_trigger->Sumw2();


  AllJetPt_raw_hist = new TH1F((TString) (desc + "_AllJetPt_raw_hist"), "", 100, 0., 500.);     AllJetPt_raw_hist->Sumw2();
  AllJetPt_hist = new TH1F((TString) (desc + "_AllJetPt_hist"), "", 100, 0., 500.);     AllJetPt_hist->Sumw2();
  AllJetPhi_hist = new TH1F((TString) (desc + "_AllJetPhi_hist"), "", 36, 0.,TMath::Pi());     AllJetPhi_hist->Sumw2();
  AllJetEta_hist = new TH1F((TString) (desc + "_AllJetEta_hist"), "", 50, -5., 5.);     AllJetEta_hist->Sumw2();

  First_AllJetPt_hist = new TH1F((TString) (desc + "_First_AllJetPt_hist"), "", 100, 0., 500.);     First_AllJetPt_hist->Sumw2();
  First_AllJetPhi_hist = new TH1F((TString) (desc + "_First_AllJetPhi_hist"), "", 36, 0.,TMath::Pi());     First_AllJetPhi_hist->Sumw2();
  First_AllJetEta_hist = new TH1F((TString) (desc + "_First_AllJetEta_hist"), "", 50, -5., 5.);     First_AllJetEta_hist->Sumw2();

  Sub_AllJetPt_hist = new TH1F((TString) (desc + "_Sub_AllJetPt_hist"), "", 100, 0., 500.);     Sub_AllJetPt_hist->Sumw2();
  Sub_AllJetPhi_hist = new TH1F((TString) (desc + "_Sub_AllJetPhi_hist"), "", 36, 0.,TMath::Pi());     Sub_AllJetPhi_hist->Sumw2();
  Sub_AllJetEta_hist = new TH1F((TString) (desc + "_Sub_AllJetEta_hist"), "", 50, -5., 5.);     Sub_AllJetEta_hist->Sumw2();
  JetPt_fraction = new TH1F((TString) (desc + "_JetPt_fraction"), "", 100,0.,5);     JetPt_fraction->Sumw2();

  Centrality = new TH1F((TString) (desc + "_Centrality"), "", 200,0.,200);     Centrality->Sumw2();
  Centrality_weighted = new TH1F((TString) (desc + "_Centrality_weighted"), "", 200,0.,200);     Centrality_weighted->Sumw2();

  track_vz = new TH1F((TString) (desc + "_track_vz"), "", 80, -20., 20.);
  track_vz_weighted = new TH1F((TString) (desc + "_track_vz_weighted"), "", 80, -20., 20.);


  for (int ibin2=0;ibin2<nPtBins;ibin2++)
  {
    fullrecjet_aftercuts[ibin2] = new TH1F((TString) (desc + "_fullrecjet_aftercuts" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 10, x_Bins_mc);    fullrecjet_aftercuts[ibin2]->Sumw2();
    fullrefjet_aftercuts[ibin2] = new TH1F((TString) (desc + "_fullrefjet_aftercuts" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 10, x_Bins_mc);    fullrefjet_aftercuts[ibin2]->Sumw2();

  }


  for (int ibin=0;ibin<nCBins;ibin++)
  {


    jet_pT_hist[ibin] = new TH1F((TString) (desc + "_jet_pT_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 25, 0., 500.); jet_pT_hist[ibin]->Sumw2();
    jet_phi_hist[ibin] = new TH1F((TString) (desc + "_jet_phi_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] ), "", 36, 0.,TMath::Pi());
    jet_phi_hist[ibin]->Sumw2();
    jet_eta_hist[ibin] = new TH1F((TString) (desc + "_jet_eta_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] ), "", 50, -5., 5.);     jet_eta_hist[ibin]->Sumw2();
    jet_corrpT_hist[ibin] = new TH1F((TString) (desc + "_jet_corrpT_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 25, 0., 500.);     jet_corrpT_hist[ibin]->Sumw2();
    LeadingJetPt_hist[ibin] = new TH1F((TString) (desc + "_LeadingJetPt_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] ), "", 25, 0., 500.);     LeadingJetPt_hist[ibin]->Sumw2();
    LeadingJetPhi_hist[ibin] = new TH1F((TString) (desc + "_LeadingJetPhi_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 72, -TMath::Pi(),TMath::Pi());     LeadingJetPhi_hist[ibin]->Sumw2();
    LeadingJetEta_hist[ibin] = new TH1F((TString) (desc + "_LeadingJetEta_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 50, -5., 5.);     LeadingJetEta_hist[ibin]->Sumw2();
    SubJetPt_hist[ibin] = new TH1F((TString) (desc + "_SubJetPt_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 25, 0., 500.);     SubJetPt_hist[ibin]->Sumw2();
    SubJetPhi_hist[ibin] = new TH1F((TString) (desc + "_SubJetPhi_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 72, -TMath::Pi(),TMath::Pi());     SubJetPhi_hist[ibin]->Sumw2();
    SubJetEta_hist[ibin] = new TH1F((TString) (desc + "_SubJetEta_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 50, -5., 5.);     SubJetEta_hist[ibin]->Sumw2();
    dPhi_leadingjet_hist[ibin] = new TH1F((TString) (desc + "_dPhi_leadingjet_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 30,0,3.14159);     dPhi_leadingjet_hist[ibin]->Sumw2();
    dPhi_subleadingjet_hist[ibin] = new TH1F((TString) (desc + "_dPhi_subleadingjet_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 30,0,3.14159);     dPhi_subleadingjet_hist[ibin]->Sumw2();
    JetShapeDiffParticles_bkg_1D_total[ibin] = new TH1F((TString) (desc + "_JetShapeDiffParticles_bkg_1D_total"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 6,0.0,0.3);     JetShapeDiffParticles_bkg_1D_total[ibin]->Sumw2();     
    JetShapeIntegratedParticles_bkg_1D_total[ibin] = new TH1F((TString) (desc + "_JetShapeIntegratedParticles_bkg_1D_total"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] ), "", 6,0.05,0.35);     JetShapeIntegratedParticles_bkg_1D_total[ibin]->Sumw2();
    JetEnergy_gen_vs_rec[ibin] = new TH2F((TString) (desc + "_JetEnergy_gen_vs_rec"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 500, 0, 1000, 500, 0, 1000);     JetEnergy_gen_vs_rec[ibin]->Sumw2();
    dPhi_hist[ibin] = new TH1F((TString) (desc + "_dPhi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 30, 0.,TMath::Pi());     dPhi_hist[ibin]->Sumw2();
    dPhi_after_hist[ibin] = new TH1F((TString) (desc + "_dPhi_after_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] ), "", 30,0,3.14159);     dPhi_after_hist[ibin]->Sumw2();
    dPhi_vector[ibin] = new TH1D((TString) (desc + "_dPhi_vector"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] ), "", 60,-2*TMath::Pi(),2*TMath::Pi());     dPhi_vector[ibin]->Sumw2();
    
    for (int ibin3=0;ibin3<nTrkPtBins;ibin3++)
    {
  	rr_gen_leadingjet[ibin][ibin3] = new TH1F((TString) (desc + "_rr_gen_leadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 50, 0., 5.);    rr_gen_leadingjet[ibin][ibin3]->Sumw2();
      	rr_reco_leadingjet[ibin][ibin3] = new TH1F((TString) (desc + "_rr_reco_leadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 50, 0., 5.);    rr_reco_leadingjet[ibin][ibin3]->Sumw2();
      	rr_reco_leadingjet_weighted[ibin][ibin3] = new TH1F((TString) (desc + "_rr_reco_leadingjet_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 50, 0., 5.);    rr_reco_leadingjet_weighted[ibin][ibin3]->Sumw2();


   }

      Aj[ibin] = new TH1F((TString) (desc + "_Aj"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 25,0.,1.25);     Aj[ibin]->Sumw2();


    for (int ibin2=0;ibin2<nPtBins;ibin2++){ 
      
      
      ThirdJetPt_hist[ibin][ibin2] = new TH1F((TString) (desc + "_ThirdJetPt_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     ThirdJetPt_hist[ibin][ibin2]->Sumw2();
      ThirdJetPhi_hist[ibin][ibin2] = new TH1F((TString) (desc + "_ThirdJetPhi_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     ThirdJetPhi_hist[ibin][ibin2]->Sumw2();
      ThirdJetEta_hist[ibin][ibin2] = new TH1F((TString) (desc + "_ThirdJetEta_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     ThirdJetEta_hist[ibin][ibin2]->Sumw2();

      //dPhi_hist[ibin][ibin2] = new TH1F((TString) (desc + "_dPhi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 30, 0.,TMath::Pi());     dPhi_hist[ibin][ibin2]->Sumw2();
      dEta_hist[ibin][ibin2] = new TH1F((TString) (desc + "_dEta_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100,-5,5);     dEta_hist[ibin][ibin2]->Sumw2();
     /// dPhi_after_hist[ibin][ibin2] = new TH1F((TString) (desc + "_dPhi_after_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 30,0,3.14159);     dPhi_after_hist[ibin][ibin2]->Sumw2();
      dPhi_jet_track[ibin][ibin2] = new TH1F((TString) (desc + "_dPhi_jet_track"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 30,0,3.14159);     dPhi_jet_track[ibin][ibin2]->Sumw2();
      dPhi_jet_track_ptweight[ibin][ibin2] = new TH1F((TString) (desc + "_dPhi_jet_track_ptweight"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 30,0,3.14159);     dPhi_jet_track_ptweight[ibin][ibin2]->Sumw2();
      // dPhi_hist[ibin][ibin2] = new TH1F((TString) (desc + "_dPhi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 6);     // dPhi_hist[ibin][ibin2]->Sumw2();
      //dPhi_after_hist[ibin][ibin2] = new TH1F((TString) (desc + "_dPhi_after_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 6);     //dPhi_after_hist[ibin][ibin2]->Sumw2();
      all_jets_corrpT[ibin][ibin2] = new TH1F((TString) (desc + "_all_jets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, 0., 500.);     all_jets_corrpT[ibin][ibin2]->Sumw2();
      all_jets_phi[ibin][ibin2] = new TH1F((TString) (desc + "_all_jets_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 72, -TMath::Pi(), TMath::Pi());     all_jets_phi[ibin][ibin2]->Sumw2();
      all_jets_eta[ibin][ibin2] = new TH1F((TString) (desc + "_all_jets_eta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     all_jets_eta[ibin][ibin2]->Sumw2();
      //// leading jet histograms
      only_leadingjets_corrpT[ibin][ibin2] = new TH1F((TString) (desc + "_only_leadingjets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     only_leadingjets_corrpT[ibin][ibin2]->Sumw2();
      only_leadingjets_phi[ibin][ibin2] = new TH1F((TString) (desc + "_only_leadingjets_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     only_leadingjets_phi[ibin][ibin2]->Sumw2();
      only_leadingjets_eta[ibin][ibin2] = new TH1F((TString) (desc + "_only_leadingjets_eta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     only_leadingjets_eta[ibin][ibin2]->Sumw2();
      all_cand_pT_hist[ibin][ibin2] = new TH1F((TString) (desc + "_all_cand_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     all_cand_pT_hist[ibin][ibin2]->Sumw2();
      all_cand_phi_hist[ibin][ibin2] = new TH1F((TString) (desc + "_all_cand_phi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     all_cand_phi_hist[ibin][ibin2]->Sumw2();
      all_cand_eta_hist[ibin][ibin2] = new TH1F((TString) (desc + "_all_cand_eta_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     all_cand_eta_hist[ibin][ibin2]->Sumw2();
      track_cand_pT_hist[ibin][ibin2] = new TH1D((TString) (desc + "_track_cand_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);     track_cand_pT_hist[ibin][ibin2]->Sumw2();
      track_cand_phi_hist[ibin][ibin2] = new TH1D((TString) (desc + "_track_cand_phi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 72, -TMath::Pi(),TMath::Pi());     track_cand_phi_hist[ibin][ibin2]->Sumw2();
      track_cand_eta_hist[ibin][ibin2] = new TH1D((TString) (desc + "_track_cand_eta_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 12, -2.4, 2.4);     track_cand_eta_hist[ibin][ibin2]->Sumw2();
      track_cand_pT_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_track_cand_pT_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     track_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();


      NumNeutral[ibin][ibin2] = new TH1F((TString) (desc + "_NumNeutral"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumNeutral[ibin][ibin2]->Sumw2();
      NumPhotons[ibin][ibin2] = new TH1F((TString) (desc + "_NumPhotons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumPhotons[ibin][ibin2]->Sumw2();
      NumAll[ibin][ibin2] = new TH1F((TString) (desc + "_NumAll"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumAll[ibin][ibin2]->Sumw2();
      NumCharged[ibin][ibin2] = new TH1F((TString) (desc + "_NumCharged"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumCharged[ibin][ibin2]->Sumw2();
      NumChargedHadrons[ibin][ibin2] = new TH1F((TString) (desc + "_NumChargedHadrons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumChargedHadrons[ibin][ibin2]->Sumw2();
      NumChargedParticles[ibin][ibin2] = new TH1F((TString) (desc + "_NumChargedParticles"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumChargedParticles[ibin][ibin2]->Sumw2();
      NumElectrons[ibin][ibin2] = new TH1F((TString) (desc + "_NumElectrons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumElectrons[ibin][ibin2]->Sumw2();
      NumMuons[ibin][ibin2] = new TH1F((TString) (desc + "_NumMuons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumMuons[ibin][ibin2]->Sumw2();

      NumNeutral_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_NumNeutral_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumNeutral_bkg[ibin][ibin2]->Sumw2();
      NumPhotons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_NumPhotons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumPhotons_bkg[ibin][ibin2]->Sumw2();
      NumAll_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_NumAll_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumAll_bkg[ibin][ibin2]->Sumw2();
      NumChargedHadrons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_NumChargedHadrons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumChargedHadrons_bkg[ibin][ibin2]->Sumw2();
      NumChargedParticles_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_NumChargedParticles_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumChargedParticles_bkg[ibin][ibin2]->Sumw2();
      NumElectrons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_NumElectrons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumElectrons_bkg[ibin][ibin2]->Sumw2();
      NumMuons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_NumMuons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumMuons_bkg[ibin][ibin2]->Sumw2();
      NumNeutral_subleadingjet_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_NumNeutral_subleadingjet_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumNeutral_subleadingjet_bkg[ibin][ibin2]->Sumw2();
      NumPhotons_subleadingjet_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_NumPhotons_subleadingjet_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumPhotons_subleadingjet_bkg[ibin][ibin2]->Sumw2();
      NumAll_subleadingjet_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_NumAll_subleadingjet_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumAll_subleadingjet_bkg[ibin][ibin2]->Sumw2();
      NumChargedHadrons_subleadingjet_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_NumChargedHadrons_subleadingjet_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumChargedHadrons_subleadingjet_bkg[ibin][ibin2]->Sumw2();

      ///===================////////======================//////////////===================//////////=================//////==========

      chargedhadrons_cand_pT_hist[ibin][ibin2] = new TH1F((TString) (desc + "_chargedhadrons_cand_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     chargedhadrons_cand_pT_hist[ibin][ibin2]->Sumw2();
      chargedhadrons_cand_phi_hist[ibin][ibin2] = new TH1F((TString) (desc + "_chargedhadrons_cand_phi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     chargedhadrons_cand_phi_hist[ibin][ibin2]->Sumw2();
      chargedhadrons_cand_eta_hist[ibin][ibin2] = new TH1F((TString) (desc + "_chargedhadrons_cand_eta_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     chargedhadrons_cand_eta_hist[ibin][ibin2]->Sumw2();

      chargedparticles_cand_pT_hist[ibin][ibin2] = new TH1F((TString) (desc + "_chargedparticles_cand_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     chargedparticles_cand_pT_hist[ibin][ibin2]->Sumw2();
      chargedparticles_cand_phi_hist[ibin][ibin2] = new TH1F((TString) (desc + "_chargedparticles_cand_phi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     chargedparticles_cand_phi_hist[ibin][ibin2]->Sumw2();
      chargedparticles_cand_eta_hist[ibin][ibin2] = new TH1F((TString) (desc + "_chargedparticles_cand_eta_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     chargedparticles_cand_eta_hist[ibin][ibin2]->Sumw2();

      electrons_cand_pT_hist[ibin][ibin2] = new TH1F((TString) (desc + "_electrons_cand_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     electrons_cand_pT_hist[ibin][ibin2]->Sumw2();
      electrons_cand_phi_hist[ibin][ibin2] = new TH1F((TString) (desc + "_electrons_cand_phi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     electrons_cand_phi_hist[ibin][ibin2]->Sumw2();
      electrons_cand_eta_hist[ibin][ibin2] = new TH1F((TString) (desc + "_electrons_cand_eta_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     electrons_cand_eta_hist[ibin][ibin2]->Sumw2();

      muons_cand_pT_hist[ibin][ibin2] = new TH1F((TString) (desc + "_muons_cand_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     muons_cand_pT_hist[ibin][ibin2]->Sumw2();
      muons_cand_phi_hist[ibin][ibin2] = new TH1F((TString) (desc + "_muons_cand_phi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     muons_cand_phi_hist[ibin][ibin2]->Sumw2();
      muons_cand_eta_hist[ibin][ibin2] = new TH1F((TString) (desc + "_muons_cand_eta_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     muons_cand_eta_hist[ibin][ibin2]->Sumw2();

      neutral_cand_pT_hist[ibin][ibin2] = new TH1F((TString) (desc + "_neutral_cand_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     neutral_cand_pT_hist[ibin][ibin2]->Sumw2();
      neutral_cand_phi_hist[ibin][ibin2] = new TH1F((TString) (desc + "_neutral_cand_phi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     neutral_cand_phi_hist[ibin][ibin2]->Sumw2();
      neutral_cand_eta_hist[ibin][ibin2] = new TH1F((TString) (desc + "_neutral_cand_eta_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     neutral_cand_eta_hist[ibin][ibin2]->Sumw2();

      photons_cand_pT_hist[ibin][ibin2] = new TH1F((TString) (desc + "_photons_cand_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     photons_cand_pT_hist[ibin][ibin2]->Sumw2();
      photons_cand_phi_hist[ibin][ibin2] = new TH1F((TString) (desc + "_photons_cand_phi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     photons_cand_phi_hist[ibin][ibin2]->Sumw2();
      photons_cand_eta_hist[ibin][ibin2] = new TH1F((TString) (desc + "_photons_cand_eta_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     photons_cand_eta_hist[ibin][ibin2]->Sumw2();

      ///===================////////======================//////////////===================//////////=================//////==========
      ///===================////////======================//////////////===================//////////=================//////==========
      ///===================////////======================//////////////===================//////////=================//////==========



      chargedhadrons_cand_pT_hist_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_chargedhadrons_cand_pT_hist_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     chargedhadrons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();

      chargedparticles_cand_pT_hist_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_chargedparticles_cand_pT_hist_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     chargedparticles_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();

      electrons_cand_pT_hist_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_electrons_cand_pT_hist_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     electrons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();

      muons_cand_pT_hist_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_muons_cand_pT_hist_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     muons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();

      neutral_cand_pT_hist_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_neutral_cand_pT_hist_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     neutral_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();

      photons_cand_pT_hist_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_photons_cand_pT_hist_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     photons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();


      /// subleading jet

      chargedhadrons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_chargedhadrons_cand_pT_hist_subleadingjet_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     chargedhadrons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Sumw2();
      neutral_cand_pT_hist_subleadingjet_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_neutral_cand_pT_hist_subleadingjet_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     neutral_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Sumw2();
      photons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_photons_cand_pT_hist_subleadingjet_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     photons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Sumw2();

      ///===================////////======================//////////////===================//////////=================//////==========



      //SumJetPtFraction_hist[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     //SumJetPtFraction_hist[ibin][ibin2]->Sumw2();


      SumJetPtFraction_hist[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100,0.,500);     SumJetPtFraction_hist[ibin][ibin2]->Sumw2();


      //  JetShapeIntegratedParticles[ibin][ibin2] = new TProfile((TString) (desc + "_JetShapeIntegratedParticles"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.05,0.35,0.,2.);     JetShapeIntegratedParticles[ibin][ibin2]->Sumw2();



      JetShapeIntegratedParticles[ibin][ibin2] = new TH1F((TString) (desc + "_JetShapeIntegratedParticles"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.05,0.35);     JetShapeIntegratedParticles[ibin][ibin2]->Sumw2();



      //JetShapeDiffParticles[ibin][ibin2] = new TProfile((TString) (desc + "_JetShapeDiffParticles"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3,0.,100.);     JetShapeDiffParticles[ibin][ibin2]->Sumw2();


      JetShapeDiffParticles[ibin][ibin2] = new TH1F((TString) (desc + "_JetShapeDiffParticles"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     JetShapeDiffParticles[ibin][ibin2]->Sumw2();




      JetShapeIntegratedParticles_bkgsub[ibin][ibin2] = new TProfile((TString) (desc + "_JetShapeIntegratedParticles_bkgsub"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.05,0.35,0.,2.);     JetShapeIntegratedParticles_bkgsub[ibin][ibin2]->Sumw2();
      JetShapeDiffParticles_bkgsub[ibin][ibin2] = new TProfile((TString) (desc + "_JetShapeDiffParticles_bkgsub"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3,0.,100.);     JetShapeDiffParticles_bkgsub[ibin][ibin2]->Sumw2();

      JetShapeIntegratedParticles_bkg[ibin][ibin2] = new TProfile((TString) (desc + "_JetShapeIntegratedParticles_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.05,0.35,0.,2.);     JetShapeIntegratedParticles_bkg[ibin][ibin2]->Sumw2();
      JetShapeDiffParticles_bkg[ibin][ibin2] = new TProfile((TString) (desc + "_JetShapeDiffParticles_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3,0.,100.);     JetShapeDiffParticles_bkg[ibin][ibin2]->Sumw2();


      radius_hist[ibin][ibin2] = new TH1F((TString) (desc + "_radius_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, 0., 5.);     radius_hist[ibin][ibin2]->Sumw2();


      radius_hist_mine[ibin][ibin2] = new TH1F((TString) (desc + "_radius_hist_mine"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     radius_hist_mine[ibin][ibin2]->Sumw2();

      radius_hist_mine2[ibin][ibin2] = new TH1F((TString) (desc + "_radius_hist_mine2"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     radius_hist_mine2[ibin][ibin2]->Sumw2();


      radius_hist_mine3[ibin][ibin2] = new TH1F((TString) (desc + "_radius_hist_mine3"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     radius_hist_mine3[ibin][ibin2]->Sumw2();




      /*   J(tShapeDiffParticles_bkg_1D[ibin][ibin2] = new TH1F((TString) (desc + "_JetShapeDiffParticles_bkg_1D"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     JetShapeDiffParticles_bkg_1D[ibin][ibin2]->Sumw2();     
	   JetShapeDiffParticles_1D[ibin][ibin2] = new TH1F((TString) (desc + "_JetShapeDiffParticles_1D"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     JetShapeDiffParticles_1D[ibin][ibin2]->Sumw2();


	   JetShapeIntegratedParticles_bkg_1D[ibin][ibin2] = new TH1F((TString) (desc + "_JetShapeIntegratedParticles_bkg_1D"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.05,0.35);     JetShapeIntegratedParticles_bkg_1D[ibin][ibin2]->Sumw2();
	   JetShapeIntegratedParticles_1D[ibin][ibin2] = new TH1F((TString) (desc + "_JetShapeIntegratedParticles_1D"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.05,0.35);     JetShapeIntegratedParticles_1D[ibin][ibin2]->Sumw2();

       */

      JetShapeDiffParticles_bkg_1D[ibin][ibin2] = new TH1F((TString) (desc + "_JetShapeDiffParticles_bkg_1D"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     JetShapeDiffParticles_bkg_1D[ibin][ibin2]->Sumw2();     
      JetShapeDiffParticles_1D[ibin][ibin2] =     new TH1F((TString) (desc + "_JetShapeDiffParticles_1D"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     JetShapeDiffParticles_1D[ibin][ibin2]->Sumw2();


      JetShapeIntegratedParticles_bkg_1D[ibin][ibin2] = new TH1F((TString) (desc + "_JetShapeIntegratedParticles_bkg_1D"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.05,0.35);     JetShapeIntegratedParticles_bkg_1D[ibin][ibin2]->Sumw2();
      JetShapeIntegratedParticles_1D[ibin][ibin2] = new TH1F((TString) (desc + "_JetShapeIntegratedParticles_1D"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.05,0.35);     JetShapeIntegratedParticles_1D[ibin][ibin2]->Sumw2();

      JetShapeDiffParticlesGen_1D[ibin][ibin2] = new TH1F((TString) (desc + "_JetShapeDiffParticlesGen_1D"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     JetShapeDiffParticlesGen_1D[ibin][ibin2]->Sumw2();

      JetShapeDiffParticlesGen_bkg_1D[ibin][ibin2] = new TH1F((TString) (desc + "_JetShapeDiffParticlesGen_bkg_1D"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);    
      JetShapeDiffParticlesGen_bkg_1D[ibin][ibin2]->Sumw2();

      //// subleading jet histograms
      only_subleadingjets_corrpT[ibin][ibin2] = new TH1F((TString) (desc + "_only_subleadingjets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     only_subleadingjets_corrpT[ibin][ibin2]->Sumw2();
      only_subleadingjets_phi[ibin][ibin2] = new TH1F((TString) (desc + "_only_subleadingjets_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     only_subleadingjets_phi[ibin][ibin2]->Sumw2();
      only_subleadingjets_eta[ibin][ibin2] = new TH1F((TString) (desc + "_only_subleadingjets_eta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     only_subleadingjets_eta[ibin][ibin2]->Sumw2();
      all_cand_pT_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_all_cand_pT_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     all_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();
      all_cand_phi_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_all_cand_phi_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     all_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();
      all_cand_eta_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_all_cand_eta_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     all_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();
      neutral_cand_pT_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_neutral_cand_pT_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     neutral_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();
      neutral_cand_phi_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_neutral_cand_phi_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     neutral_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();
      neutral_cand_eta_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_neutral_cand_eta_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     neutral_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();
      photons_cand_pT_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_photons_cand_pT_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     photons_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();
      photons_cand_phi_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_photons_cand_phi_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     photons_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();
      photons_cand_eta_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_photons_cand_eta_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     photons_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();

      chargedhadrons_cand_pT_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_chargedhadrons_cand_pT_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     chargedhadrons_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();
      chargedhadrons_cand_phi_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_chargedhadrons_cand_phi_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     chargedhadrons_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();
      chargedhadrons_cand_eta_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_chargedhadrons_cand_eta_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     chargedhadrons_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();


      NumNeutral_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_NumNeutral_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumNeutral_subleadingjet[ibin][ibin2]->Sumw2();
      NumPhotons_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_NumPhotons_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumPhotons_subleadingjet[ibin][ibin2]->Sumw2();
      NumAll_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_NumAll_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumAll_subleadingjet[ibin][ibin2]->Sumw2();

      NumChargedHadrons_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_NumChargedHadrons_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumChargedHadrons_subleadingjet[ibin][ibin2]->Sumw2();



      //SumJetPtFraction_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100,0.,300);     //SumJetPtFraction_hist_subleadingjet[ibin][ibin2]->Sumw2();


      SumJetPtFraction_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_subleadingjet[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_photons[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_subleadingjet_photons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_subleadingjet_photons[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_neutrals[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_subleadingjet_neutrals"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_subleadingjet_neutrals[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_chargedhadrons[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_subleadingjet_chargedhadrons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_subleadingjet_chargedhadrons[ibin][ibin2]->Sumw2();

      SumJetPtFraction_hist_leadingjet_photons[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_photons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_photons[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_neutrals[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_neutrals"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_neutrals[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_charged[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_charged"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_charged[ibin][ibin2]->Sumw2();

      SumJetPtFraction_hist_leadingjet_chargedhadrons[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_chargedhadrons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_chargedhadrons[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_chargedparticles[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_chargedparticles"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_chargedparticles[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_electrons[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_electrons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_electrons[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_muons[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_muons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_muons[ibin][ibin2]->Sumw2();



      /// bkg -- sumpt fraction plots 

      SumJetPtFraction_hist_leadingjet_photons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_photons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_photons_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_neutrals_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_neutrals_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_neutrals_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_chargedhadrons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_chargedhadrons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_chargedhadrons_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_chargedparticles_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_chargedparticles_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_chargedparticles_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_electrons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_electrons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_electrons_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_muons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_muons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_muons_bkg[ibin][ibin2]->Sumw2();



      SumJetPtFraction_hist_subleadingjet_photons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_subleadingjet_photons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_subleadingjet_photons_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_neutrals_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_subleadingjet_neutrals_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_subleadingjet_neutrals_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_chargedhadrons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_subleadingjet_chargedhadrons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_subleadingjet_chargedhadrons_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_chargedparticles_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_subleadingjet_chargedparticles_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_subleadingjet_chargedparticles_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_electrons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_subleadingjet_electrons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_subleadingjet_electrons_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_muons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_subleadingjet_muons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_subleadingjet_muons_bkg[ibin][ibin2]->Sumw2();


      radius_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_radius_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, 0., 5.);     radius_hist_subleadingjet[ibin][ibin2]->Sumw2();
      Centrality_hist[ibin][ibin2] = new TH1F((TString) (desc + "_Centrality_hist_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 200, 0., 200.);     Centrality_hist[ibin][ibin2]->Sumw2();
      //Aj[ibin][ibin2] = new TH1F((TString) (desc + "_Aj"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 25,0.,1.25);     Aj[ibin][ibin2]->Sumw2();

      dN_tracks[ibin][ibin2] = new TH1F((TString) (desc + "_dN_tracks"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 300,0.,10.);     dN_tracks[ibin][ibin2]->Sumw2();
      dN_chargedhadrons[ibin][ibin2] = new TH1F((TString) (desc + "_dN_chargedhadrons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 300,0.,10.);     dN_chargedhadrons[ibin][ibin2]->Sumw2();
      dN_chargedparticles[ibin][ibin2] = new TH1F((TString) (desc + "_dN_chargedparticles"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 300,0.,10.);     dN_chargedparticles[ibin][ibin2]->Sumw2();
      dN_electrons[ibin][ibin2] = new TH1F((TString) (desc + "_dN_electrons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 300,0.,10.);     dN_electrons[ibin][ibin2]->Sumw2();
      dN_muons[ibin][ibin2] = new TH1F((TString) (desc + "_dN_muons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 300,0.,10.);     dN_muons[ibin][ibin2]->Sumw2();
      dN_neutrals[ibin][ibin2] = new TH1F((TString) (desc + "_dN_neutrals"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 300,0.,10.);     dN_neutrals[ibin][ibin2]->Sumw2();
      dN_photons[ibin][ibin2] = new TH1F((TString) (desc + "_dN_photons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 300,0.,10.);     dN_photons[ibin][ibin2]->Sumw2();

      JetPt_fraction_hist[ibin][ibin2] = new TH1F((TString) (desc + "_JetPt_fraction_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100,0.,5);     JetPt_fraction_hist[ibin][ibin2]->Sumw2();



      track_cand_pT_hist_weighted[ibin][ibin2] = new TH1D((TString) (desc + "_track_cand_pT_hist_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);     track_cand_pT_hist_weighted[ibin][ibin2]->Sumw2();

        track_cand_eta_hist_weighted[ibin][ibin2] = new TH1D((TString) (desc + "_track_cand_eta_hist_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 12, -2.4, 2.4);     track_cand_eta_hist_weighted[ibin][ibin2]->Sumw2();

	        track_cand_phi_hist_weighted[ibin][ibin2] = new TH1D((TString) (desc + "_track_cand_phi_hist_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 72, -TMath::Pi(), TMath::Pi());     track_cand_phi_hist_weighted[ibin][ibin2]->Sumw2();



///////// gen level eta-phi-pt spectras

          track_gen_pT_hist[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);   track_gen_pT_hist[ibin][ibin2]->Sumw2();

          track_gen_eta_hist[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_eta_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 12, -2.4, 2.4);   track_gen_eta_hist[ibin][ibin2]->Sumw2();

          track_gen_phi_hist[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_phi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 72, -TMath::Pi
              (), TMath::Pi());   track_gen_phi_hist[ibin][ibin2]->Sumw2();

////
          track_gen_pT_pythia[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_pT_pythia"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);   track_gen_pT_pythia[ibin][ibin2]->Sumw2();

          track_gen_eta_pythia[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_eta_pythia"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 12, -2.4, 2.4);   track_gen_eta_pythia[ibin][ibin2]->Sumw2();

          track_gen_phi_pythia[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_phi_pythia"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 72, -TMath::Pi
              (), TMath::Pi());   track_gen_phi_pythia[ibin][ibin2]->Sumw2();

          ////

          track_gen_pT_hydjet[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_pT_hydjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);   track_gen_pT_hydjet[ibin][ibin2]->Sumw2();

          track_gen_eta_hydjet[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_eta_hydjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 12, -2.4, 2.4);   track_gen_eta_hydjet[ibin][ibin2]->Sumw2();

          track_gen_phi_hydjet[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_phi_hydjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 72, -TMath::Pi(), TMath::Pi());   track_gen_phi_hydjet[ibin][ibin2]->Sumw2();

        
        
        
        track_gen_pT_hydjet2[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_pT_hydjet2"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);   track_gen_pT_hydjet2[ibin][ibin2]->Sumw2();

        track_gen_pT_hydjet3[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_pT_hydjet3"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);   track_gen_pT_hydjet3[ibin][ibin2]->Sumw2();

        track_gen_pT_hydjet4[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_pT_hydjet4"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);   track_gen_pT_hydjet4[ibin][ibin2]->Sumw2();

        track_gen_pT_hydjet5[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_pT_hydjet5"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);   track_gen_pT_hydjet5[ibin][ibin2]->Sumw2();

       track_gen_pT_hydjet6[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_pT_hydjet6"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);   track_gen_pT_hydjet6[ibin][ibin2]->Sumw2();

        track_gen_pT_hydjet7[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_pT_hydjet7"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);   track_gen_pT_hydjet7[ibin][ibin2]->Sumw2();

        track_gen_pT_hydjet8[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_pT_hydjet8"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);   track_gen_pT_hydjet8[ibin][ibin2]->Sumw2();

        track_gen_pT_hydjet9[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_pT_hydjet9"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);   track_gen_pT_hydjet9[ibin][ibin2]->Sumw2();

        track_gen_pT_hydjet10[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_pT_hydjet10"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);   track_gen_pT_hydjet10[ibin][ibin2]->Sumw2();

        track_gen_pT_hydjet11[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_pT_hydjet11"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);   track_gen_pT_hydjet11[ibin][ibin2]->Sumw2();

        
//// bkg pt spectras


EtaRef_bkg_pt[ibin][ibin2] = new TH1F((TString) (desc + "_EtaRef_bkg_pt"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk); EtaRef_bkg_pt[ibin][ibin2]->Sumw2();

EtaRef_bkg_pt_weighted[ibin][ibin2] = new TH1F((TString) (desc + "_EtaRef_bkg_pt_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk); EtaRef_bkg_pt_weighted[ibin][ibin2]->Sumw2();
    
		track_bkg_pT_hist_weighted[ibin][ibin2] = new TH1D((TString) (desc + "_track_bkg_pT_hist_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);     track_bkg_pT_hist_weighted[ibin][ibin2]->Sumw2();


		track_bkg_pT_hist[ibin][ibin2] = new TH1D((TString) (desc + "_track_bkg_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);     track_bkg_pT_hist[ibin][ibin2]->Sumw2();

		track_gen_bkg_pT_hist[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_bkg_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);   track_gen_bkg_pT_hist[ibin][ibin2]->Sumw2();


  /*    track_bkg_pT_hist_weighted[ibin][ibin2] = new TH1F((TString) (desc + "_track_bkg_pT_hist_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);     track_bkg_pT_hist_weighted[ibin][ibin2]->Sumw2();


      track_bkg_pT_hist[ibin][ibin2] = new TH1F((TString) (desc + "_track_bkg_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);     track_bkg_pT_hist[ibin][ibin2]->Sumw2();

                track_gen_bkg_pT_hist[ibin][ibin2] = new TH1F((TString) (desc + "_track_gen_bkg_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);   track_gen_bkg_pT_hist[ibin][ibin2]->Sumw2();

*/

      track_bkg_pT_hist_weighted_oocone[ibin][ibin2] = new TH1F((TString) (desc + "_track_bkg_pT_hist_weighted_oocone"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);     track_bkg_pT_hist_weighted_oocone[ibin][ibin2]->Sumw2();

            track_bkg_pT_hist_oocone[ibin][ibin2] = new TH1F((TString) (desc + "_track_bkg_pT_hist_oocone"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);     track_bkg_pT_hist_oocone[ibin][ibin2]->Sumw2();

	                    track_gen_bkg_pT_hist_oocone[ibin][ibin2] = new TH1F((TString) (desc + "_track_gen_bkg_pT_hist_oocone"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);   track_gen_bkg_pT_hist_oocone[ibin][ibin2]->Sumw2();




			    genparticles_pT_hist[ibin][ibin2] = new TH1F((TString) (desc + "_genparticles_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);     genparticles_pT_hist[ibin][ibin2]->Sumw2();


			    genparticles_bkg_pT_hist[ibin][ibin2] = new TH1F((TString) (desc + "_genparticles_bkg_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);     genparticles_bkg_pT_hist[ibin][ibin2]->Sumw2();


  gensube_hist[ibin][ibin2] = new TH1F((TString) (desc + "_gensube_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0, 100);  gensube_hist[ibin][ibin2]->Sumw2();



			    ///// bkg eta distributions

			    track_bkg_eta_hist_weighted[ibin][ibin2] = new TH1D((TString) (desc + "_track_bkg_eta_hist_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 12, -2.4, 2.4);     track_bkg_eta_hist_weighted[ibin][ibin2]->Sumw2();


			    track_bkg_eta_hist[ibin][ibin2] = new TH1D((TString) (desc + "_track_bkg_eta_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 12, -2.4, 2.4);     track_bkg_eta_hist[ibin][ibin2]->Sumw2();

			    track_gen_bkg_eta_hist[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_bkg_eta_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 12, -2.4, 2.4);   track_gen_bkg_eta_hist[ibin][ibin2]->Sumw2();




      JetShapeDiffParticles_1D_check_error[ibin][ibin2] = new TH1F((TString) (desc + "_JetShapeDiffParticles_1D_check_error"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, 0., 5);  JetShapeDiffParticles_1D_check_error[ibin][ibin2]->Sumw2();





      Multiplicity_vs_radius[ibin][ibin2] = new TH1F((TString) (desc + "_Multiplicity_vs_radius"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     Multiplicity_vs_radius[ibin][ibin2]->Sumw2();

      Multiplicity_vs_radius_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_Multiplicity_vs_radius_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     Multiplicity_vs_radius_bkg[ibin][ibin2]->Sumw2();


      Multiplicity_vs_radius_gen[ibin][ibin2] = new TH1F((TString) (desc + "_Multiplicity_vs_radius_gen"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     Multiplicity_vs_radius_gen[ibin][ibin2]->Sumw2();

      Multiplicity_vs_radius_gen_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_Multiplicity_vs_radius_gen_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     Multiplicity_vs_radius_gen_bkg[ibin][ibin2]->Sumw2();



      Multiplicity_vs_radius_gen_all[ibin][ibin2] = new TH1F((TString) (desc + "_Multiplicity_vs_radius_gen_all"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     Multiplicity_vs_radius_gen_all[ibin][ibin2]->Sumw2();

      Multiplicity_vs_radius_gen_all_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_Multiplicity_vs_radius_gen_all_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     Multiplicity_vs_radius_gen_all_bkg[ibin][ibin2]->Sumw2();




      Multiplicity_vs_radius_ptwt[ibin][ibin2] = new TH1F((TString) (desc + "_Multiplicity_vs_radius_ptwt"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     Multiplicity_vs_radius_ptwt[ibin][ibin2]->Sumw2();

      Multiplicity_vs_radius_bkg_ptwt[ibin][ibin2] = new TH1F((TString) (desc + "_Multiplicity_vs_radius_bkg_ptwt"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     Multiplicity_vs_radius_bkg_ptwt[ibin][ibin2]->Sumw2();


      Multiplicity_vs_radius_gen_ptwt[ibin][ibin2] = new TH1F((TString) (desc + "_Multiplicity_vs_radius_gen_ptwt"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     Multiplicity_vs_radius_gen_ptwt[ibin][ibin2]->Sumw2();

      Multiplicity_vs_radius_gen_bkg_ptwt[ibin][ibin2] = new TH1F((TString) (desc + "_Multiplicity_vs_radius_gen_bkg_ptwt"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     Multiplicity_vs_radius_gen_bkg_ptwt[ibin][ibin2]->Sumw2();


      SumPt_only[ibin][ibin2] = new TH1F((TString) (desc + "_SumPt_only"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 500, 0., 10.);     SumPt_only[ibin][ibin2]->Sumw2();



      JetEnergy_resolution[ibin][ibin2] = new TH1F((TString) (desc + "_JetEnergy_resolution"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -2., 2.);     JetEnergy_resolution[ibin][ibin2]->Sumw2();

      JetEnergy_ratio[ibin][ibin2] = new TH1F((TString) (desc + "_JetEnergy_ratio"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 75, 0., 5.);     JetEnergy_ratio[ibin][ibin2]->Sumw2();

      JetEnergy_ratio_new[ibin][ibin2] = new TH1F((TString) (desc + "_JetEnergy_ratio_new"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 75, 0., 5.);     JetEnergy_ratio_new[ibin][ibin2]->Sumw2();

      JetShapeDiffParticles_1D_EMix[ibin][ibin2] = new TH1F((TString) (desc + "_JetShapeDiffParticles_1D_EMix"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     JetShapeDiffParticles_1D_EMix[ibin][ibin2]->Sumw2();

      emix_jetpt[ibin][ibin2] = new TH1F((TString) (desc + "_emix_jetpt"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, 0., 500.);     emix_jetpt[ibin][ibin2]->Sumw2(); // pt1 bin loop

      emixed_jetpt[ibin][ibin2] = new TH1F((TString) (desc + "_emixed_jetpt"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, 0., 500.);     emixed_jetpt[ibin][ibin2]->Sumw2(); // pt1 bin loop

      emixed_jeteta[ibin][ibin2] = new TH1F((TString) (desc + "_emixed_jeteta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5, 5);     emixed_jeteta[ibin][ibin2]->Sumw2(); // pt1 bin loop


      emixed_jetphi[ibin][ibin2] = new TH1F((TString) (desc + "_emixed_jetphi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 72, -TMath::Pi(), TMath::Pi());     emixed_jetphi[ibin][ibin2]->Sumw2(); // pt1 bin loop


      for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){

	/*hJetTrackSignalBackground[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackSignalBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 200,-5,5,150,-TMath::Pi(),2*TMath::Pi());     hJetTrackSignalBackground[ibin][ibin2][ibin3]->Sumw2();

	hJetTrackBackground[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 200,-5,5,150,-TMath::Pi(),2*TMath::Pi());  hJetTrackBackground[ibin][ibin2][ibin3]->Sumw2();

*/


        hJetTrackSignalBackground[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackSignalBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,75,-TMath::Pi(),2*TMath::Pi());     hJetTrackSignalBackground[ibin][ibin2][ibin3]->Sumw2();

        hJetTrackBackground[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,75,-TMath::Pi(),2*TMath::Pi());  hJetTrackBackground[ibin][ibin2][ibin3]->Sumw2();



	/// without corrections

  hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackSignalBackground_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,75,-TMath::Pi(),2*TMath::Pi());     hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Sumw2();
      
        hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackBackground_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,75,-TMath::Pi(),2*TMath::Pi());  hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3]->Sumw2();


          hJetGenTrackSignalBackground[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetGenTrackSignalBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,75,-TMath::Pi(),2*TMath::Pi());     hJetGenTrackSignalBackground[ibin][ibin2][ibin3]->Sumw2();


            hJetGenTrackBackground[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetGenTrackBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,75,-TMath::Pi(),2*TMath::Pi());     hJetGenTrackBackground[ibin][ibin2][ibin3]->Sumw2();


 hJetGenTrackSignalBackground_pythia[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetGenTrackSignalBackground_pythia"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,75,-TMath::Pi(),2*TMath::Pi());     hJetGenTrackSignalBackground_pythia[ibin][ibin2][ibin3]->Sumw2();


 hJetGenTrackSignalBackground_hydjet[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetGenTrackSignalBackground_hydjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,75,-TMath::Pi(),2*TMath::Pi());     hJetGenTrackSignalBackground_hydjet[ibin][ibin2][ibin3]->Sumw2();


 hJetGenTrackSignalBackground_hydro[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetGenTrackSignalBackground_hydro"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,75,-TMath::Pi(),2*TMath::Pi());     hJetGenTrackSignalBackground_hydro[ibin][ibin2][ibin3]->Sumw2();


 hJetGenTrackSignalBackground_minijets[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetGenTrackSignalBackground_minijets"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,75,-TMath::Pi(),2*TMath::Pi());     hJetGenTrackSignalBackground_minijets[ibin][ibin2][ibin3]->Sumw2();



       hJetGenTrackBackground_pythia[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetGenTrackBackground_pythia"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,75,-TMath::Pi(),2*TMath::Pi());     hJetGenTrackBackground_pythia[ibin][ibin2][ibin3]->Sumw2();


             hJetGenTrackBackground_hydjet[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetGenTrackBackground_hydjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,75,-TMath::Pi(),2*TMath::Pi());     hJetGenTrackBackground_hydjet[ibin][ibin2][ibin3]->Sumw2();



/*

	hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackSignalBackground_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 205,-5,5,150,-TMath::Pi(),2*TMath::Pi());     hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Sumw2();

	hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackBackground_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 205,-5,5,150,-TMath::Pi(),2*TMath::Pi());  hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3]->Sumw2();


	hJetGenTrackSignalBackground[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetGenTrackSignalBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 205,-5,5,150,-TMath::Pi(),2*TMath::Pi());     hJetGenTrackSignalBackground[ibin][ibin2][ibin3]->Sumw2();


	hJetGenTrackBackground[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetGenTrackBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 205,-5,5,150,-TMath::Pi(),2*TMath::Pi());     hJetGenTrackBackground[ibin][ibin2][ibin3]->Sumw2();

*/




  Ngen[ibin][ibin2][ibin3] = new TH1D((TString) (desc + "_Ngen"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100, 0., 100.);     Ngen[ibin][ibin2][ibin3]->Sumw2();

	Ntracks[ibin][ibin2][ibin3] = new TH1D((TString) (desc + "_Ntracks"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 1000, 0., 5000.);     Ntracks[ibin][ibin2][ibin3]->Sumw2();


  Ntracks_awayside[ibin][ibin2][ibin3] = new TH1D((TString) (desc + "_Ntracks_awayside"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 1000, 0., 5000.);     Ntracks_awayside[ibin][ibin2][ibin3]->Sumw2();


	TrkPt[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_TrkPt"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100, 0., 100.);     TrkPt[ibin][ibin2][ibin3]->Sumw2();

	TrkEta[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_TrkEta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 12, -2.4, 2.4 );     TrkEta[ibin][ibin2][ibin3]->Sumw2();

	TrkPhi[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_TrkPhi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 72, -TMath::Pi(), TMath::Pi());     TrkPhi[ibin][ibin2][ibin3]->Sumw2();



TrkPt_weighted[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_TrkPt_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100, 0., 100.);     TrkPt[ibin][ibin2][ibin3]->Sumw2();

  TrkEta_weighted[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_TrkEta_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 12, -2.4, 2.4);     TrkEta[ibin][ibin2][ibin3]->Sumw2();

    TrkPhi_weighted[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_TrkPhi_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 72, -TMath::Pi(), TMath::Pi());     TrkPhi[ibin][ibin2][ibin3]->Sumw2();







       GenTrkPt[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_GenTrkPt"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100, 0., 100.);     GenTrkPt[ibin][ibin2][ibin3]->Sumw2();
     
             GenTrkEta[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_GenTrkEta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 12, -2.4, 2.4);     GenTrkEta[ibin][ibin2][ibin3]->Sumw2();

	             GenTrkPhi[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_GenTrkPhi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 72, -TMath::Pi(), TMath::Pi());     GenTrkPhi[ibin][ibin2][ibin3]->Sumw2();

///// bkg eta, phi for inclusive tracks

    GenTrkEta_bkg[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_GenTrkEta_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 12, -2.4, 2.4);     GenTrkEta_bkg[ibin][ibin2][ibin3]->Sumw2();

                   GenTrkPhi_bkg[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_GenTrkPhi_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 72, -TMath::Pi(), TMath::Pi());     GenTrkPhi_bkg[ibin][ibin2][ibin3]->Sumw2();

 TrkEta_bkg[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_TrkEta_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 12, -2.4, 2.4);     TrkEta_bkg[ibin][ibin2][ibin3]->Sumw2();

   TrkPhi_bkg[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_TrkPhi_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 72, -TMath::Pi(), TMath::Pi());     TrkPhi_bkg[ibin][ibin2][ibin3]->Sumw2();


     
     
      } /// ibin3




    } // pt bin loop
  } // centrality bin loop
} // hist class loop

void hist_class::NormalizeByCounts()
{


  /*	for(int centi = 0; centi < nCBins; centi++) {
	for(int pti = 0; pti < nPtBins; pti++) {

	if( sig_count[centi][pti] < 0.5e-18 ) continue;
	JetShapeDiffParticles_1D[centi][pti]->Scale(1. / sig_count[centi][pti]);
	JetShapeDiffParticlesGen_1D[centi][pti]->Scale(1. / sig_count[centi][pti]);
	Multiplicity_vs_radius[centi][pti]->Scale(1. / sig_count[centi][pti]);
	Multiplicity_vs_radius_gen[centi][pti]->Scale(1. / sig_count[centi][pti]);

	if( bkg_count[centi][pti] < 0.5e-18  ) continue;
	JetShapeDiffParticles_bkg_1D[centi][pti]->Scale(1. / bkg_count[centi][pti]);
	JetShapeDiffParticlesGen_bkg_1D[centi][pti]->Scale(1. / bkg_count[centi][pti]);
	Multiplicity_vs_radius_bkg[centi][pti]->Scale(1. / bkg_count[centi][pti]);
	Multiplicity_vs_radius_gen_bkg[centi][pti]->Scale(1. / bkg_count[centi][pti]);

	}

	for(int pti_sub = 0; pti_sub < nPtBins_sub; pti_sub++) {

	if( sig_count_sub[centi][pti_sub] < 0.5e-18 ) continue;
	JetShapeDiffParticles_subleadingjet_1D[centi][pti_sub]->Scale(1. / sig_count_sub[centi][pti_sub]);

	if( bkg_count_sub[centi][pti_sub] < 0.5e-18  ) continue;
	JetShapeDiffParticles_subleadingjet_bkg_1D[centi][pti_sub]->Scale(1. / bkg_count_sub[centi][pti_sub]);

	}


	}
   */	

}





void hist_class::AddHists(hist_class *more_hists, float wt)
{


  std::cout << "AddHists" << std::endl;
  for(int centi = 0; centi < nCBins; centi++) 
  {
    for(int pti = 0; pti < nPtBins; pti++) 
    {
      std::cout << "centi: " << centi << ", ptio : " << pti << ", wt: " << wt << ", add signal: " << more_hists->sig_count[centi][pti] << ", bkg: " << more_hists->bkg_count[centi][pti] << ", result signal: " << sig_count[centi][pti] << ", result_bkg: " << bkg_count[centi][pti] << "\n";
      sig_count[centi][pti] += more_hists->sig_count[centi][pti]*wt;
      bkg_count[centi][pti] += more_hists->bkg_count[centi][pti]*wt;

    }

    for(int pti_sub = 0; pti_sub < nPtBins_sub; pti_sub++) 
    {
      std::cout << "centi: " << centi << ", ptio_sub : " << pti_sub << ", wt: " << wt << ", add signal: " << more_hists->sig_count_sub[centi][pti_sub] << ", bkg: " << more_hists->bkg_count_sub[centi][pti_sub] << ", result signal: " << sig_count_sub[centi][pti_sub] << ", result_bkg: " << bkg_count[centi][pti_sub] << "\n";

      sig_count_sub[centi][pti_sub] += more_hists->sig_count_sub[centi][pti_sub]*wt;
      bkg_count_sub[centi][pti_sub] += more_hists->bkg_count_sub[centi][pti_sub]*wt;

    }


  }


  //    JetEnergy_gen_vs_rec->Sumw2();   more_hists->JetEnergy_gen_vs_rec->Sumw2();
  //      JetEnergy_gen_vs_rec->Add(more_hists->JetEnergy_gen_vs_rec, wt);

  alljets_rec->Sumw2();   more_hists->alljets_rec->Sumw2();
  alljets_rec->Add(more_hists->alljets_rec, wt);

  NumberOfMatches->Sumw2();   more_hists->NumberOfMatches->Sumw2();
  NumberOfMatches->Add(more_hists->NumberOfMatches, wt);

  NEvents->Sumw2();   more_hists->NEvents->Sumw2();
  NEvents->Add(more_hists->NEvents, wt);

  
  
  
  
  NEvents_test->Add(more_hists->NEvents_test, wt);

  NEvents_after_noise->Add(more_hists->NEvents_after_noise, wt);
  NEvents_after_spike->Add(more_hists->NEvents_after_spike, wt);

  NEvents_after_dphi->Sumw2();   more_hists->NEvents_after_dphi->Sumw2();
  NEvents_after_dphi->Add(more_hists->NEvents_after_dphi, wt);

  NEvents_before_dphi->Sumw2();   more_hists->NEvents_before_dphi->Sumw2();
  NEvents_before_dphi->Add(more_hists->NEvents_before_dphi, wt);

  NEvents_dijets->Sumw2();   more_hists->NEvents_dijets->Sumw2();
  NEvents_dijets->Add(more_hists->NEvents_dijets, wt);

  NEvents_after_trigger->Sumw2();   more_hists->NEvents_after_trigger->Sumw2();
  NEvents_after_trigger->Add(more_hists->NEvents_after_trigger, wt);

  AllJetPt_raw_hist->Sumw2();   more_hists->AllJetPt_raw_hist->Sumw2();
  AllJetPt_raw_hist->Add(more_hists->AllJetPt_raw_hist, wt);
  AllJetPt_hist->Sumw2();   more_hists->AllJetPt_hist->Sumw2();
  AllJetPt_hist->Add(more_hists->AllJetPt_hist, wt);
  AllJetPhi_hist->Sumw2();   more_hists->AllJetPhi_hist->Sumw2();
  AllJetPhi_hist->Add(more_hists->AllJetPhi_hist, wt);
  AllJetEta_hist->Sumw2();   more_hists->AllJetEta_hist->Sumw2();
  AllJetEta_hist->Add(more_hists->AllJetEta_hist, wt);

  First_AllJetPt_hist->Sumw2();   more_hists->First_AllJetPt_hist->Sumw2();
  First_AllJetPt_hist->Add(more_hists->First_AllJetPt_hist, wt);

  First_AllJetPhi_hist->Sumw2();   more_hists->First_AllJetPhi_hist->Sumw2();
  First_AllJetPhi_hist->Add(more_hists->First_AllJetPhi_hist, wt);

  First_AllJetEta_hist->Sumw2();   more_hists->First_AllJetEta_hist->Sumw2();
  First_AllJetEta_hist->Add(more_hists->First_AllJetEta_hist, wt);

  Sub_AllJetPt_hist->Sumw2();   more_hists->Sub_AllJetPt_hist->Sumw2();
  Sub_AllJetPt_hist->Add(more_hists->Sub_AllJetPt_hist, wt);

  Sub_AllJetPhi_hist->Sumw2();   more_hists->Sub_AllJetPhi_hist->Sumw2();
  Sub_AllJetPhi_hist->Add(more_hists->Sub_AllJetPhi_hist, wt);

  Sub_AllJetEta_hist->Sumw2();   more_hists->Sub_AllJetEta_hist->Sumw2();
  Sub_AllJetEta_hist->Add(more_hists->Sub_AllJetEta_hist, wt);

  JetPt_fraction->Sumw2();   more_hists->JetPt_fraction->Sumw2();
  JetPt_fraction->Add(more_hists->JetPt_fraction, wt);


  Centrality->Sumw2();   more_hists->Centrality->Sumw2();
  Centrality->Add(more_hists->Centrality, wt);


  Centrality_weighted->Sumw2();   more_hists->Centrality_weighted->Sumw2();
  Centrality_weighted->Add(more_hists->Centrality_weighted, wt);


  track_vz->Sumw2();   more_hists->track_vz->Sumw2();
  track_vz->Add(more_hists->track_vz, wt);

  track_vz_weighted->Sumw2();   more_hists->track_vz_weighted->Sumw2();
  track_vz_weighted->Add(more_hists->track_vz_weighted, wt);

  for (int ibin2=0;ibin2<nPtBins;ibin2++)
  {
    fullrecjet_aftercuts[ibin2]->Sumw2();   more_hists->fullrecjet_aftercuts[ibin2]->Sumw2();
    fullrecjet_aftercuts[ibin2]->Add(more_hists->fullrecjet_aftercuts[ibin2], wt);

    fullrefjet_aftercuts[ibin2]->Sumw2();   more_hists->fullrefjet_aftercuts[ibin2]->Sumw2();
    fullrefjet_aftercuts[ibin2]->Add(more_hists->fullrefjet_aftercuts[ibin2], wt);


  }

  for (int ibin=0;ibin<nCBins;ibin++){

    jet_pT_hist[ibin]->Sumw2();   more_hists->jet_pT_hist[ibin]->Sumw2();
    jet_pT_hist[ibin]->Add(more_hists->jet_pT_hist[ibin], wt);
    jet_phi_hist[ibin]->Sumw2();   more_hists->jet_phi_hist[ibin]->Sumw2();
    jet_phi_hist[ibin]->Add(more_hists->jet_phi_hist[ibin], wt);
    jet_eta_hist[ibin]->Sumw2();   more_hists->jet_eta_hist[ibin]->Sumw2();
    jet_eta_hist[ibin]->Add(more_hists->jet_eta_hist[ibin], wt);
    jet_corrpT_hist[ibin]->Sumw2();   more_hists->jet_corrpT_hist[ibin]->Sumw2();
    jet_corrpT_hist[ibin]->Add(more_hists->jet_corrpT_hist[ibin], wt);
    LeadingJetPt_hist[ibin]->Sumw2();   more_hists->LeadingJetPt_hist[ibin]->Sumw2();
    LeadingJetPt_hist[ibin]->Add(more_hists->LeadingJetPt_hist[ibin], wt);
    LeadingJetPhi_hist[ibin]->Sumw2();   more_hists->LeadingJetPhi_hist[ibin]->Sumw2();
    LeadingJetPhi_hist[ibin]->Add(more_hists->LeadingJetPhi_hist[ibin], wt);
    LeadingJetEta_hist[ibin]->Sumw2();   more_hists->LeadingJetEta_hist[ibin]->Sumw2();
    LeadingJetEta_hist[ibin]->Add(more_hists->LeadingJetEta_hist[ibin], wt);
    SubJetPt_hist[ibin]->Sumw2();   more_hists->SubJetPt_hist[ibin]->Sumw2();
    SubJetPt_hist[ibin]->Add(more_hists->SubJetPt_hist[ibin], wt);
    SubJetPhi_hist[ibin]->Sumw2();   more_hists->SubJetPhi_hist[ibin]->Sumw2();
    SubJetPhi_hist[ibin]->Add(more_hists->SubJetPhi_hist[ibin], wt);
    SubJetEta_hist[ibin]->Sumw2();   more_hists->SubJetEta_hist[ibin]->Sumw2();
    SubJetEta_hist[ibin]->Add(more_hists->SubJetEta_hist[ibin], wt);

    dPhi_leadingjet_hist[ibin]->Sumw2();   more_hists->dPhi_leadingjet_hist[ibin]->Sumw2();
    dPhi_leadingjet_hist[ibin]->Add(more_hists->dPhi_leadingjet_hist[ibin], wt);
    dPhi_subleadingjet_hist[ibin]->Sumw2();   more_hists->dPhi_subleadingjet_hist[ibin]->Sumw2();
    dPhi_subleadingjet_hist[ibin]->Add(more_hists->dPhi_subleadingjet_hist[ibin], wt);

    JetShapeDiffParticles_bkg_1D_total[ibin]->Sumw2();   more_hists->JetShapeDiffParticles_bkg_1D_total[ibin]->Sumw2();
    JetShapeDiffParticles_bkg_1D_total[ibin]->Add(more_hists->JetShapeDiffParticles_bkg_1D_total[ibin], wt);

    JetShapeIntegratedParticles_bkg_1D_total[ibin]->Sumw2();   more_hists->JetShapeIntegratedParticles_bkg_1D_total[ibin]->Sumw2();
    JetShapeIntegratedParticles_bkg_1D_total[ibin]->Add(more_hists->JetShapeIntegratedParticles_bkg_1D_total[ibin], wt);


    JetEnergy_gen_vs_rec[ibin]->Sumw2();   more_hists->JetEnergy_gen_vs_rec[ibin]->Sumw2();
    JetEnergy_gen_vs_rec[ibin]->Add(more_hists->JetEnergy_gen_vs_rec[ibin], wt);

    for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){

      rr_gen_leadingjet[ibin][ibin3]->Sumw2();   more_hists->rr_gen_leadingjet[ibin][ibin3]->Sumw2();
      rr_gen_leadingjet[ibin][ibin3]->Add(more_hists->rr_gen_leadingjet[ibin][ibin3], wt);


      rr_reco_leadingjet[ibin][ibin3]->Sumw2();   more_hists->rr_reco_leadingjet[ibin][ibin3]->Sumw2();
      rr_reco_leadingjet[ibin][ibin3]->Add(more_hists->rr_reco_leadingjet[ibin][ibin3], wt);

      rr_reco_leadingjet_weighted[ibin][ibin3]->Sumw2();   more_hists->rr_reco_leadingjet_weighted[ibin][ibin3]->Sumw2();
      rr_reco_leadingjet_weighted[ibin][ibin3]->Add(more_hists->rr_reco_leadingjet_weighted[ibin][ibin3], wt);

    }

    dPhi_hist[ibin]->Sumw2();   more_hists->dPhi_hist[ibin]->Sumw2();
    dPhi_hist[ibin]->Add(more_hists->dPhi_hist[ibin], wt);

    dPhi_after_hist[ibin]->Sumw2();   more_hists->dPhi_after_hist[ibin]->Sumw2();
    dPhi_after_hist[ibin]->Add(more_hists->dPhi_after_hist[ibin], wt);


    dPhi_vector[ibin]->Sumw2();   more_hists->dPhi_vector[ibin]->Sumw2();
    dPhi_vector[ibin]->Add(more_hists->dPhi_vector[ibin], wt);

    Aj[ibin]->Sumw2();   more_hists->Aj[ibin]->Sumw2();
    Aj[ibin]->Add(more_hists->Aj[ibin], wt);




    for (int ibin2=0;ibin2<nPtBins;ibin2++)
    { 
      
      
      ThirdJetPt_hist[ibin][ibin2]->Sumw2();   more_hists->ThirdJetPt_hist[ibin][ibin2]->Sumw2();
      ThirdJetPt_hist[ibin][ibin2]->Add(more_hists->ThirdJetPt_hist[ibin][ibin2], wt);
      ThirdJetPhi_hist[ibin][ibin2]->Sumw2();   more_hists->ThirdJetPhi_hist[ibin][ibin2]->Sumw2();
      ThirdJetPhi_hist[ibin][ibin2]->Add(more_hists->ThirdJetPhi_hist[ibin][ibin2], wt);
      ThirdJetEta_hist[ibin][ibin2]->Sumw2();   more_hists->ThirdJetEta_hist[ibin][ibin2]->Sumw2();
      ThirdJetEta_hist[ibin][ibin2]->Add(more_hists->ThirdJetEta_hist[ibin][ibin2], wt);
      //dPhi_hist[ibin][ibin2]->Sumw2();   more_hists->dPhi_hist[ibin][ibin2]->Sumw2();
      //dPhi_hist[ibin][ibin2]->Add(more_hists->dPhi_hist[ibin][ibin2], wt);

      dEta_hist[ibin][ibin2]->Sumw2();   more_hists->dEta_hist[ibin][ibin2]->Sumw2();
      dEta_hist[ibin][ibin2]->Add(more_hists->dEta_hist[ibin][ibin2], wt);


      //dPhi_after_hist[ibin][ibin2]->Sumw2();   more_hists->dPhi_after_hist[ibin][ibin2]->Sumw2();
      //dPhi_after_hist[ibin][ibin2]->Add(more_hists->dPhi_after_hist[ibin][ibin2], wt);



      all_jets_corrpT[ibin][ibin2]->Sumw2();   more_hists->all_jets_corrpT[ibin][ibin2]->Sumw2();
      all_jets_corrpT[ibin][ibin2]->Add(more_hists->all_jets_corrpT[ibin][ibin2], wt);
      all_jets_phi[ibin][ibin2]->Sumw2();   more_hists->all_jets_phi[ibin][ibin2]->Sumw2();
      all_jets_phi[ibin][ibin2]->Add(more_hists->all_jets_phi[ibin][ibin2], wt);
      all_jets_eta[ibin][ibin2]->Sumw2();   more_hists->all_jets_eta[ibin][ibin2]->Sumw2();
      all_jets_eta[ibin][ibin2]->Add(more_hists->all_jets_eta[ibin][ibin2], wt);

      //// leading jet histograms
      only_leadingjets_corrpT[ibin][ibin2]->Sumw2();   more_hists->only_leadingjets_corrpT[ibin][ibin2]->Sumw2();
      only_leadingjets_corrpT[ibin][ibin2]->Add(more_hists->only_leadingjets_corrpT[ibin][ibin2], wt);
      only_leadingjets_phi[ibin][ibin2]->Sumw2();   more_hists->only_leadingjets_phi[ibin][ibin2]->Sumw2();
      only_leadingjets_phi[ibin][ibin2]->Add(more_hists->only_leadingjets_phi[ibin][ibin2], wt);
      only_leadingjets_eta[ibin][ibin2]->Sumw2();   more_hists->only_leadingjets_eta[ibin][ibin2]->Sumw2();
      only_leadingjets_eta[ibin][ibin2]->Add(more_hists->only_leadingjets_eta[ibin][ibin2], wt);
      all_cand_pT_hist[ibin][ibin2]->Sumw2();   more_hists->all_cand_pT_hist[ibin][ibin2]->Sumw2();
      all_cand_pT_hist[ibin][ibin2]->Add(more_hists->all_cand_pT_hist[ibin][ibin2], wt);
      all_cand_phi_hist[ibin][ibin2]->Sumw2();   more_hists->all_cand_phi_hist[ibin][ibin2]->Sumw2();
      all_cand_phi_hist[ibin][ibin2]->Add(more_hists->all_cand_phi_hist[ibin][ibin2], wt);
      all_cand_eta_hist[ibin][ibin2]->Sumw2();   more_hists->all_cand_eta_hist[ibin][ibin2]->Sumw2();
      all_cand_eta_hist[ibin][ibin2]->Add(more_hists->all_cand_eta_hist[ibin][ibin2], wt);

      track_cand_pT_hist[ibin][ibin2]->Sumw2();   more_hists->track_cand_pT_hist[ibin][ibin2]->Sumw2();
      track_cand_pT_hist[ibin][ibin2]->Add(more_hists->track_cand_pT_hist[ibin][ibin2], wt);
      track_cand_phi_hist[ibin][ibin2]->Sumw2();   more_hists->track_cand_phi_hist[ibin][ibin2]->Sumw2();
      track_cand_phi_hist[ibin][ibin2]->Add(more_hists->track_cand_phi_hist[ibin][ibin2], wt);
      track_cand_eta_hist[ibin][ibin2]->Sumw2();   more_hists->track_cand_eta_hist[ibin][ibin2]->Sumw2();
      track_cand_eta_hist[ibin][ibin2]->Add(more_hists->track_cand_eta_hist[ibin][ibin2], wt);


      track_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->track_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();
      track_cand_pT_hist_subleadingjet[ibin][ibin2]->Add(more_hists->track_cand_pT_hist_subleadingjet[ibin][ibin2], wt);

      neutral_cand_pT_hist[ibin][ibin2]->Sumw2();   more_hists->neutral_cand_pT_hist[ibin][ibin2]->Sumw2();
      neutral_cand_pT_hist[ibin][ibin2]->Add(more_hists->neutral_cand_pT_hist[ibin][ibin2], wt);
      neutral_cand_phi_hist[ibin][ibin2]->Sumw2();   more_hists->neutral_cand_phi_hist[ibin][ibin2]->Sumw2();
      neutral_cand_phi_hist[ibin][ibin2]->Add(more_hists->neutral_cand_phi_hist[ibin][ibin2], wt);
      neutral_cand_eta_hist[ibin][ibin2]->Sumw2();   more_hists->neutral_cand_eta_hist[ibin][ibin2]->Sumw2();
      neutral_cand_eta_hist[ibin][ibin2]->Add(more_hists->neutral_cand_eta_hist[ibin][ibin2], wt);
      photons_cand_pT_hist[ibin][ibin2]->Sumw2();   more_hists->photons_cand_pT_hist[ibin][ibin2]->Sumw2();
      photons_cand_pT_hist[ibin][ibin2]->Add(more_hists->photons_cand_pT_hist[ibin][ibin2], wt);
      photons_cand_phi_hist[ibin][ibin2]->Sumw2();   more_hists->photons_cand_phi_hist[ibin][ibin2]->Sumw2();
      photons_cand_phi_hist[ibin][ibin2]->Add(more_hists->photons_cand_phi_hist[ibin][ibin2], wt);
      photons_cand_eta_hist[ibin][ibin2]->Sumw2();   more_hists->photons_cand_eta_hist[ibin][ibin2]->Sumw2();
      photons_cand_eta_hist[ibin][ibin2]->Add(more_hists->photons_cand_eta_hist[ibin][ibin2], wt);


      NumNeutral[ibin][ibin2]->Sumw2();   more_hists->NumNeutral[ibin][ibin2]->Sumw2();
      NumNeutral[ibin][ibin2]->Add(more_hists->NumNeutral[ibin][ibin2], wt);
      NumPhotons[ibin][ibin2]->Sumw2();   more_hists->NumPhotons[ibin][ibin2]->Sumw2();
      NumPhotons[ibin][ibin2]->Add(more_hists->NumPhotons[ibin][ibin2], wt);
      NumAll[ibin][ibin2]->Sumw2();   more_hists->NumAll[ibin][ibin2]->Sumw2();
      NumAll[ibin][ibin2]->Add(more_hists->NumAll[ibin][ibin2], wt);
      NumCharged[ibin][ibin2]->Sumw2();   more_hists->NumCharged[ibin][ibin2]->Sumw2();
      NumCharged[ibin][ibin2]->Add(more_hists->NumCharged[ibin][ibin2], wt);
      NumChargedHadrons[ibin][ibin2]->Sumw2();   more_hists->NumChargedHadrons[ibin][ibin2]->Sumw2();
      NumChargedHadrons[ibin][ibin2]->Add(more_hists->NumChargedHadrons[ibin][ibin2], wt);
      NumChargedParticles[ibin][ibin2]->Sumw2();   more_hists->NumChargedParticles[ibin][ibin2]->Sumw2();
      NumChargedParticles[ibin][ibin2]->Add(more_hists->NumChargedParticles[ibin][ibin2], wt);
      NumElectrons[ibin][ibin2]->Sumw2();   more_hists->NumElectrons[ibin][ibin2]->Sumw2();
      NumElectrons[ibin][ibin2]->Add(more_hists->NumElectrons[ibin][ibin2], wt);
      NumMuons[ibin][ibin2]->Sumw2();   more_hists->NumMuons[ibin][ibin2]->Sumw2();
      NumMuons[ibin][ibin2]->Add(more_hists->NumMuons[ibin][ibin2], wt);


      NumNeutral_bkg[ibin][ibin2]->Sumw2();   more_hists->NumNeutral_bkg[ibin][ibin2]->Sumw2();
      NumNeutral_bkg[ibin][ibin2]->Add(more_hists->NumNeutral_bkg[ibin][ibin2], wt);
      NumPhotons_bkg[ibin][ibin2]->Sumw2();   more_hists->NumPhotons_bkg[ibin][ibin2]->Sumw2();
      NumPhotons_bkg[ibin][ibin2]->Add(more_hists->NumPhotons_bkg[ibin][ibin2], wt);
      NumAll_bkg[ibin][ibin2]->Sumw2();   more_hists->NumAll_bkg[ibin][ibin2]->Sumw2();
      NumAll_bkg[ibin][ibin2]->Add(more_hists->NumAll_bkg[ibin][ibin2], wt);
      NumChargedHadrons_bkg[ibin][ibin2]->Sumw2();   more_hists->NumChargedHadrons_bkg[ibin][ibin2]->Sumw2();
      NumChargedHadrons_bkg[ibin][ibin2]->Add(more_hists->NumChargedHadrons_bkg[ibin][ibin2], wt);
      NumChargedParticles_bkg[ibin][ibin2]->Sumw2();   more_hists->NumChargedParticles_bkg[ibin][ibin2]->Sumw2();
      NumChargedParticles_bkg[ibin][ibin2]->Add(more_hists->NumChargedParticles_bkg[ibin][ibin2], wt);
      NumElectrons_bkg[ibin][ibin2]->Sumw2();   more_hists->NumElectrons_bkg[ibin][ibin2]->Sumw2();
      NumElectrons_bkg[ibin][ibin2]->Add(more_hists->NumElectrons_bkg[ibin][ibin2], wt);
      NumMuons_bkg[ibin][ibin2]->Sumw2();   more_hists->NumMuons_bkg[ibin][ibin2]->Sumw2();
      NumMuons_bkg[ibin][ibin2]->Add(more_hists->NumMuons_bkg[ibin][ibin2], wt);


      chargedhadrons_cand_pT_hist[ibin][ibin2]->Sumw2();   more_hists->chargedhadrons_cand_pT_hist[ibin][ibin2]->Sumw2();
      chargedhadrons_cand_pT_hist[ibin][ibin2]->Add(more_hists->chargedhadrons_cand_pT_hist[ibin][ibin2], wt);
      chargedhadrons_cand_phi_hist[ibin][ibin2]->Sumw2();   more_hists->chargedhadrons_cand_phi_hist[ibin][ibin2]->Sumw2();
      chargedhadrons_cand_phi_hist[ibin][ibin2]->Add(more_hists->chargedhadrons_cand_phi_hist[ibin][ibin2], wt);
      chargedhadrons_cand_eta_hist[ibin][ibin2]->Sumw2();   more_hists->chargedhadrons_cand_eta_hist[ibin][ibin2]->Sumw2();
      chargedhadrons_cand_eta_hist[ibin][ibin2]->Add(more_hists->chargedhadrons_cand_eta_hist[ibin][ibin2], wt);

      chargedparticles_cand_pT_hist[ibin][ibin2]->Sumw2();   more_hists->chargedparticles_cand_pT_hist[ibin][ibin2]->Sumw2();
      chargedparticles_cand_pT_hist[ibin][ibin2]->Add(more_hists->chargedparticles_cand_pT_hist[ibin][ibin2], wt);
      chargedparticles_cand_phi_hist[ibin][ibin2]->Sumw2();   more_hists->chargedparticles_cand_phi_hist[ibin][ibin2]->Sumw2();
      chargedparticles_cand_phi_hist[ibin][ibin2]->Add(more_hists->chargedparticles_cand_phi_hist[ibin][ibin2], wt);
      chargedparticles_cand_eta_hist[ibin][ibin2]->Sumw2();   more_hists->chargedparticles_cand_eta_hist[ibin][ibin2]->Sumw2();
      chargedparticles_cand_eta_hist[ibin][ibin2]->Add(more_hists->chargedparticles_cand_eta_hist[ibin][ibin2], wt);

      electrons_cand_pT_hist[ibin][ibin2]->Sumw2();   more_hists->electrons_cand_pT_hist[ibin][ibin2]->Sumw2();
      electrons_cand_pT_hist[ibin][ibin2]->Add(more_hists->electrons_cand_pT_hist[ibin][ibin2], wt);
      electrons_cand_phi_hist[ibin][ibin2]->Sumw2();   more_hists->electrons_cand_phi_hist[ibin][ibin2]->Sumw2();
      electrons_cand_phi_hist[ibin][ibin2]->Add(more_hists->electrons_cand_phi_hist[ibin][ibin2], wt);
      electrons_cand_eta_hist[ibin][ibin2]->Sumw2();   more_hists->electrons_cand_eta_hist[ibin][ibin2]->Sumw2();
      electrons_cand_eta_hist[ibin][ibin2]->Add(more_hists->electrons_cand_eta_hist[ibin][ibin2], wt);

      muons_cand_pT_hist[ibin][ibin2]->Sumw2();   more_hists->muons_cand_pT_hist[ibin][ibin2]->Sumw2();
      muons_cand_pT_hist[ibin][ibin2]->Add(more_hists->muons_cand_pT_hist[ibin][ibin2], wt);
      muons_cand_phi_hist[ibin][ibin2]->Sumw2();   more_hists->muons_cand_phi_hist[ibin][ibin2]->Sumw2();
      muons_cand_phi_hist[ibin][ibin2]->Add(more_hists->muons_cand_phi_hist[ibin][ibin2], wt);
      muons_cand_eta_hist[ibin][ibin2]->Sumw2();   more_hists->muons_cand_eta_hist[ibin][ibin2]->Sumw2();
      muons_cand_eta_hist[ibin][ibin2]->Add(more_hists->muons_cand_eta_hist[ibin][ibin2], wt);
      /////


      muons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();   more_hists->muons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();
      muons_cand_pT_hist_bkg[ibin][ibin2]->Add(more_hists->muons_cand_pT_hist_bkg[ibin][ibin2], wt);
      electrons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();   more_hists->electrons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();
      electrons_cand_pT_hist_bkg[ibin][ibin2]->Add(more_hists->electrons_cand_pT_hist_bkg[ibin][ibin2], wt);
      chargedparticles_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();   more_hists->chargedparticles_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();
      chargedparticles_cand_pT_hist_bkg[ibin][ibin2]->Add(more_hists->chargedparticles_cand_pT_hist_bkg[ibin][ibin2], wt);
      chargedhadrons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();   more_hists->chargedhadrons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();
      chargedhadrons_cand_pT_hist_bkg[ibin][ibin2]->Add(more_hists->chargedhadrons_cand_pT_hist_bkg[ibin][ibin2], wt);
      neutral_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();   more_hists->neutral_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();
      neutral_cand_pT_hist_bkg[ibin][ibin2]->Add(more_hists->neutral_cand_pT_hist_bkg[ibin][ibin2], wt);
      photons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();   more_hists->photons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();
      photons_cand_pT_hist_bkg[ibin][ibin2]->Add(more_hists->photons_cand_pT_hist_bkg[ibin][ibin2], wt);


      SumJetPtFraction_hist[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist[ibin][ibin2], wt);

      JetShapeIntegratedParticles[ibin][ibin2]->Sumw2();   more_hists->JetShapeIntegratedParticles[ibin][ibin2]->Sumw2();
      JetShapeIntegratedParticles[ibin][ibin2]->Add(more_hists->JetShapeIntegratedParticles[ibin][ibin2], wt);
      JetShapeDiffParticles[ibin][ibin2]->Sumw2();   more_hists->JetShapeDiffParticles[ibin][ibin2]->Sumw2();
      JetShapeDiffParticles[ibin][ibin2]->Add(more_hists->JetShapeDiffParticles[ibin][ibin2], wt);

      JetShapeIntegratedParticles_bkgsub[ibin][ibin2]->Sumw2();   more_hists->JetShapeIntegratedParticles_bkgsub[ibin][ibin2]->Sumw2();
      JetShapeIntegratedParticles_bkgsub[ibin][ibin2]->Add(more_hists->JetShapeIntegratedParticles_bkgsub[ibin][ibin2], wt);

      JetShapeDiffParticles_bkgsub[ibin][ibin2]->Sumw2();   more_hists->JetShapeDiffParticles_bkgsub[ibin][ibin2]->Sumw2();
      JetShapeDiffParticles_bkgsub[ibin][ibin2]->Add(more_hists->JetShapeDiffParticles_bkgsub[ibin][ibin2], wt);


      JetShapeIntegratedParticles_bkg[ibin][ibin2]->Sumw2();   more_hists->JetShapeIntegratedParticles_bkg[ibin][ibin2]->Sumw2();
      JetShapeIntegratedParticles_bkg[ibin][ibin2]->Add(more_hists->JetShapeIntegratedParticles_bkg[ibin][ibin2], wt);

      JetShapeDiffParticles_bkg[ibin][ibin2]->Sumw2();   more_hists->JetShapeDiffParticles_bkg[ibin][ibin2]->Sumw2();
      JetShapeDiffParticles_bkg[ibin][ibin2]->Add(more_hists->JetShapeDiffParticles_bkg[ibin][ibin2], wt);


      JetShapeIntegratedParticles_bkg_1D[ibin][ibin2]->Sumw2();   more_hists->JetShapeIntegratedParticles_bkg_1D[ibin][ibin2]->Sumw2();
      JetShapeIntegratedParticles_bkg_1D[ibin][ibin2]->Add(more_hists->JetShapeIntegratedParticles_bkg_1D[ibin][ibin2], wt);
      JetShapeDiffParticles_bkg_1D[ibin][ibin2]->Sumw2();   more_hists->JetShapeDiffParticles_bkg_1D[ibin][ibin2]->Sumw2();
      JetShapeDiffParticles_bkg_1D[ibin][ibin2]->Add(more_hists->JetShapeDiffParticles_bkg_1D[ibin][ibin2], wt);

      JetShapeIntegratedParticles_1D[ibin][ibin2]->Sumw2();   more_hists->JetShapeIntegratedParticles_1D[ibin][ibin2]->Sumw2();
      JetShapeIntegratedParticles_1D[ibin][ibin2]->Add(more_hists->JetShapeIntegratedParticles_1D[ibin][ibin2], wt);
      JetShapeDiffParticles_1D[ibin][ibin2]->Sumw2();   more_hists->JetShapeDiffParticles_1D[ibin][ibin2]->Sumw2();
      JetShapeDiffParticles_1D[ibin][ibin2]->Add(more_hists->JetShapeDiffParticles_1D[ibin][ibin2], wt);


      JetShapeDiffParticlesGen_1D[ibin][ibin2]->Sumw2();   more_hists->JetShapeDiffParticlesGen_1D[ibin][ibin2]->Sumw2();
      JetShapeDiffParticlesGen_1D[ibin][ibin2]->Add(more_hists->JetShapeDiffParticlesGen_1D[ibin][ibin2], wt);

      JetShapeDiffParticlesGen_bkg_1D[ibin][ibin2]->Sumw2();   more_hists->JetShapeDiffParticlesGen_bkg_1D[ibin][ibin2]->Sumw2();
      JetShapeDiffParticlesGen_bkg_1D[ibin][ibin2]->Add(more_hists->JetShapeDiffParticlesGen_bkg_1D[ibin][ibin2], wt);



      radius_hist[ibin][ibin2]->Sumw2();   more_hists->radius_hist[ibin][ibin2]->Sumw2();
      radius_hist[ibin][ibin2]->Add(more_hists->radius_hist[ibin][ibin2], wt);

      radius_hist_mine[ibin][ibin2]->Sumw2();   more_hists->radius_hist_mine[ibin][ibin2]->Sumw2();
      radius_hist_mine[ibin][ibin2]->Add(more_hists->radius_hist_mine[ibin][ibin2], wt);

      radius_hist_mine2[ibin][ibin2]->Sumw2();   more_hists->radius_hist_mine2[ibin][ibin2]->Sumw2();
      radius_hist_mine2[ibin][ibin2]->Add(more_hists->radius_hist_mine2[ibin][ibin2], wt);


      radius_hist_mine3[ibin][ibin2]->Sumw2();   more_hists->radius_hist_mine3[ibin][ibin2]->Sumw2();
      radius_hist_mine3[ibin][ibin2]->Add(more_hists->radius_hist_mine3[ibin][ibin2], wt);






      //// subleading jet histograms
      only_subleadingjets_corrpT[ibin][ibin2]->Sumw2();   more_hists->only_subleadingjets_corrpT[ibin][ibin2]->Sumw2();
      only_subleadingjets_corrpT[ibin][ibin2]->Add(more_hists->only_subleadingjets_corrpT[ibin][ibin2], wt);
      only_subleadingjets_phi[ibin][ibin2]->Sumw2();   more_hists->only_subleadingjets_phi[ibin][ibin2]->Sumw2();
      only_subleadingjets_phi[ibin][ibin2]->Add(more_hists->only_subleadingjets_phi[ibin][ibin2], wt);
      only_subleadingjets_eta[ibin][ibin2]->Sumw2();   more_hists->only_subleadingjets_eta[ibin][ibin2]->Sumw2();
      only_subleadingjets_eta[ibin][ibin2]->Add(more_hists->only_subleadingjets_eta[ibin][ibin2], wt);
      all_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->all_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();
      all_cand_pT_hist_subleadingjet[ibin][ibin2]->Add(more_hists->all_cand_pT_hist_subleadingjet[ibin][ibin2], wt);
      all_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->all_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();
      all_cand_phi_hist_subleadingjet[ibin][ibin2]->Add(more_hists->all_cand_phi_hist_subleadingjet[ibin][ibin2], wt);
      all_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->all_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();
      all_cand_eta_hist_subleadingjet[ibin][ibin2]->Add(more_hists->all_cand_eta_hist_subleadingjet[ibin][ibin2], wt);
      neutral_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->neutral_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();
      neutral_cand_pT_hist_subleadingjet[ibin][ibin2]->Add(more_hists->neutral_cand_pT_hist_subleadingjet[ibin][ibin2], wt);
      neutral_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->neutral_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();
      neutral_cand_phi_hist_subleadingjet[ibin][ibin2]->Add(more_hists->neutral_cand_phi_hist_subleadingjet[ibin][ibin2], wt);
      neutral_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->neutral_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();
      neutral_cand_eta_hist_subleadingjet[ibin][ibin2]->Add(more_hists->neutral_cand_eta_hist_subleadingjet[ibin][ibin2], wt);


      photons_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->photons_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();
      photons_cand_pT_hist_subleadingjet[ibin][ibin2]->Add(more_hists->photons_cand_pT_hist_subleadingjet[ibin][ibin2], wt);
      photons_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->photons_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();
      photons_cand_phi_hist_subleadingjet[ibin][ibin2]->Add(more_hists->photons_cand_phi_hist_subleadingjet[ibin][ibin2], wt);
      photons_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->photons_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();
      photons_cand_eta_hist_subleadingjet[ibin][ibin2]->Add(more_hists->photons_cand_eta_hist_subleadingjet[ibin][ibin2], wt);

      chargedhadrons_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->chargedhadrons_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();
      chargedhadrons_cand_pT_hist_subleadingjet[ibin][ibin2]->Add(more_hists->chargedhadrons_cand_pT_hist_subleadingjet[ibin][ibin2], wt);
      chargedhadrons_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->chargedhadrons_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();
      chargedhadrons_cand_phi_hist_subleadingjet[ibin][ibin2]->Add(more_hists->chargedhadrons_cand_phi_hist_subleadingjet[ibin][ibin2], wt);
      chargedhadrons_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->chargedhadrons_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();
      chargedhadrons_cand_eta_hist_subleadingjet[ibin][ibin2]->Add(more_hists->chargedhadrons_cand_eta_hist_subleadingjet[ibin][ibin2], wt);


      chargedhadrons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Sumw2();   more_hists->chargedhadrons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Sumw2();
      chargedhadrons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Add(more_hists->chargedhadrons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2], wt);

      neutral_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Sumw2();   more_hists->neutral_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Sumw2();
      neutral_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Add(more_hists->neutral_cand_pT_hist_subleadingjet_bkg[ibin][ibin2], wt);

      photons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Sumw2();   more_hists->photons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Sumw2();
      photons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Add(more_hists->photons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2], wt);

      NumNeutral_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->NumNeutral_subleadingjet[ibin][ibin2]->Sumw2();
      NumNeutral_subleadingjet[ibin][ibin2]->Add(more_hists->NumNeutral_subleadingjet[ibin][ibin2], wt);
      NumPhotons_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->NumPhotons_subleadingjet[ibin][ibin2]->Sumw2();
      NumPhotons_subleadingjet[ibin][ibin2]->Add(more_hists->NumPhotons_subleadingjet[ibin][ibin2], wt);
      NumAll_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->NumAll_subleadingjet[ibin][ibin2]->Sumw2();
      NumAll_subleadingjet[ibin][ibin2]->Add(more_hists->NumAll_subleadingjet[ibin][ibin2], wt);
      NumChargedHadrons_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->NumChargedHadrons_subleadingjet[ibin][ibin2]->Sumw2();
      NumChargedHadrons_subleadingjet[ibin][ibin2]->Add(more_hists->NumChargedHadrons_subleadingjet[ibin][ibin2], wt);

      NumNeutral_subleadingjet_bkg[ibin][ibin2]->Sumw2();   more_hists->NumNeutral_subleadingjet_bkg[ibin][ibin2]->Sumw2();
      NumNeutral_subleadingjet_bkg[ibin][ibin2]->Add(more_hists->NumNeutral_subleadingjet_bkg[ibin][ibin2], wt);
      NumPhotons_subleadingjet_bkg[ibin][ibin2]->Sumw2();   more_hists->NumPhotons_subleadingjet_bkg[ibin][ibin2]->Sumw2();
      NumPhotons_subleadingjet_bkg[ibin][ibin2]->Add(more_hists->NumPhotons_subleadingjet_bkg[ibin][ibin2], wt);
      NumAll_subleadingjet_bkg[ibin][ibin2]->Sumw2();   more_hists->NumAll_subleadingjet_bkg[ibin][ibin2]->Sumw2();
      NumAll_subleadingjet_bkg[ibin][ibin2]->Add(more_hists->NumAll_subleadingjet_bkg[ibin][ibin2], wt);
      NumChargedHadrons_subleadingjet_bkg[ibin][ibin2]->Sumw2();   more_hists->NumChargedHadrons_subleadingjet_bkg[ibin][ibin2]->Sumw2();
      NumChargedHadrons_subleadingjet_bkg[ibin][ibin2]->Add(more_hists->NumChargedHadrons_subleadingjet_bkg[ibin][ibin2], wt);




      SumJetPtFraction_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_subleadingjet[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_subleadingjet[ibin][ibin2], wt);

      SumJetPtFraction_hist_subleadingjet_photons[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_subleadingjet_photons[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_photons[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_subleadingjet_photons[ibin][ibin2], wt);

      SumJetPtFraction_hist_subleadingjet_chargedhadrons[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_subleadingjet_chargedhadrons[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_chargedhadrons[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_subleadingjet_chargedhadrons[ibin][ibin2], wt);

      SumJetPtFraction_hist_subleadingjet_neutrals[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_subleadingjet_neutrals[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_neutrals[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_subleadingjet_neutrals[ibin][ibin2], wt);



      /// bkg..
      SumJetPtFraction_hist_subleadingjet_photons_bkg[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_subleadingjet_photons_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_photons_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_subleadingjet_photons_bkg[ibin][ibin2], wt);

      SumJetPtFraction_hist_subleadingjet_chargedhadrons_bkg[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_subleadingjet_chargedhadrons_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_chargedhadrons_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_subleadingjet_chargedhadrons_bkg[ibin][ibin2], wt);

      SumJetPtFraction_hist_subleadingjet_neutrals_bkg[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_subleadingjet_neutrals_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_neutrals_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_subleadingjet_neutrals_bkg[ibin][ibin2], wt);

      SumJetPtFraction_hist_subleadingjet_chargedparticles_bkg[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_subleadingjet_chargedparticles_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_chargedparticles_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_subleadingjet_chargedparticles_bkg[ibin][ibin2], wt);


      SumJetPtFraction_hist_subleadingjet_electrons_bkg[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_subleadingjet_electrons_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_electrons_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_subleadingjet_electrons_bkg[ibin][ibin2], wt);

      SumJetPtFraction_hist_subleadingjet_muons_bkg[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_subleadingjet_muons_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_muons_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_subleadingjet_muons_bkg[ibin][ibin2], wt);


      SumJetPtFraction_hist_leadingjet_photons[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_photons[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_photons[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_photons[ibin][ibin2], wt);
      SumJetPtFraction_hist_leadingjet_neutrals[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_neutrals[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_neutrals[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_neutrals[ibin][ibin2], wt);

      SumJetPtFraction_hist_leadingjet_charged[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_charged[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_charged[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_charged[ibin][ibin2], wt);

      ///

      SumJetPtFraction_hist_leadingjet_chargedhadrons[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_chargedhadrons[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_chargedhadrons[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_chargedhadrons[ibin][ibin2], wt);

      SumJetPtFraction_hist_leadingjet_chargedparticles[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_chargedparticles[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_chargedparticles[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_chargedparticles[ibin][ibin2], wt);

      SumJetPtFraction_hist_leadingjet_electrons[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_electrons[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_electrons[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_electrons[ibin][ibin2], wt);

      SumJetPtFraction_hist_leadingjet_muons[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_muons[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_muons[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_muons[ibin][ibin2], wt);
      ////

      /// starts -- bkg sumpt fraction plots
      SumJetPtFraction_hist_leadingjet_photons_bkg[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_photons_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_photons_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_photons_bkg[ibin][ibin2], wt);

      SumJetPtFraction_hist_leadingjet_neutrals_bkg[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_neutrals_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_neutrals_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_neutrals_bkg[ibin][ibin2], wt);

      SumJetPtFraction_hist_leadingjet_chargedhadrons_bkg[ibin][ibin2]->Sumw2();more_hists->SumJetPtFraction_hist_leadingjet_chargedhadrons_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_chargedhadrons_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_chargedhadrons_bkg[ibin][ibin2], wt);

      SumJetPtFraction_hist_leadingjet_chargedparticles_bkg[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_chargedparticles_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_chargedparticles_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_chargedparticles_bkg[ibin][ibin2], wt);

      SumJetPtFraction_hist_leadingjet_electrons_bkg[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_electrons_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_electrons_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_electrons_bkg[ibin][ibin2], wt);

      SumJetPtFraction_hist_leadingjet_muons_bkg[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_muons_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_muons_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_muons_bkg[ibin][ibin2], wt);
      /// end -- bkg sumpt fraction plots




      radius_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->radius_hist_subleadingjet[ibin][ibin2]->Sumw2();
      radius_hist_subleadingjet[ibin][ibin2]->Add(more_hists->radius_hist_subleadingjet[ibin][ibin2], wt);
      Centrality_hist[ibin][ibin2]->Sumw2();   more_hists->Centrality_hist[ibin][ibin2]->Sumw2();
      Centrality_hist[ibin][ibin2]->Add(more_hists->Centrality_hist[ibin][ibin2], wt);


      //Aj[ibin][ibin2]->Sumw2();   more_hists->Aj[ibin][ibin2]->Sumw2();
     /// Aj[ibin][ibin2]->Add(more_hists->Aj[ibin][ibin2], wt);


      dN_tracks[ibin][ibin2]->Sumw2();   more_hists->dN_tracks[ibin][ibin2]->Sumw2();
      dN_tracks[ibin][ibin2]->Add(more_hists->dN_tracks[ibin][ibin2], wt);
      dN_chargedhadrons[ibin][ibin2]->Sumw2();   more_hists->dN_chargedhadrons[ibin][ibin2]->Sumw2();
      dN_chargedhadrons[ibin][ibin2]->Add(more_hists->dN_chargedhadrons[ibin][ibin2], wt);
      dN_chargedparticles[ibin][ibin2]->Sumw2();   more_hists->dN_chargedparticles[ibin][ibin2]->Sumw2();
      dN_chargedparticles[ibin][ibin2]->Add(more_hists->dN_chargedparticles[ibin][ibin2], wt);
      dN_electrons[ibin][ibin2]->Sumw2();   more_hists->dN_electrons[ibin][ibin2]->Sumw2();
      dN_electrons[ibin][ibin2]->Add(more_hists->dN_electrons[ibin][ibin2], wt);
      dN_muons[ibin][ibin2]->Sumw2();   more_hists->dN_muons[ibin][ibin2]->Sumw2();
      dN_muons[ibin][ibin2]->Add(more_hists->dN_muons[ibin][ibin2], wt);
      dN_neutrals[ibin][ibin2]->Sumw2();   more_hists->dN_neutrals[ibin][ibin2]->Sumw2();
      dN_neutrals[ibin][ibin2]->Add(more_hists->dN_neutrals[ibin][ibin2], wt);
      dN_photons[ibin][ibin2]->Sumw2();   more_hists->dN_photons[ibin][ibin2]->Sumw2();
      dN_photons[ibin][ibin2]->Add(more_hists->dN_photons[ibin][ibin2], wt);

      JetPt_fraction_hist[ibin][ibin2]->Sumw2();   more_hists->JetPt_fraction_hist[ibin][ibin2]->Sumw2();
      JetPt_fraction_hist[ibin][ibin2]->Add(more_hists->JetPt_fraction_hist[ibin][ibin2], wt);

      dPhi_jet_track[ibin][ibin2]->Sumw2();   more_hists->dPhi_jet_track[ibin][ibin2]->Sumw2();
      dPhi_jet_track[ibin][ibin2]->Add(more_hists->dPhi_jet_track[ibin][ibin2], wt);

      dPhi_jet_track_ptweight[ibin][ibin2]->Sumw2();   more_hists->dPhi_jet_track_ptweight[ibin][ibin2]->Sumw2();
      dPhi_jet_track_ptweight[ibin][ibin2]->Add(more_hists->dPhi_jet_track_ptweight[ibin][ibin2], wt);


      track_cand_pT_hist_weighted[ibin][ibin2]->Sumw2();   more_hists->track_cand_pT_hist_weighted[ibin][ibin2]->Sumw2();
      track_cand_pT_hist_weighted[ibin][ibin2]->Add(more_hists->track_cand_pT_hist_weighted[ibin][ibin2], wt);

      track_cand_eta_hist_weighted[ibin][ibin2]->Sumw2();   more_hists->track_cand_eta_hist_weighted[ibin][ibin2]->Sumw2();
      track_cand_eta_hist_weighted[ibin][ibin2]->Add(more_hists->track_cand_eta_hist_weighted[ibin][ibin2], wt);

      track_cand_phi_hist_weighted[ibin][ibin2]->Sumw2();   more_hists->track_cand_phi_hist_weighted[ibin][ibin2]->Sumw2();
      track_cand_phi_hist_weighted[ibin][ibin2]->Add(more_hists->track_cand_phi_hist_weighted[ibin][ibin2], wt);


      track_bkg_pT_hist_weighted[ibin][ibin2]->Sumw2();   more_hists->track_bkg_pT_hist_weighted[ibin][ibin2]->Sumw2();
      track_bkg_pT_hist_weighted[ibin][ibin2]->Add(more_hists->track_bkg_pT_hist_weighted[ibin][ibin2], wt);

      track_bkg_pT_hist[ibin][ibin2]->Sumw2();   more_hists->track_bkg_pT_hist[ibin][ibin2]->Sumw2();
      track_bkg_pT_hist[ibin][ibin2]->Add(more_hists->track_bkg_pT_hist[ibin][ibin2], wt);
      
      track_gen_bkg_pT_hist[ibin][ibin2]->Sumw2();   more_hists->track_gen_bkg_pT_hist[ibin][ibin2]->Sumw2();
      track_gen_bkg_pT_hist[ibin][ibin2]->Add(more_hists->track_gen_bkg_pT_hist[ibin][ibin2], wt);


      track_bkg_pT_hist_weighted_oocone[ibin][ibin2]->Sumw2();   more_hists->track_bkg_pT_hist_weighted_oocone[ibin][ibin2]->Sumw2();
      track_bkg_pT_hist_weighted_oocone[ibin][ibin2]->Add(more_hists->track_bkg_pT_hist_weighted_oocone[ibin][ibin2], wt);

      track_bkg_pT_hist_oocone[ibin][ibin2]->Sumw2();   more_hists->track_bkg_pT_hist_oocone[ibin][ibin2]->Sumw2();
      track_bkg_pT_hist_oocone[ibin][ibin2]->Add(more_hists->track_bkg_pT_hist_oocone[ibin][ibin2], wt);

      track_gen_bkg_pT_hist_oocone[ibin][ibin2]->Sumw2();   more_hists->track_gen_bkg_pT_hist_oocone[ibin][ibin2]->Sumw2();
      track_gen_bkg_pT_hist_oocone[ibin][ibin2]->Add(more_hists->track_gen_bkg_pT_hist_oocone[ibin][ibin2], wt);

      EtaRef_bkg_pt[ibin][ibin2]->Sumw2(); more_hists->EtaRef_bkg_pt[ibin][ibin2]->Sumw2();
      EtaRef_bkg_pt[ibin][ibin2]->Add(more_hists->EtaRef_bkg_pt[ibin][ibin2], wt);

      EtaRef_bkg_pt_weighted[ibin][ibin2]->Sumw2(); more_hists->EtaRef_bkg_pt_weighted[ibin][ibin2]->Sumw2();
      EtaRef_bkg_pt_weighted[ibin][ibin2]->Add(more_hists->EtaRef_bkg_pt_weighted[ibin][ibin2], wt);


      gensube_hist[ibin][ibin2]->Sumw2();   more_hists->gensube_hist[ibin][ibin2]->Sumw2();
      gensube_hist[ibin][ibin2]->Add(more_hists->gensube_hist[ibin][ibin2], wt);

      genparticles_pT_hist[ibin][ibin2]->Sumw2();   more_hists->genparticles_pT_hist[ibin][ibin2]->Sumw2();
      genparticles_pT_hist[ibin][ibin2]->Add(more_hists->genparticles_pT_hist[ibin][ibin2], wt);

      genparticles_bkg_pT_hist[ibin][ibin2]->Sumw2();   more_hists->genparticles_bkg_pT_hist[ibin][ibin2]->Sumw2();
      genparticles_bkg_pT_hist[ibin][ibin2]->Add(more_hists->genparticles_bkg_pT_hist[ibin][ibin2], wt);


      track_bkg_eta_hist_weighted[ibin][ibin2]->Sumw2();   more_hists->track_bkg_eta_hist_weighted[ibin][ibin2]->Sumw2();
      track_bkg_eta_hist_weighted[ibin][ibin2]->Add(more_hists->track_bkg_eta_hist_weighted[ibin][ibin2], wt);

      track_bkg_eta_hist[ibin][ibin2]->Sumw2();   more_hists->track_bkg_eta_hist[ibin][ibin2]->Sumw2();
      track_bkg_eta_hist[ibin][ibin2]->Add(more_hists->track_bkg_eta_hist[ibin][ibin2], wt);

      track_gen_bkg_eta_hist[ibin][ibin2]->Sumw2();   more_hists->track_gen_bkg_eta_hist[ibin][ibin2]->Sumw2();
      track_gen_bkg_eta_hist[ibin][ibin2]->Add(more_hists->track_gen_bkg_eta_hist[ibin][ibin2], wt);





      Multiplicity_vs_radius[ibin][ibin2]->Sumw2();   more_hists->Multiplicity_vs_radius[ibin][ibin2]->Sumw2();
      Multiplicity_vs_radius[ibin][ibin2]->Add(more_hists->Multiplicity_vs_radius[ibin][ibin2], wt);


      Multiplicity_vs_radius_bkg[ibin][ibin2]->Sumw2();   more_hists->Multiplicity_vs_radius_bkg[ibin][ibin2]->Sumw2();
      Multiplicity_vs_radius_bkg[ibin][ibin2]->Add(more_hists->Multiplicity_vs_radius_bkg[ibin][ibin2], wt);

      Multiplicity_vs_radius_gen[ibin][ibin2]->Sumw2();   more_hists->Multiplicity_vs_radius_gen[ibin][ibin2]->Sumw2();
      Multiplicity_vs_radius_gen[ibin][ibin2]->Add(more_hists->Multiplicity_vs_radius_gen[ibin][ibin2], wt);


      Multiplicity_vs_radius_gen_bkg[ibin][ibin2]->Sumw2();   more_hists->Multiplicity_vs_radius_gen_bkg[ibin][ibin2]->Sumw2();
      Multiplicity_vs_radius_gen_bkg[ibin][ibin2]->Add(more_hists->Multiplicity_vs_radius_gen_bkg[ibin][ibin2], wt);


      Multiplicity_vs_radius_gen_all[ibin][ibin2]->Sumw2();   more_hists->Multiplicity_vs_radius_gen_all[ibin][ibin2]->Sumw2();
      Multiplicity_vs_radius_gen_all[ibin][ibin2]->Add(more_hists->Multiplicity_vs_radius_gen_all[ibin][ibin2], wt);


      Multiplicity_vs_radius_gen_all_bkg[ibin][ibin2]->Sumw2();   more_hists->Multiplicity_vs_radius_gen_all_bkg[ibin][ibin2]->Sumw2();
      Multiplicity_vs_radius_gen_all_bkg[ibin][ibin2]->Add(more_hists->Multiplicity_vs_radius_gen_all_bkg[ibin][ibin2], wt);


      Multiplicity_vs_radius_ptwt[ibin][ibin2]->Sumw2();   more_hists->Multiplicity_vs_radius_ptwt[ibin][ibin2]->Sumw2();
      Multiplicity_vs_radius_ptwt[ibin][ibin2]->Add(more_hists->Multiplicity_vs_radius_ptwt[ibin][ibin2], wt);


      Multiplicity_vs_radius_bkg_ptwt[ibin][ibin2]->Sumw2();   more_hists->Multiplicity_vs_radius_bkg_ptwt[ibin][ibin2]->Sumw2();
      Multiplicity_vs_radius_bkg_ptwt[ibin][ibin2]->Add(more_hists->Multiplicity_vs_radius_bkg_ptwt[ibin][ibin2], wt);

      Multiplicity_vs_radius_gen_ptwt[ibin][ibin2]->Sumw2();   more_hists->Multiplicity_vs_radius_gen_ptwt[ibin][ibin2]->Sumw2();
      Multiplicity_vs_radius_gen_ptwt[ibin][ibin2]->Add(more_hists->Multiplicity_vs_radius_gen_ptwt[ibin][ibin2], wt);

      Multiplicity_vs_radius_gen_bkg_ptwt[ibin][ibin2]->Sumw2();   more_hists->Multiplicity_vs_radius_gen_bkg_ptwt[ibin][ibin2]->Sumw2();
      Multiplicity_vs_radius_gen_bkg_ptwt[ibin][ibin2]->Add(more_hists->Multiplicity_vs_radius_gen_bkg_ptwt[ibin][ibin2], wt);



      JetShapeDiffParticles_1D_check_error[ibin][ibin2]->Sumw2();   more_hists->JetShapeDiffParticles_1D_check_error[ibin][ibin2]->Sumw2();
      JetShapeDiffParticles_1D_check_error[ibin][ibin2]->Add(more_hists->JetShapeDiffParticles_1D_check_error[ibin][ibin2], wt);



      SumPt_only[ibin][ibin2]->Sumw2();   more_hists->SumPt_only[ibin][ibin2]->Sumw2();
      SumPt_only[ibin][ibin2]->Add(more_hists->SumPt_only[ibin][ibin2], wt);

      JetEnergy_resolution[ibin][ibin2]->Sumw2();   more_hists->JetEnergy_resolution[ibin][ibin2]->Sumw2();
      JetEnergy_resolution[ibin][ibin2]->Add(more_hists->JetEnergy_resolution[ibin][ibin2], wt);

      JetEnergy_ratio[ibin][ibin2]->Sumw2();   more_hists->JetEnergy_ratio[ibin][ibin2]->Sumw2();
      JetEnergy_ratio[ibin][ibin2]->Add(more_hists->JetEnergy_ratio[ibin][ibin2], wt);


      JetEnergy_ratio_new[ibin][ibin2]->Sumw2();   more_hists->JetEnergy_ratio_new[ibin][ibin2]->Sumw2();
      JetEnergy_ratio_new[ibin][ibin2]->Add(more_hists->JetEnergy_ratio_new[ibin][ibin2], wt);


      track_gen_pT_hist[ibin][ibin2]->Sumw2();   more_hists->track_gen_pT_hist[ibin][ibin2]->Sumw2();
      track_gen_pT_hist[ibin][ibin2]->Add(more_hists->track_gen_pT_hist[ibin][ibin2], wt);
      track_gen_eta_hist[ibin][ibin2]->Sumw2();   more_hists->track_gen_eta_hist[ibin][ibin2]->Sumw2();
      track_gen_eta_hist[ibin][ibin2]->Add(more_hists->track_gen_eta_hist[ibin][ibin2], wt);
      track_gen_phi_hist[ibin][ibin2]->Sumw2();   more_hists->track_gen_phi_hist[ibin][ibin2]->Sumw2();
      track_gen_phi_hist[ibin][ibin2]->Add(more_hists->track_gen_phi_hist[ibin][ibin2], wt);
      //
      track_gen_pT_pythia[ibin][ibin2]->Sumw2();   more_hists->track_gen_pT_pythia[ibin][ibin2]->Sumw2();
      track_gen_pT_pythia[ibin][ibin2]->Add(more_hists->track_gen_pT_pythia[ibin][ibin2], wt);
      track_gen_eta_pythia[ibin][ibin2]->Sumw2();   more_hists->track_gen_eta_pythia[ibin][ibin2]->Sumw2();
      track_gen_eta_pythia[ibin][ibin2]->Add(more_hists->track_gen_eta_pythia[ibin][ibin2], wt);
      track_gen_phi_pythia[ibin][ibin2]->Sumw2();   more_hists->track_gen_phi_pythia[ibin][ibin2]->Sumw2();
      track_gen_phi_pythia[ibin][ibin2]->Add(more_hists->track_gen_phi_pythia[ibin][ibin2], wt);

      track_gen_pT_hydjet[ibin][ibin2]->Sumw2();   more_hists->track_gen_pT_hydjet[ibin][ibin2]->Sumw2();
      track_gen_pT_hydjet[ibin][ibin2]->Add(more_hists->track_gen_pT_hydjet[ibin][ibin2], wt);
      track_gen_eta_hydjet[ibin][ibin2]->Sumw2();   more_hists->track_gen_eta_hydjet[ibin][ibin2]->Sumw2();
      track_gen_eta_hydjet[ibin][ibin2]->Add(more_hists->track_gen_eta_hydjet[ibin][ibin2], wt);
      track_gen_phi_hydjet[ibin][ibin2]->Sumw2();   more_hists->track_gen_phi_hydjet[ibin][ibin2]->Sumw2();
      track_gen_phi_hydjet[ibin][ibin2]->Add(more_hists->track_gen_phi_hydjet[ibin][ibin2], wt);

        
        
        track_gen_pT_hydjet2[ibin][ibin2]->Sumw2();   more_hists->track_gen_pT_hydjet2[ibin][ibin2]->Sumw2();
        track_gen_pT_hydjet2[ibin][ibin2]->Add(more_hists->track_gen_pT_hydjet2[ibin][ibin2], wt);
        
        track_gen_pT_hydjet3[ibin][ibin2]->Sumw2();   more_hists->track_gen_pT_hydjet3[ibin][ibin2]->Sumw2();
        track_gen_pT_hydjet3[ibin][ibin2]->Add(more_hists->track_gen_pT_hydjet3[ibin][ibin2], wt);
        
        track_gen_pT_hydjet4[ibin][ibin2]->Sumw2();   more_hists->track_gen_pT_hydjet4[ibin][ibin2]->Sumw2();
        track_gen_pT_hydjet4[ibin][ibin2]->Add(more_hists->track_gen_pT_hydjet4[ibin][ibin2], wt);
        
        track_gen_pT_hydjet5[ibin][ibin2]->Sumw2();   more_hists->track_gen_pT_hydjet5[ibin][ibin2]->Sumw2();
        track_gen_pT_hydjet5[ibin][ibin2]->Add(more_hists->track_gen_pT_hydjet5[ibin][ibin2], wt);
        
        track_gen_pT_hydjet6[ibin][ibin2]->Sumw2();   more_hists->track_gen_pT_hydjet6[ibin][ibin2]->Sumw2();
        track_gen_pT_hydjet6[ibin][ibin2]->Add(more_hists->track_gen_pT_hydjet6[ibin][ibin2], wt);
        
        track_gen_pT_hydjet7[ibin][ibin2]->Sumw2();   more_hists->track_gen_pT_hydjet7[ibin][ibin2]->Sumw2();
        track_gen_pT_hydjet7[ibin][ibin2]->Add(more_hists->track_gen_pT_hydjet7[ibin][ibin2], wt);
        
        track_gen_pT_hydjet8[ibin][ibin2]->Sumw2();   more_hists->track_gen_pT_hydjet8[ibin][ibin2]->Sumw2();
        track_gen_pT_hydjet8[ibin][ibin2]->Add(more_hists->track_gen_pT_hydjet8[ibin][ibin2], wt);
        
        track_gen_pT_hydjet9[ibin][ibin2]->Sumw2();   more_hists->track_gen_pT_hydjet9[ibin][ibin2]->Sumw2();
        track_gen_pT_hydjet9[ibin][ibin2]->Add(more_hists->track_gen_pT_hydjet9[ibin][ibin2], wt);
        
        track_gen_pT_hydjet10[ibin][ibin2]->Sumw2();   more_hists->track_gen_pT_hydjet10[ibin][ibin2]->Sumw2();
        track_gen_pT_hydjet10[ibin][ibin2]->Add(more_hists->track_gen_pT_hydjet10[ibin][ibin2], wt);
        
        track_gen_pT_hydjet11[ibin][ibin2]->Sumw2();   more_hists->track_gen_pT_hydjet11[ibin][ibin2]->Sumw2();
        track_gen_pT_hydjet11[ibin][ibin2]->Add(more_hists->track_gen_pT_hydjet11[ibin][ibin2], wt);

      JetShapeDiffParticles_1D_EMix[ibin][ibin2]->Sumw2();   more_hists->JetShapeDiffParticles_1D_EMix[ibin][ibin2]->Sumw2();
      JetShapeDiffParticles_1D_EMix[ibin][ibin2]->Add(more_hists->JetShapeDiffParticles_1D_EMix[ibin][ibin2], wt);



      emix_jetpt[ibin][ibin2]->Sumw2();   more_hists->emix_jetpt[ibin][ibin2]->Sumw2();
      emix_jetpt[ibin][ibin2]->Add(more_hists->emix_jetpt[ibin][ibin2], wt);

      emixed_jetpt[ibin][ibin2]->Sumw2();   more_hists->emixed_jetpt[ibin][ibin2]->Sumw2();
      emixed_jetpt[ibin][ibin2]->Add(more_hists->emixed_jetpt[ibin][ibin2], wt);

      emixed_jeteta[ibin][ibin2]->Sumw2();   more_hists->emixed_jeteta[ibin][ibin2]->Sumw2();
      emixed_jeteta[ibin][ibin2]->Add(more_hists->emixed_jeteta[ibin][ibin2], wt);


      emixed_jetphi[ibin][ibin2]->Sumw2();   more_hists->emixed_jetphi[ibin][ibin2]->Sumw2();
      emixed_jetphi[ibin][ibin2]->Add(more_hists->emixed_jetphi[ibin][ibin2], wt);


      for (int ibin3=0;ibin3<nTrkPtBins;ibin3++)
      {

	hJetTrackSignalBackground[ibin][ibin2][ibin3]->Sumw2();   more_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3]->Sumw2();
	hJetTrackSignalBackground[ibin][ibin2][ibin3]->Add(more_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3], wt);

	hJetTrackBackground[ibin][ibin2][ibin3]->Sumw2();   more_hists->hJetTrackBackground[ibin][ibin2][ibin3]->Sumw2();
	hJetTrackBackground[ibin][ibin2][ibin3]->Add(more_hists->hJetTrackBackground[ibin][ibin2][ibin3], wt);

  	Ngen[ibin][ibin2][ibin3]->Sumw2();   more_hists->Ngen[ibin][ibin2][ibin3]->Sumw2();
  	Ngen[ibin][ibin2][ibin3]->Add(more_hists->Ngen[ibin][ibin2][ibin3], wt);

  	Ntracks[ibin][ibin2][ibin3]->Sumw2();   more_hists->Ntracks[ibin][ibin2][ibin3]->Sumw2();
  	Ntracks[ibin][ibin2][ibin3]->Add(more_hists->Ntracks[ibin][ibin2][ibin3], wt);

  	Ntracks_awayside[ibin][ibin2][ibin3]->Sumw2();   more_hists->Ntracks_awayside[ibin][ibin2][ibin3]->Sumw2();
  	Ntracks_awayside[ibin][ibin2][ibin3]->Add(more_hists->Ntracks_awayside[ibin][ibin2][ibin3], wt);

	TrkPt[ibin][ibin2][ibin3]->Sumw2();   more_hists->TrkPt[ibin][ibin2][ibin3]->Sumw2();
	TrkPt[ibin][ibin2][ibin3]->Add(more_hists->TrkPt[ibin][ibin2][ibin3], wt);

	TrkEta[ibin][ibin2][ibin3]->Sumw2();   more_hists->TrkEta[ibin][ibin2][ibin3]->Sumw2();
	TrkEta[ibin][ibin2][ibin3]->Add(more_hists->TrkEta[ibin][ibin2][ibin3], wt);

	TrkPhi[ibin][ibin2][ibin3]->Sumw2();   more_hists->TrkPhi[ibin][ibin2][ibin3]->Sumw2();
	TrkPhi[ibin][ibin2][ibin3]->Add(more_hists->TrkPhi[ibin][ibin2][ibin3], wt);


	  TrkPt_weighted[ibin][ibin2][ibin3]->Sumw2();   more_hists->TrkPt_weighted[ibin][ibin2][ibin3]->Sumw2();
 	 TrkPt_weighted[ibin][ibin2][ibin3]->Add(more_hists->TrkPt_weighted[ibin][ibin2][ibin3], wt);

	  TrkEta_weighted[ibin][ibin2][ibin3]->Sumw2();   more_hists->TrkEta_weighted[ibin][ibin2][ibin3]->Sumw2();
  	TrkEta_weighted[ibin][ibin2][ibin3]->Add(more_hists->TrkEta_weighted[ibin][ibin2][ibin3], wt);

 	 TrkPhi_weighted[ibin][ibin2][ibin3]->Sumw2();   more_hists->TrkPhi_weighted[ibin][ibin2][ibin3]->Sumw2();
  	TrkPhi_weighted[ibin][ibin2][ibin3]->Add(more_hists->TrkPhi_weighted[ibin][ibin2][ibin3], wt);



	hJetGenTrackSignalBackground[ibin][ibin2][ibin3]->Sumw2();   more_hists->hJetGenTrackSignalBackground[ibin][ibin2][ibin3]->Sumw2();
	hJetGenTrackSignalBackground[ibin][ibin2][ibin3]->Add(more_hists->hJetGenTrackSignalBackground[ibin][ibin2][ibin3], wt);

	hJetGenTrackBackground[ibin][ibin2][ibin3]->Sumw2();   more_hists->hJetGenTrackBackground[ibin][ibin2][ibin3]->Sumw2();
	hJetGenTrackBackground[ibin][ibin2][ibin3]->Add(more_hists->hJetGenTrackBackground[ibin][ibin2][ibin3], wt);

	GenTrkEta[ibin][ibin2][ibin3]->Sumw2();   more_hists->GenTrkEta[ibin][ibin2][ibin3]->Sumw2();
	GenTrkEta[ibin][ibin2][ibin3]->Add(more_hists->GenTrkEta[ibin][ibin2][ibin3], wt);

	GenTrkPhi[ibin][ibin2][ibin3]->Sumw2();   more_hists->GenTrkPhi[ibin][ibin2][ibin3]->Sumw2();
	GenTrkPhi[ibin][ibin2][ibin3]->Add(more_hists->GenTrkPhi[ibin][ibin2][ibin3], wt);

  	GenTrkPt[ibin][ibin2][ibin3]->Sumw2();   more_hists->GenTrkPt[ibin][ibin2][ibin3]->Sumw2();
  	GenTrkPt[ibin][ibin2][ibin3]->Add(more_hists->GenTrkPt[ibin][ibin2][ibin3], wt);


  /// bkg for inc tracks

  GenTrkEta_bkg[ibin][ibin2][ibin3]->Sumw2();   more_hists->GenTrkEta_bkg[ibin][ibin2][ibin3]->Sumw2();
  GenTrkEta_bkg[ibin][ibin2][ibin3]->Add(more_hists->GenTrkEta_bkg[ibin][ibin2][ibin3], wt);

  GenTrkPhi_bkg[ibin][ibin2][ibin3]->Sumw2();   more_hists->GenTrkPhi_bkg[ibin][ibin2][ibin3]->Sumw2();
  GenTrkPhi_bkg[ibin][ibin2][ibin3]->Add(more_hists->GenTrkPhi_bkg[ibin][ibin2][ibin3], wt);

  TrkEta_bkg[ibin][ibin2][ibin3]->Sumw2();   more_hists->TrkEta_bkg[ibin][ibin2][ibin3]->Sumw2();
  TrkEta_bkg[ibin][ibin2][ibin3]->Add(more_hists->TrkEta_bkg[ibin][ibin2][ibin3], wt);

  TrkPhi_bkg[ibin][ibin2][ibin3]->Sumw2();   more_hists->TrkPhi_bkg[ibin][ibin2][ibin3]->Sumw2();
  TrkPhi_bkg[ibin][ibin2][ibin3]->Add(more_hists->TrkPhi_bkg[ibin][ibin2][ibin3], wt);



  hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Sumw2();   more_hists->hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Sumw2();
  hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Add(more_hists->hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3], wt);

  hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3]->Sumw2();   more_hists->hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3]->Sumw2();
  hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3]->Add(more_hists->hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3], wt);

  hJetGenTrackSignalBackground_pythia[ibin][ibin2][ibin3]->Sumw2();   more_hists->hJetGenTrackSignalBackground_pythia[ibin][ibin2][ibin3]->Sumw2();
  hJetGenTrackSignalBackground_pythia[ibin][ibin2][ibin3]->Add(more_hists->hJetGenTrackSignalBackground_pythia[ibin][ibin2][ibin3], wt);

  hJetGenTrackSignalBackground_hydjet[ibin][ibin2][ibin3]->Sumw2();   more_hists->hJetGenTrackSignalBackground_hydjet[ibin][ibin2][ibin3]->Sumw2();
  hJetGenTrackSignalBackground_hydjet[ibin][ibin2][ibin3]->Add(more_hists->hJetGenTrackSignalBackground_hydjet[ibin][ibin2][ibin3], wt);

  hJetGenTrackSignalBackground_hydro[ibin][ibin2][ibin3]->Sumw2();   more_hists->hJetGenTrackSignalBackground_hydro[ibin][ibin2][ibin3]->Sumw2();
  hJetGenTrackSignalBackground_hydro[ibin][ibin2][ibin3]->Add(more_hists->hJetGenTrackSignalBackground_hydro[ibin][ibin2][ibin3], wt);

  hJetGenTrackSignalBackground_minijets[ibin][ibin2][ibin3]->Sumw2();   more_hists->hJetGenTrackSignalBackground_minijets[ibin][ibin2][ibin3]->Sumw2();
  hJetGenTrackSignalBackground_minijets[ibin][ibin2][ibin3]->Add(more_hists->hJetGenTrackSignalBackground_minijets[ibin][ibin2][ibin3], wt);


  hJetGenTrackBackground_pythia[ibin][ibin2][ibin3]->Sumw2();   more_hists->hJetGenTrackBackground_pythia[ibin][ibin2][ibin3]->Sumw2();
  hJetGenTrackBackground_pythia[ibin][ibin2][ibin3]->Add(more_hists->hJetGenTrackBackground_pythia[ibin][ibin2][ibin3], wt);

  hJetGenTrackBackground_hydjet[ibin][ibin2][ibin3]->Sumw2();   more_hists->hJetGenTrackBackground_hydjet[ibin][ibin2][ibin3]->Sumw2();
  hJetGenTrackBackground_hydjet[ibin][ibin2][ibin3]->Add(more_hists->hJetGenTrackBackground_hydjet[ibin][ibin2][ibin3], wt);






      } /// ibin3



    }



  }
}


void hist_class::Delete()
{

  //  delete JetEnergy_gen_vs_rec;
  delete NumberOfMatches;
  delete NEvents;
  delete NEvents_test;
  delete NEvents_after_noise;
  delete NEvents_after_spike;
  delete NEvents_after_dphi;
  delete NEvents_before_dphi;
  delete NEvents_dijets;
  delete NEvents_after_trigger;
  delete alljets_rec;

  delete AllJetPt_raw_hist;
  delete AllJetPhi_hist;
  delete AllJetEta_hist;
  delete AllJetPt_hist;
  delete First_AllJetPhi_hist;
  delete First_AllJetEta_hist;
  delete First_AllJetPt_hist;
  delete Sub_AllJetPhi_hist;
  delete Sub_AllJetEta_hist;
  delete Sub_AllJetPt_hist;
  delete JetPt_fraction;
  delete Centrality;
  delete Centrality_weighted;

  delete track_vz;
  delete track_vz_weighted;

  for (int ibin2=0;ibin2<nPtBins;ibin2++){
    delete fullrecjet_aftercuts[ibin2];
    delete fullrefjet_aftercuts[ibin2];
  }



  for (int ibin=0;ibin<nCBins;ibin++){

    delete jet_pT_hist[ibin];
    delete jet_phi_hist[ibin];
    delete jet_eta_hist[ibin];
    delete jet_corrpT_hist[ibin];
    delete LeadingJetPt_hist[ibin];
    delete LeadingJetPhi_hist[ibin];
    delete LeadingJetEta_hist[ibin];
    delete SubJetPt_hist[ibin];
    delete SubJetPhi_hist[ibin];
    delete SubJetEta_hist[ibin];

    delete dPhi_leadingjet_hist[ibin];
    delete dPhi_subleadingjet_hist[ibin];

    delete JetShapeDiffParticles_bkg_1D_total[ibin];
    delete JetShapeIntegratedParticles_bkg_1D_total[ibin];

    delete JetEnergy_gen_vs_rec[ibin];


    for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){

      delete rr_gen_leadingjet[ibin][ibin3];
      delete rr_reco_leadingjet[ibin][ibin3];
      delete rr_reco_leadingjet_weighted[ibin][ibin3];
    }



    delete dPhi_hist[ibin];
    delete dPhi_after_hist[ibin];
    delete Aj[ibin];
    delete dPhi_vector[ibin];


    for (int ibin2=0;ibin2<nPtBins;ibin2++){


      delete ThirdJetPt_hist[ibin][ibin2];
      delete ThirdJetPhi_hist[ibin][ibin2];
      delete ThirdJetEta_hist[ibin][ibin2];
      //delete dPhi_hist[ibin][ibin2];
      delete dEta_hist[ibin][ibin2];
      //delete dPhi_after_hist[ibin][ibin2];

      delete EtaRef_bkg_pt[ibin][ibin2];
      delete EtaRef_bkg_pt_weighted[ibin][ibin2];



      delete all_jets_corrpT[ibin][ibin2];
      delete all_jets_phi[ibin][ibin2];
      delete all_jets_eta[ibin][ibin2];

      delete only_leadingjets_corrpT[ibin][ibin2];
      delete only_leadingjets_phi[ibin][ibin2];
      delete only_leadingjets_eta[ibin][ibin2];
      delete neutral_cand_pT_hist[ibin][ibin2];
      delete neutral_cand_phi_hist[ibin][ibin2];
      delete neutral_cand_eta_hist[ibin][ibin2];
      delete photons_cand_pT_hist[ibin][ibin2];
      delete photons_cand_phi_hist[ibin][ibin2];
      delete photons_cand_eta_hist[ibin][ibin2];
      delete NumNeutral[ibin][ibin2];
      delete NumPhotons[ibin][ibin2];
      delete NumAll[ibin][ibin2];
      delete NumCharged[ibin][ibin2];

      delete NumAll_bkg[ibin][ibin2];
      delete NumNeutral_bkg[ibin][ibin2];
      delete NumPhotons_bkg[ibin][ibin2];
      delete NumChargedHadrons_bkg[ibin][ibin2];
      delete NumChargedParticles_bkg[ibin][ibin2];
      delete NumElectrons_bkg[ibin][ibin2];
      delete NumMuons_bkg[ibin][ibin2];

      delete radius_hist[ibin][ibin2];
      delete radius_hist_mine[ibin][ibin2];
      delete radius_hist_mine2[ibin][ibin2];
      delete radius_hist_mine3[ibin][ibin2];

      delete JetShapeIntegratedParticles[ibin][ibin2];
      delete JetShapeDiffParticles[ibin][ibin2];
      delete JetShapeIntegratedParticles_bkgsub[ibin][ibin2];
      delete JetShapeDiffParticles_bkgsub[ibin][ibin2];
      delete JetShapeIntegratedParticles_bkg[ibin][ibin2];
      delete JetShapeDiffParticles_bkg[ibin][ibin2];



      delete SumJetPtFraction_hist[ibin][ibin2];
      delete all_cand_pT_hist[ibin][ibin2];
      delete all_cand_phi_hist[ibin][ibin2];
      delete all_cand_eta_hist[ibin][ibin2];

      delete track_cand_pT_hist[ibin][ibin2];
      delete track_cand_phi_hist[ibin][ibin2];
      delete track_cand_eta_hist[ibin][ibin2];
      delete track_cand_pT_hist_subleadingjet[ibin][ibin2];

      delete only_subleadingjets_corrpT[ibin][ibin2];
      delete only_subleadingjets_phi[ibin][ibin2];
      delete only_subleadingjets_eta[ibin][ibin2];
      delete neutral_cand_pT_hist_subleadingjet[ibin][ibin2];
      delete neutral_cand_phi_hist_subleadingjet[ibin][ibin2];
      delete neutral_cand_eta_hist_subleadingjet[ibin][ibin2];
      delete photons_cand_pT_hist_subleadingjet[ibin][ibin2];
      delete photons_cand_phi_hist_subleadingjet[ibin][ibin2];
      delete photons_cand_eta_hist_subleadingjet[ibin][ibin2];
      delete NumNeutral_subleadingjet[ibin][ibin2];
      delete NumPhotons_subleadingjet[ibin][ibin2];
      delete NumAll_subleadingjet[ibin][ibin2];
      delete NumChargedHadrons_subleadingjet[ibin][ibin2];



      delete JetShapeDiffParticles_bkg_1D[ibin][ibin2];
      delete JetShapeDiffParticles_1D[ibin][ibin2];

      delete JetShapeIntegratedParticles_bkg_1D[ibin][ibin2];
      delete JetShapeIntegratedParticles_1D[ibin][ibin2];

      delete JetShapeDiffParticlesGen_bkg_1D[ibin][ibin2];
      delete JetShapeDiffParticlesGen_1D[ibin][ibin2];



      delete SumJetPtFraction_hist_subleadingjet[ibin][ibin2];
      delete SumJetPtFraction_hist_subleadingjet_photons[ibin][ibin2];
      delete SumJetPtFraction_hist_subleadingjet_neutrals[ibin][ibin2];
      delete SumJetPtFraction_hist_subleadingjet_chargedhadrons[ibin][ibin2];

      delete SumJetPtFraction_hist_leadingjet_photons[ibin][ibin2];
      delete SumJetPtFraction_hist_leadingjet_neutrals[ibin][ibin2];
      delete SumJetPtFraction_hist_leadingjet_charged[ibin][ibin2];


      delete all_cand_pT_hist_subleadingjet[ibin][ibin2];
      delete all_cand_phi_hist_subleadingjet[ibin][ibin2];
      delete all_cand_eta_hist_subleadingjet[ibin][ibin2];
      delete radius_hist_subleadingjet[ibin][ibin2];
      delete Centrality_hist[ibin][ibin2];
      //delete Aj[ibin][ibin2];

      delete SumJetPtFraction_hist_leadingjet_chargedhadrons[ibin][ibin2];
      delete SumJetPtFraction_hist_leadingjet_chargedparticles[ibin][ibin2];
      delete SumJetPtFraction_hist_leadingjet_electrons[ibin][ibin2];
      delete SumJetPtFraction_hist_leadingjet_muons[ibin][ibin2];

      delete NumChargedHadrons[ibin][ibin2];
      delete NumChargedParticles[ibin][ibin2];
      delete NumElectrons[ibin][ibin2];
      delete NumMuons[ibin][ibin2];

      delete chargedhadrons_cand_pT_hist[ibin][ibin2];
      delete chargedhadrons_cand_phi_hist[ibin][ibin2];
      delete chargedhadrons_cand_eta_hist[ibin][ibin2];

      delete chargedparticles_cand_pT_hist[ibin][ibin2];
      delete chargedparticles_cand_phi_hist[ibin][ibin2];
      delete chargedparticles_cand_eta_hist[ibin][ibin2];

      delete electrons_cand_pT_hist[ibin][ibin2];
      delete electrons_cand_phi_hist[ibin][ibin2];
      delete electrons_cand_eta_hist[ibin][ibin2];

      delete muons_cand_pT_hist[ibin][ibin2];
      delete muons_cand_phi_hist[ibin][ibin2];
      delete muons_cand_eta_hist[ibin][ibin2];

      delete dN_tracks[ibin][ibin2];
      delete dN_chargedhadrons[ibin][ibin2];
      delete dN_chargedparticles[ibin][ibin2];
      delete dN_electrons[ibin][ibin2];
      delete dN_muons[ibin][ibin2];
      delete dN_neutrals[ibin][ibin2];
      delete dN_photons[ibin][ibin2];

      delete chargedhadrons_cand_pT_hist_bkg[ibin][ibin2];
      delete chargedparticles_cand_pT_hist_bkg[ibin][ibin2];
      delete electrons_cand_pT_hist_bkg[ibin][ibin2];
      delete muons_cand_pT_hist_bkg[ibin][ibin2];
      delete neutral_cand_pT_hist_bkg[ibin][ibin2];
      delete photons_cand_pT_hist_bkg[ibin][ibin2];

      delete SumJetPtFraction_hist_leadingjet_photons_bkg[ibin][ibin2];
      delete SumJetPtFraction_hist_leadingjet_neutrals_bkg[ibin][ibin2];
      delete SumJetPtFraction_hist_leadingjet_chargedhadrons_bkg[ibin][ibin2];
      delete SumJetPtFraction_hist_leadingjet_chargedparticles_bkg[ibin][ibin2];
      delete SumJetPtFraction_hist_leadingjet_electrons_bkg[ibin][ibin2];
      delete SumJetPtFraction_hist_leadingjet_muons_bkg[ibin][ibin2];


      /// subleading jet
      delete chargedhadrons_cand_pT_hist_subleadingjet[ibin][ibin2];
      delete chargedhadrons_cand_phi_hist_subleadingjet[ibin][ibin2];
      delete chargedhadrons_cand_eta_hist_subleadingjet[ibin][ibin2];

      delete chargedhadrons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2];
      delete neutral_cand_pT_hist_subleadingjet_bkg[ibin][ibin2];
      delete photons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2];


      delete SumJetPtFraction_hist_subleadingjet_photons_bkg[ibin][ibin2];
      delete SumJetPtFraction_hist_subleadingjet_neutrals_bkg[ibin][ibin2];
      delete SumJetPtFraction_hist_subleadingjet_chargedhadrons_bkg[ibin][ibin2];
      delete SumJetPtFraction_hist_subleadingjet_chargedparticles_bkg[ibin][ibin2];
      delete SumJetPtFraction_hist_subleadingjet_electrons_bkg[ibin][ibin2];
      delete SumJetPtFraction_hist_subleadingjet_muons_bkg[ibin][ibin2];


      delete NumAll_subleadingjet_bkg[ibin][ibin2];
      delete NumNeutral_subleadingjet_bkg[ibin][ibin2];
      delete NumPhotons_subleadingjet_bkg[ibin][ibin2];
      delete NumChargedHadrons_subleadingjet_bkg[ibin][ibin2];

      delete JetPt_fraction_hist[ibin][ibin2];
      delete dPhi_jet_track[ibin][ibin2];
      delete dPhi_jet_track_ptweight[ibin][ibin2];
      delete track_cand_pT_hist_weighted[ibin][ibin2];
      delete track_cand_eta_hist_weighted[ibin][ibin2];
      delete track_cand_phi_hist_weighted[ibin][ibin2];

      delete track_gen_pT_hist[ibin][ibin2];
      delete track_gen_eta_hist[ibin][ibin2];
      delete track_gen_phi_hist[ibin][ibin2];

      delete track_gen_pT_pythia[ibin][ibin2];
      delete track_gen_eta_pythia[ibin][ibin2];
      delete track_gen_phi_pythia[ibin][ibin2];
      
      delete track_gen_pT_hydjet[ibin][ibin2];
      delete track_gen_eta_hydjet[ibin][ibin2];
      delete track_gen_phi_hydjet[ibin][ibin2];

        delete track_gen_pT_hydjet2[ibin][ibin2];
        delete track_gen_pT_hydjet3[ibin][ibin2];
        delete track_gen_pT_hydjet4[ibin][ibin2];
        delete track_gen_pT_hydjet5[ibin][ibin2];
        delete track_gen_pT_hydjet6[ibin][ibin2];
        delete track_gen_pT_hydjet7[ibin][ibin2];
        delete track_gen_pT_hydjet8[ibin][ibin2];
        delete track_gen_pT_hydjet9[ibin][ibin2];
        delete track_gen_pT_hydjet10[ibin][ibin2];
        delete track_gen_pT_hydjet11[ibin][ibin2];


      delete SumPt_only[ibin][ibin2];

      delete JetEnergy_resolution[ibin][ibin2];
      delete JetEnergy_ratio[ibin][ibin2];

      delete JetEnergy_ratio_new[ibin][ibin2];

      delete track_bkg_pT_hist[ibin][ibin2];
      delete track_bkg_pT_hist_weighted[ibin][ibin2];
      delete track_gen_bkg_pT_hist[ibin][ibin2];

      delete track_bkg_pT_hist_oocone[ibin][ibin2];
      delete track_bkg_pT_hist_weighted_oocone[ibin][ibin2];
      delete track_gen_bkg_pT_hist_oocone[ibin][ibin2];

      delete genparticles_pT_hist[ibin][ibin2];
      delete genparticles_bkg_pT_hist[ibin][ibin2];
      delete gensube_hist[ibin][ibin2];

      delete track_bkg_eta_hist[ibin][ibin2];
      delete track_bkg_eta_hist_weighted[ibin][ibin2];
      delete track_gen_bkg_eta_hist[ibin][ibin2];



      delete Multiplicity_vs_radius[ibin][ibin2];
      delete Multiplicity_vs_radius_bkg[ibin][ibin2];

      delete Multiplicity_vs_radius_gen[ibin][ibin2];
      delete Multiplicity_vs_radius_gen_bkg[ibin][ibin2];

      delete Multiplicity_vs_radius_gen_all[ibin][ibin2];
      delete Multiplicity_vs_radius_gen_all_bkg[ibin][ibin2];


      delete Multiplicity_vs_radius_ptwt[ibin][ibin2];
      delete Multiplicity_vs_radius_bkg_ptwt[ibin][ibin2];

      delete Multiplicity_vs_radius_gen_ptwt[ibin][ibin2];
      delete Multiplicity_vs_radius_gen_bkg_ptwt[ibin][ibin2];
      delete JetShapeDiffParticles_1D_check_error[ibin][ibin2];
      delete JetShapeDiffParticles_1D_EMix[ibin][ibin2];


      delete emix_jetpt[ibin][ibin2];
      delete emixed_jetpt[ibin][ibin2];
      delete emixed_jeteta[ibin][ibin2];
      delete emixed_jetphi[ibin][ibin2];


     for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){
	delete hJetTrackSignalBackground[ibin][ibin2][ibin3];
	delete hJetTrackBackground[ibin][ibin2][ibin3];
	delete hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3];
	delete hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3];

	delete Ntracks[ibin][ibin2][ibin3];
  delete Ngen[ibin][ibin2][ibin3];
  delete Ntracks_awayside[ibin][ibin2][ibin3];

	delete TrkPt[ibin][ibin2][ibin3];
	delete TrkEta[ibin][ibin2][ibin3];
	delete TrkPhi[ibin][ibin2][ibin3];

  delete TrkPt_weighted[ibin][ibin2][ibin3];
  delete TrkEta_weighted[ibin][ibin2][ibin3];
  delete TrkPhi_weighted[ibin][ibin2][ibin3];


	delete GenTrkPt[ibin][ibin2][ibin3];
  delete GenTrkEta[ibin][ibin2][ibin3];
  delete GenTrkPhi[ibin][ibin2][ibin3];


  delete GenTrkEta_bkg[ibin][ibin2][ibin3];
  delete GenTrkPhi_bkg[ibin][ibin2][ibin3];
  delete TrkEta_bkg[ibin][ibin2][ibin3];
  delete TrkPhi_bkg[ibin][ibin2][ibin3];


  delete hJetGenTrackSignalBackground[ibin][ibin2][ibin3];
  delete hJetGenTrackBackground[ibin][ibin2][ibin3];

  delete hJetGenTrackSignalBackground_pythia[ibin][ibin2][ibin3];
  delete hJetGenTrackSignalBackground_hydjet[ibin][ibin2][ibin3];

  delete hJetGenTrackSignalBackground_hydro[ibin][ibin2][ibin3];
  delete hJetGenTrackSignalBackground_minijets[ibin][ibin2][ibin3];

  delete hJetGenTrackBackground_pythia[ibin][ibin2][ibin3];
  delete hJetGenTrackBackground_hydjet[ibin][ibin2][ibin3];


      } /// ibin3


    } // ibin2


  } // ibin
}

void hist_class::Write()
{

  TString parti_str = "";
  if( parti >= 0 ) 
  {
    parti_str += "_part";
    parti_str +=  parti;
  }

  TString pT_str = "";
  if( trkPtCut >= 0.49 && trkPtCut < 1.5 ) pT_str = "trkPtCut1";
  else if( trkPtCut >= 1.5 && trkPtCut < 2.5 ) pT_str = "trkPtCut2";
  else if( trkPtCut >= 2.5 && trkPtCut < 3.5 ) pT_str = "trkPtCut3";
  else if( trkPtCut >= 3.5 && trkPtCut < 4.5 ) pT_str = "trkPtCut4";
  else assert(0);  
 // TString out_name = (TString) ("root_output/" + dataset_type_strs[dataset_type_code] + "_" + pT_str + parti_str + ".root");
 // if( test_it ) out_name = (TString) ("root_output/" + dataset_type_strs[dataset_type_code] + "_" + pT_str + parti_str + "_test.root");
    
    
    TString out_name = (TString) ("root_output/" + dataset_type_strs[dataset_type_code] + "_" + npart_type_strs[npart] + "_" + pT_str + ".root");
    if( test_it ) out_name = (TString) ("root_output/" + dataset_type_strs[dataset_type_code] + "_" + npart_type_strs[npart] + "_" + pT_str + "_test.root");

    
  TFile *out_file = new TFile(out_name, "RECREATE");
  //std::cout << "parti: " << parti << ", parti_str: " << parti_str << ", out_name: " << out_name << std::endl;


  //     JetEnergy_gen_vs_rec->Write();

  alljets_rec->Write();
  NumberOfMatches->Write();
  NEvents->Write();
  NEvents_test->Write();
  NEvents_after_noise->Write();
  NEvents_after_spike->Write();
  NEvents_after_trigger->Write();
  NEvents_after_dphi->Write();
  NEvents_before_dphi->Write();
  NEvents_dijets->Write();

  AllJetPt_raw_hist->Write();
  AllJetPhi_hist->Write();
  AllJetEta_hist->Write();
  AllJetPt_hist->Write();

  First_AllJetPhi_hist->Write();
  First_AllJetEta_hist->Write();
  First_AllJetPt_hist->Write();

  Sub_AllJetPhi_hist->Write();
  Sub_AllJetEta_hist->Write();
  Sub_AllJetPt_hist->Write();
  JetPt_fraction->Write();
  Centrality->Write();
  Centrality_weighted->Write();

  track_vz->Write();
  track_vz_weighted->Write();

  for (int ibin2=0;ibin2<nPtBins;ibin2++)
  {
    fullrecjet_aftercuts[ibin2]->Write();
    fullrefjet_aftercuts[ibin2]->Write();
  }


  for (int ibin=0;ibin<nCBins;ibin++)
  {

    jet_pT_hist[ibin]->Write();
    jet_phi_hist[ibin]->Write();
    jet_eta_hist[ibin]->Write();
    jet_corrpT_hist[ibin]->Write();
    LeadingJetPt_hist[ibin]->Write();
    LeadingJetPhi_hist[ibin]->Write();
    LeadingJetEta_hist[ibin]->Write();
    SubJetPt_hist[ibin]->Write();
    SubJetPhi_hist[ibin]->Write();
    SubJetEta_hist[ibin]->Write();
    dPhi_leadingjet_hist[ibin]->Write();
    dPhi_subleadingjet_hist[ibin]->Write();
    JetShapeDiffParticles_bkg_1D_total[ibin]->Write();
    JetShapeIntegratedParticles_bkg_1D_total[ibin]->Write();
    JetEnergy_gen_vs_rec[ibin]->Write();


 	   for (int ibin3=0;ibin3<nTrkPtBins;ibin3++)
    	{
      	rr_gen_leadingjet[ibin][ibin3]->Write();
      	rr_reco_leadingjet[ibin][ibin3]->Write();
      	rr_reco_leadingjet_weighted[ibin][ibin3]->Write();
    	}


    dPhi_hist[ibin]->Write();
    dPhi_after_hist[ibin]->Write();
    Aj[ibin]->Write();
    dPhi_vector[ibin]->Write();

	for (int ibin2=0;ibin2<nPtBins;ibin2++)
	{
        EtaRef_bkg_pt[ibin][ibin2]->Write();
        EtaRef_bkg_pt_weighted[ibin][ibin2]->Write();

        ThirdJetPt_hist[ibin][ibin2]->Write();
        ThirdJetPhi_hist[ibin][ibin2]->Write();
        ThirdJetEta_hist[ibin][ibin2]->Write();

      	radius_hist[ibin][ibin2]->Write();
      	radius_hist_mine[ibin][ibin2]->Write();
      	radius_hist_mine2[ibin][ibin2]->Write();
      	radius_hist_mine3[ibin][ibin2]->Write();


      //dPhi_hist[ibin][ibin2]->Write();
      	dEta_hist[ibin][ibin2]->Write();
      //dPhi_after_hist[ibin][ibin2]->Write();


      	all_jets_corrpT[ibin][ibin2]->Write();
      	all_jets_phi[ibin][ibin2]->Write();
      	all_jets_eta[ibin][ibin2]->Write();
      	only_leadingjets_corrpT[ibin][ibin2]->Write();
      	only_leadingjets_phi[ibin][ibin2]->Write();
      	only_leadingjets_eta[ibin][ibin2]->Write();
      	only_subleadingjets_corrpT[ibin][ibin2]->Write();
      	only_subleadingjets_phi[ibin][ibin2]->Write();
      	only_subleadingjets_eta[ibin][ibin2]->Write();


      	neutral_cand_pT_hist[ibin][ibin2]->Write();
      	neutral_cand_phi_hist[ibin][ibin2]->Write();
      	neutral_cand_eta_hist[ibin][ibin2]->Write();
      photons_cand_pT_hist[ibin][ibin2]->Write();
      photons_cand_phi_hist[ibin][ibin2]->Write();
      photons_cand_eta_hist[ibin][ibin2]->Write();
      NumNeutral[ibin][ibin2]->Write();
      NumPhotons[ibin][ibin2]->Write();
      NumAll[ibin][ibin2]->Write();
      NumCharged[ibin][ibin2]->Write();

      NumAll_bkg[ibin][ibin2]->Write();
      NumNeutral_bkg[ibin][ibin2]->Write();
      NumPhotons_bkg[ibin][ibin2]->Write();
      NumChargedHadrons_bkg[ibin][ibin2]->Write();
      NumChargedParticles_bkg[ibin][ibin2]->Write();
      NumElectrons_bkg[ibin][ibin2]->Write();
      NumMuons_bkg[ibin][ibin2]->Write();

      JetShapeIntegratedParticles[ibin][ibin2]->Write();
      JetShapeDiffParticles[ibin][ibin2]->Write();
      JetShapeIntegratedParticles_bkgsub[ibin][ibin2]->Write();
      JetShapeDiffParticles_bkgsub[ibin][ibin2]->Write();
      JetShapeIntegratedParticles_bkg[ibin][ibin2]->Write();
      JetShapeDiffParticles_bkg[ibin][ibin2]->Write();


      JetShapeDiffParticles_bkg_1D[ibin][ibin2]->Write();
      JetShapeDiffParticles_1D[ibin][ibin2]->Write();
      JetShapeIntegratedParticles_bkg_1D[ibin][ibin2]->Write();
      JetShapeIntegratedParticles_1D[ibin][ibin2]->Write();

      JetShapeDiffParticlesGen_bkg_1D[ibin][ibin2]->Write();
      JetShapeDiffParticlesGen_1D[ibin][ibin2]->Write();

      SumJetPtFraction_hist[ibin][ibin2]->Write();
      all_cand_pT_hist[ibin][ibin2]->Write();
      all_cand_phi_hist[ibin][ibin2]->Write();
      all_cand_eta_hist[ibin][ibin2]->Write();

      track_cand_pT_hist[ibin][ibin2]->Write();
      track_cand_phi_hist[ibin][ibin2]->Write();
      track_cand_eta_hist[ibin][ibin2]->Write();
      track_cand_pT_hist_subleadingjet[ibin][ibin2]->Write();

      neutral_cand_pT_hist_subleadingjet[ibin][ibin2]->Write();
      neutral_cand_phi_hist_subleadingjet[ibin][ibin2]->Write();
      neutral_cand_eta_hist_subleadingjet[ibin][ibin2]->Write();
      photons_cand_pT_hist_subleadingjet[ibin][ibin2]->Write();
      photons_cand_phi_hist_subleadingjet[ibin][ibin2]->Write();
      photons_cand_eta_hist_subleadingjet[ibin][ibin2]->Write();
      NumNeutral_subleadingjet[ibin][ibin2]->Write();
      NumPhotons_subleadingjet[ibin][ibin2]->Write();
      NumAll_subleadingjet[ibin][ibin2]->Write();
      NumChargedHadrons_subleadingjet[ibin][ibin2]->Write();



      SumJetPtFraction_hist_subleadingjet[ibin][ibin2]->Write();
      SumJetPtFraction_hist_subleadingjet_photons[ibin][ibin2]->Write();
      SumJetPtFraction_hist_subleadingjet_neutrals[ibin][ibin2]->Write();
      SumJetPtFraction_hist_subleadingjet_chargedhadrons[ibin][ibin2]->Write();

      SumJetPtFraction_hist_leadingjet_photons[ibin][ibin2]->Write();
      SumJetPtFraction_hist_leadingjet_neutrals[ibin][ibin2]->Write();
      SumJetPtFraction_hist_leadingjet_charged[ibin][ibin2]->Write();

      SumJetPtFraction_hist_leadingjet_chargedhadrons[ibin][ibin2]->Write();
      SumJetPtFraction_hist_leadingjet_chargedparticles[ibin][ibin2]->Write();
      SumJetPtFraction_hist_leadingjet_electrons[ibin][ibin2]->Write();
      SumJetPtFraction_hist_leadingjet_muons[ibin][ibin2]->Write();

      NumChargedHadrons[ibin][ibin2]->Write();
      NumChargedParticles[ibin][ibin2]->Write();
      NumElectrons[ibin][ibin2]->Write();
      NumMuons[ibin][ibin2]->Write();

      chargedhadrons_cand_pT_hist[ibin][ibin2]->Write();
      chargedhadrons_cand_phi_hist[ibin][ibin2]->Write();
      chargedhadrons_cand_eta_hist[ibin][ibin2]->Write();

      chargedparticles_cand_pT_hist[ibin][ibin2]->Write();
      chargedparticles_cand_phi_hist[ibin][ibin2]->Write();
      chargedparticles_cand_eta_hist[ibin][ibin2]->Write();

      electrons_cand_pT_hist[ibin][ibin2]->Write();
      electrons_cand_phi_hist[ibin][ibin2]->Write();
      electrons_cand_eta_hist[ibin][ibin2]->Write();

      muons_cand_pT_hist[ibin][ibin2]->Write();
      muons_cand_phi_hist[ibin][ibin2]->Write();
      muons_cand_eta_hist[ibin][ibin2]->Write();


      all_cand_pT_hist_subleadingjet[ibin][ibin2]->Write();
      all_cand_phi_hist_subleadingjet[ibin][ibin2]->Write();
      all_cand_eta_hist_subleadingjet[ibin][ibin2]->Write();
      radius_hist_subleadingjet[ibin][ibin2]->Write();
      Centrality_hist[ibin][ibin2]->Write();
     /// Aj[ibin][ibin2]->Write();


      dN_tracks[ibin][ibin2]->Write();
      dN_chargedhadrons[ibin][ibin2]->Write();
      dN_chargedparticles[ibin][ibin2]->Write();
      dN_electrons[ibin][ibin2]->Write();
      dN_muons[ibin][ibin2]->Write();
      dN_neutrals[ibin][ibin2]->Write();
      dN_photons[ibin][ibin2]->Write();

      chargedhadrons_cand_pT_hist_bkg[ibin][ibin2]->Write();
      chargedparticles_cand_pT_hist_bkg[ibin][ibin2]->Write();
      electrons_cand_pT_hist_bkg[ibin][ibin2]->Write();
      muons_cand_pT_hist_bkg[ibin][ibin2]->Write();
      neutral_cand_pT_hist_bkg[ibin][ibin2]->Write();
      photons_cand_pT_hist_bkg[ibin][ibin2]->Write();

      SumJetPtFraction_hist_leadingjet_photons_bkg[ibin][ibin2]->Write();
      SumJetPtFraction_hist_leadingjet_neutrals_bkg[ibin][ibin2]->Write();
      SumJetPtFraction_hist_leadingjet_chargedhadrons_bkg[ibin][ibin2]->Write();
      SumJetPtFraction_hist_leadingjet_chargedparticles_bkg[ibin][ibin2]->Write();
      SumJetPtFraction_hist_leadingjet_electrons_bkg[ibin][ibin2]->Write();
      SumJetPtFraction_hist_leadingjet_muons_bkg[ibin][ibin2]->Write();


      /// subleading jet plots
      chargedhadrons_cand_pT_hist_subleadingjet[ibin][ibin2]->Write();
      chargedhadrons_cand_phi_hist_subleadingjet[ibin][ibin2]->Write();
      chargedhadrons_cand_eta_hist_subleadingjet[ibin][ibin2]->Write();

      chargedhadrons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Write();
      neutral_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Write();
      photons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Write();


      SumJetPtFraction_hist_subleadingjet_photons_bkg[ibin][ibin2]->Write();
      SumJetPtFraction_hist_subleadingjet_neutrals_bkg[ibin][ibin2]->Write();
      SumJetPtFraction_hist_subleadingjet_chargedhadrons_bkg[ibin][ibin2]->Write();
      SumJetPtFraction_hist_subleadingjet_chargedparticles_bkg[ibin][ibin2]->Write();
      SumJetPtFraction_hist_subleadingjet_electrons_bkg[ibin][ibin2]->Write();
      SumJetPtFraction_hist_subleadingjet_muons_bkg[ibin][ibin2]->Write();


      NumAll_subleadingjet_bkg[ibin][ibin2]->Write();
      NumNeutral_subleadingjet_bkg[ibin][ibin2]->Write();
      NumPhotons_subleadingjet_bkg[ibin][ibin2]->Write();
      NumChargedHadrons_subleadingjet_bkg[ibin][ibin2]->Write();

      JetPt_fraction_hist[ibin][ibin2]->Write();
      dPhi_jet_track[ibin][ibin2]->Write();
      dPhi_jet_track_ptweight[ibin][ibin2]->Write();
      
      track_cand_pT_hist_weighted[ibin][ibin2]->Write();
      track_cand_eta_hist_weighted[ibin][ibin2]->Write();
      track_cand_phi_hist_weighted[ibin][ibin2]->Write();

      track_gen_pT_hist[ibin][ibin2]->Write();
      track_gen_eta_hist[ibin][ibin2]->Write();
      track_gen_phi_hist[ibin][ibin2]->Write();
      track_gen_bkg_pT_hist[ibin][ibin2]->Write();
      track_gen_bkg_pT_hist_oocone[ibin][ibin2]->Write();

      track_gen_pT_pythia[ibin][ibin2]->Write();
      track_gen_eta_pythia[ibin][ibin2]->Write();
      track_gen_phi_pythia[ibin][ibin2]->Write();

      track_gen_pT_hydjet[ibin][ibin2]->Write();
      track_gen_eta_hydjet[ibin][ibin2]->Write();
      track_gen_phi_hydjet[ibin][ibin2]->Write();


      track_gen_pT_hydjet2[ibin][ibin2]->Write();
      track_gen_pT_hydjet3[ibin][ibin2]->Write();
      track_gen_pT_hydjet4[ibin][ibin2]->Write();
      track_gen_pT_hydjet5[ibin][ibin2]->Write();
      track_gen_pT_hydjet6[ibin][ibin2]->Write();
      track_gen_pT_hydjet7[ibin][ibin2]->Write();
      track_gen_pT_hydjet8[ibin][ibin2]->Write();
      track_gen_pT_hydjet9[ibin][ibin2]->Write();
      track_gen_pT_hydjet10[ibin][ibin2]->Write();
      track_gen_pT_hydjet11[ibin][ibin2]->Write();

      genparticles_pT_hist[ibin][ibin2]->Write();
      genparticles_bkg_pT_hist[ibin][ibin2]->Write();
      gensube_hist[ibin][ibin2]->Write();

      track_bkg_eta_hist[ibin][ibin2]->Write();
      track_bkg_eta_hist_weighted[ibin][ibin2]->Write();
      track_gen_bkg_eta_hist[ibin][ibin2]->Write();


      SumPt_only[ibin][ibin2]->Write();

      JetEnergy_resolution[ibin][ibin2]->Write();
      JetEnergy_ratio[ibin][ibin2]->Write();
      JetEnergy_ratio_new[ibin][ibin2]->Write();

      track_bkg_pT_hist[ibin][ibin2]->Write();
      track_bkg_pT_hist_weighted[ibin][ibin2]->Write();

      track_bkg_pT_hist_oocone[ibin][ibin2]->Write();
      track_bkg_pT_hist_weighted_oocone[ibin][ibin2]->Write();


      Multiplicity_vs_radius[ibin][ibin2]->Write();
      Multiplicity_vs_radius_bkg[ibin][ibin2]->Write();

      Multiplicity_vs_radius_gen[ibin][ibin2]->Write();
      Multiplicity_vs_radius_gen_bkg[ibin][ibin2]->Write();

      Multiplicity_vs_radius_gen_all[ibin][ibin2]->Write();
      Multiplicity_vs_radius_gen_all_bkg[ibin][ibin2]->Write();

      Multiplicity_vs_radius_ptwt[ibin][ibin2]->Write();
      Multiplicity_vs_radius_bkg_ptwt[ibin][ibin2]->Write();

      Multiplicity_vs_radius_gen_ptwt[ibin][ibin2]->Write();
      Multiplicity_vs_radius_gen_bkg_ptwt[ibin][ibin2]->Write();
      JetShapeDiffParticles_1D_check_error[ibin][ibin2]->Write();
      JetShapeDiffParticles_1D_EMix[ibin][ibin2]->Write();
      emix_jetpt[ibin][ibin2]->Write();
      emixed_jetpt[ibin][ibin2]->Write();
      emixed_jeteta[ibin][ibin2]->Write();
      emixed_jetphi[ibin][ibin2]->Write();


      for (int ibin3=0;ibin3<nTrkPtBins;ibin3++)
      {

	hJetTrackSignalBackground[ibin][ibin2][ibin3]->Write();
	hJetTrackBackground[ibin][ibin2][ibin3]->Write();

	hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Write();
	hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3]->Write();

	Ntracks[ibin][ibin2][ibin3]->Write();
  	Ngen[ibin][ibin2][ibin3]->Write();
  	Ntracks_awayside[ibin][ibin2][ibin3]->Write();

	TrkPt[ibin][ibin2][ibin3]->Write();
	TrkEta[ibin][ibin2][ibin3]->Write();
  	TrkPhi[ibin][ibin2][ibin3]->Write();

  	TrkPt_weighted[ibin][ibin2][ibin3]->Write();
  	TrkEta_weighted[ibin][ibin2][ibin3]->Write();
  	TrkPhi_weighted[ibin][ibin2][ibin3]->Write();


	GenTrkPt[ibin][ibin2][ibin3]->Write();
 	GenTrkEta[ibin][ibin2][ibin3]->Write();
 	GenTrkPhi[ibin][ibin2][ibin3]->Write();

  	GenTrkEta_bkg[ibin][ibin2][ibin3]->Write();
  	GenTrkPhi_bkg[ibin][ibin2][ibin3]->Write();
  	TrkEta_bkg[ibin][ibin2][ibin3]->Write();
  	TrkPhi_bkg[ibin][ibin2][ibin3]->Write();

	
      hJetGenTrackSignalBackground[ibin][ibin2][ibin3]->Write();
      hJetGenTrackBackground[ibin][ibin2][ibin3]->Write();

      hJetGenTrackSignalBackground_pythia[ibin][ibin2][ibin3]->Write();
      hJetGenTrackSignalBackground_hydjet[ibin][ibin2][ibin3]->Write();
      hJetGenTrackSignalBackground_hydro[ibin][ibin2][ibin3]->Write();
      hJetGenTrackSignalBackground_minijets[ibin][ibin2][ibin3]->Write();

      hJetGenTrackBackground_pythia[ibin][ibin2][ibin3]->Write();
      hJetGenTrackBackground_hydjet[ibin][ibin2][ibin3]->Write();


      } /// ibin3
    } /// ptbin
  }  //centralitybin
  out_file->Close();
}



double mydeltaR(double eta1, double phi1, double eta2, double phi2);
double mydeltaR(double eta1, double phi1, double eta2, double phi2){
  double deltaEta = fabs(eta1-eta2);
  double deltaPhi = fabs(phi1-phi2);
  if (deltaPhi>TMath::TwoPi())   std::cout << "mydeltaR calculation. Why is deltaPhi larger than 2pi?" << std::endl;
  else if (deltaPhi>TMath::Pi()) deltaPhi = TMath::TwoPi()-deltaPhi;
  return pow(deltaEta*deltaEta+deltaPhi*deltaPhi,0.5);
}


void StudyFiles(std::vector<TString> file_names, int foi, hist_class *my_hists);
//void GetFilesOfFileNames(std::vector<TString> &files_of_file_names, std::vector<float> &Xsections); ////getting the data/simulation based on whatever dataset you pick (arg)

void GetFilesOfFileNames(std::vector<TString> &files_of_file_names, std::vector<float> &Xsections, std::vector<double> &assumed_n_evt); ////getting the data/simulation based on whatever dataset you pick (arg);

void ReadFileList(std::vector<TString> &my_file_names, TString file_of_names, bool debug=false);
float GetDatasetWeight(double n_evt_raw, double Xsection);//
//void GetBkgShape(TVector3 highest_jet_vec, TVector3 second_highest_jet_vec, std::vector<TVector3>& bkg_particles);  
//double GetBkgShape(TVector3 highest_jet_vec, TVector3 second_highest_jet_vec, std::vector<TVector3>& bkg_particles);  
void GetBkgShape(TVector3 highest_jet_vec, TVector3 second_highest_jet_vec, TVector3& bkg_dir, TVector3& bkg_dir2);

int main(int argc, char *argv[])
{
  gROOT->ProcessLine("#include <vector>");

    assert(argc == 4);
    dataset_type_code = atoi(argv[1]);    //// pick datasets you want to run over (0 is data Pb, 1 is Hydjet30, etc)
    trkPtCut = atof(argv[2]); /// cut on associated pT
    npart = atoi(argv[3]);  /// part of 10000 events 
    assert(npart >= 0 && npart < e_n_npart_types);
    // assert(npart >= 0);
    assert(trkPtCut > 0. && trkPtCut < 5.);
    std::cout << "Running with trkPtCut " << trkPtCut << " and npart: " << npart << std::endl;
    
//assert(("Length can't possibly be negative! Tell jsmith",  2>1));
  if(dataset_type_code == e_Data2011 ) is_data = true;

  else if( dataset_type_code == e_HydJet30 || dataset_type_code == e_HydJet50 || dataset_type_code == e_HydJet80|| dataset_type_code == e_HydJet100|| dataset_type_code == e_HydJet120|| dataset_type_code == e_HydJet170|| dataset_type_code == e_HydJet200 || dataset_type_code == e_HydJet250 || dataset_type_code == e_HydJet300) is_data =false;
  else assert(0);

  std::vector<TString> files_of_file_names;   files_of_file_names.clear();
  std::vector<float> Xsections;   Xsections.clear();
  std::vector<double> assumed_n_evt;   assumed_n_evt.clear();

  ////// Get list of files for each dataset ... where do I look for the dataset
  std::cout << "GetFilesOf" << std::endl;
  //  GetFilesOfFileNames(files_of_file_names, Xsections);
  GetFilesOfFileNames(files_of_file_names, Xsections, assumed_n_evt);

  //// Total histograms, will add to this for each dataset
  hist_class *hists = new hist_class((TString) ("hists"), is_data);
  //// Loop over each dataset (dataset = collection of similar files of event types)
  for(int foi = 0; foi < (int) files_of_file_names.size(); foi++) {
    TString dataset_str = "sim_dataset_";   dataset_str += foi;
    hist_class *these_hists = new hist_class((TString) ("hists_" + dataset_str), is_data);
    std::cout << "Got hist class for foi " << foi << std::endl;
    std::vector<TString> file_names;   file_names.clear();

    ////// Collect the file names for this dataset
    ReadFileList( file_names, files_of_file_names.at(foi), true);

    ////// Pass these files to StudyFiles function to fill histograms
    StudyFiles(file_names, foi, these_hists);


    float dataset_scale_val = 1.;
    if( !is_data ) {    ////// determine the proper weighting for this dataset (how much you've simulated and likelihood of event)
      assert(foi < (int) Xsections.size() );
      ////////////////////      dataset_scale_val = GetDatasetWeight(these_hists->n_evt_raw, Xsections.at(foi));

      //dataset_scale_val = GetDatasetWeight(assumed_n_evt.at(foi), Xsections.at(foi));


      std::cout << "foi: " << foi << ", xsec: " << Xsections.at(foi) << ", n_evt_raw: " << these_hists->n_evt_raw << ", scale: " << dataset_scale_val << "\n";
    }
    /////// Add the histograms for this dataset to the full histogram class that you care about
    std::cout << "Addhists" << std::endl;
    hists->AddHists(these_hists, dataset_scale_val);
    std::cout << "Delete" << std::endl;
    these_hists->Delete();
    std::cout << "Delete" << std::endl;
    delete these_hists;
  }
  hists->NormalizeByCounts();
  hists->Write();
}



void StudyFiles(std::vector<TString> file_names, int foi, hist_class *my_hists)
{


  double weight=1.;



  //////////###### PTHAT SAMPLES ###########///////////////
  TFile * wtfile_vtx_mc;


  int pthat =0;
  int pthatmax =0;

  if( dataset_type_code == e_HydJet50 ) {
    pthat = 50;
    pthatmax=80;
  }
  if( dataset_type_code == e_HydJet80 ) {
    pthat = 80;
    pthatmax=120;

  }
  if( dataset_type_code == e_HydJet100 ) {
    pthat = 100;
    pthatmax=120;
  }
  if( dataset_type_code == e_HydJet120 ) {
    pthat = 120;
    pthatmax=300;
  }
  if( dataset_type_code == e_HydJet170 ) {
    pthat = 170;
    pthatmax=200;
  }
  if( dataset_type_code == e_HydJet200 ) {
    pthat = 200;
    pthatmax=250;
  }
  if( dataset_type_code == e_HydJet250 ) {
    pthat = 250;
    pthatmax=300;
  }
  if( dataset_type_code == e_HydJet300 ) {
    pthat = 300;
    pthatmax=999;
  }


  ///========= vertex weight =======////

  TH1F* hWeight_vtx;
  TH1F* hWeight_MC_vtx;

  TH1F* hWeight_cent;
  TH1F* hWeight_MC_cent;

  TH1F * ptwt[nCBins];
  TH1F * ptwt_mc[nCBins];


  ////===========  data vertex =========//////
  if( !is_data ) {
    TFile * wtfile_vtx = TFile::Open("pbpb_data_vtx_cent.root", "readonly");
    hWeight_vtx = (TH1F*) ((TH1F*)wtfile_vtx->Get("hists_track_vz"))->Clone("hists_track_vz_clone");
    hWeight_cent = (TH1F*) ((TH1F*)wtfile_vtx->Get("hists_Centrality"))->Clone("hists_centrality_clone");

    hWeight_vtx->Scale(1./hWeight_vtx->Integral());
    hWeight_cent->Scale(1./hWeight_cent->Integral());

    // wtfile_vtx->Close();

    /// pthat80-100-120
    //    TFile * wtfile_vtx_mc = TFile::Open("Hydjet_vtx_cent_reweight.root", "readonly");

    //// pthat80 abd pthat120
    TFile * wtfile_vtx_mc = TFile::Open("Hydjet_vtx_cent_reweight_new.root", "readonly");
    hWeight_MC_vtx = (TH1F*) ((TH1F*)wtfile_vtx_mc->Get("Vertex_histo"))->Clone("hists_track_vz_clone");
    hWeight_MC_cent = (TH1F*) ((TH1F*)wtfile_vtx_mc->Get("Centrality_histo"))->Clone("hists_centrality_clone");


    hWeight_MC_vtx->Scale(1./hWeight_MC_vtx->Integral());
    hWeight_MC_cent->Scale(1./hWeight_MC_cent->Integral());

    TFile * ptwtfile = TFile::Open("pythiahydjet_mc_data_rewfiles/JetTrigMix_inclusive_finaldata_finaltrkcorr.root", "readonly");
    TFile * ptwtfile_mc = TFile::Open("pythiahydjet_mc_data_rewfiles/PYTHIAHYDJET_Merged_RecJetPtSpectrum_Inclusive.root", "readonly");

    for(int ibin = 0 ; ibin <nCBins; ibin++){
      ptwt[ibin] = (TH1F*) ptwtfile->Get("hists_all_jets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_Pt100_Pt300");
      ptwt_mc[ibin] = (TH1F*) ptwtfile_mc->Get("Hydjet_Merged_JetPt_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_Pt100_Pt300");

//      std::cout << " myname " << ptwt[ibin]->GetName() << std::endl;
 ///     std::cout << " myname_mc " << ptwt_mc[ibin]->GetName() << std::endl;

      
      ptwt[ibin]->Scale(1./ptwt[ibin]->Integral());
      ptwt_mc[ibin]->Scale(1./ptwt_mc[ibin]->Integral());
      ptwt[ibin]->Divide(ptwt_mc[ibin]);
    }



  }



  std::cout << "I am working\n";



  ////==  NEW TRACKING EFF CORRECTIONS ==///

    int npt=29;
      double ptmin[]={0.5,0.5,0.5,0.5,0.5,0.55,0.55,0.55,0.55,0.55,0.65,0.65,0.65,0.65,0.65,0.8,0.8,0.8,0.8,0.8, 1, 1, 1, 1, 1, 3, 3, 3, 8};
      double ptmax[]={ 0.55,0.55,0.55,0.55,0.55,0.65,0.65,0.65,0.65,0.65,0.8,0.8,0.8,0.8,0.8, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 8, 8, 8,300};

      int cent_min[]={ 0, 20, 40, 60,100,0, 20, 40, 60,100,0, 20, 40, 60,100,0, 20, 40, 60,100, 0,20,40, 60,100, 0,20, 40, 0};
      int cent_max[]={ 20, 40, 60,100,200,20, 40, 60,100,200,20, 40, 60,100,200,20, 40, 60,100,200,20,40,60,100,200,20,40,200,200};

      //getting histograms for track efficiency correction 
      TFile *f_eff[npt];
      TProfile *p_eff_cent[npt];
      TProfile2D *p_eff_accept[npt];
      TProfile *p_eff_pt[npt];
      TProfile *p_eff_rmin[npt];
      TFile *f_fake[npt];
      TProfile *p_fake_cent[npt];
      TProfile2D *p_fake_accept[npt];
      TProfile *p_fake_pt[npt];
      TProfile *p_fake_rmin[npt];
//cout << "MZ is working" << endl;
      for(int ipt=0; ipt<npt;ipt++){
//MZ//Fixed file path//        f_eff[ipt]= new TFile(Form("../../TrackCorrectionTables_FixedByMIT/akVs3Calo_20140920/eff/eff_pt%d_%d_cent%d_%d.root",(int)(ptmin[ipt]*100),(int)(ptmax[ipt]*100),(int)(0.5*cent_min[ipt]),(int)(0.5*cent_max[ipt])));
        f_eff[ipt]= new TFile(Form("../TrackCorrectionTables/akVs3Calo_20140920/eff/eff_pt%d_%d_cent%d_%d.root",(int)(ptmin[ipt]*100),(int)(ptmax[ipt]*100),(int)(0.5*cent_min[ipt]),(int)(0.5*cent_max[ipt])));

        p_eff_cent[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_cent");
        p_eff_pt[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_pt");
        p_eff_accept[ipt]=(TProfile2D*)f_eff[ipt]->Get("p_eff_acceptance");
        p_eff_rmin[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_rmin");

//MZ//Fixed file path below//        f_fake[ipt]= new TFile(Form("../../TrackCorrectionTables_FixedByMIT/akVs3Calo_20140920/fake/fake_pt%d_%d_cent%d_%d.root",(int)(ptmin[ipt]*100),(int)(ptmax[ipt]*100),(int)(0.5*cent_min[ipt]),(int)(0.5*cent_max[ipt])));
        f_fake[ipt]= new TFile(Form("../TrackCorrectionTables/akVs3Calo_20140920/fake/fake_pt%d_%d_cent%d_%d.root",(int)(ptmin[ipt]*100),(int)(ptmax[ipt]*100),(int)(0.5*cent_min[ipt]),(int)(0.5*cent_max[ipt])));

        p_fake_cent[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_cent");
        p_fake_pt[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_pt");
        p_fake_accept[ipt]=(TProfile2D*)f_fake[ipt]->Get("p_fake_acceptance");
        p_fake_rmin[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_rmin");
      }

      /////////////////////

 /// assert(parti <= (int) file_names.size() );
 ///cout << "MZ I am working " << endl;
  for(int fi = 0; fi < (int) file_names.size(); fi++) {
    //if( parti >= 0 && parti != fi ) continue;
    TFile *my_file = TFile::Open(file_names.at(fi));
    std::cout << "file it" << ", file_name: " << file_names.at(fi) << ", number " << fi << " of " << file_names.size() << std::endl;
    if(my_file->IsZombie()) {
      std::cout << "Is zombie" << std::endl;
    }

    TTree *inp_tree;
    inp_tree = (TTree*)  my_file->Get("mixing_tree");
    mixing_tree *my_primary = new mixing_tree(inp_tree);

    std::cout << "Got CT" << std::endl;

    int n_evt = my_primary->fChain->GetEntriesFast();


    TF1 *fcen = new TF1("fcen","[0]*exp([1]+[2]*x+[3]*x*x+[4]*x*x*x)",0,40);
    TF1 * fVz = new TF1("fVx","[0]+[1]*x+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,3)+[4]*TMath::Power(x,4)", -15., 15.);


    ///==========================   Event Loop starts ===================================
    ///==========================   Event Loop starts ===================================

      
 //cout << "MZ I am working " << endl;
      int min_ev_num = 0;      int max_ev_num = 0;
      min_ev_num = npart*10000;
      max_ev_num = (npart+1)*10000;
      
      if( npart == e_AllParts ) {
          min_ev_num = 0;
          max_ev_num = my_primary->fChain->GetEntriesFast();
          std::cout << "Run over all events from " << min_ev_num << " to " << max_ev_num << "\n";
      }
      //   std::cout << "min_ev_num " << min_ev_num << ", max_ev_num: " << max_ev_num << ", out of possible " << n_evt << "\n";
      for(int evi = min_ev_num; evi < max_ev_num; evi++) {
          if( evi > n_evt ) break;


          my_primary->fChain->GetEntry(evi);
          my_hists->NEvents->Fill(my_primary->fChain->GetEntry(evi));
          my_hists->n_evt_raw++;
          
      if (evi%1000==0) std::cout << " I am running on file " << fi << " of " << ((int) file_names.size()) << ", evi: " << evi << " of " << n_evt <<  "  pthat : " << pthat << "  pthatmax : " << pthatmax  << std::endl;


 //cout << "MZ I am working " << endl;
      if( !is_data ) {
        double pthat_new = my_primary->pthat;
        if (pthat_new > pthatmax) continue;
      }




      //! Centrality reweighting function
      fcen->SetParameters(1.98261e-02,5.55963e+00,-1.34951e-01,1.70895e-03,-9.28386e-05);
      //! vertex z reqeighting
      fVz->SetParameters(7.62788e-01,-1.13228e-02,5.85199e-03,-3.04550e-04,4.43440e-05);

      Int_t hiBin = my_primary->hiBin;
      Float_t vz  = my_primary->vz->at(0);

      //      if(fabs(vz) > 15.)continue;
  //    double wvz=1;
   //   double wcen=1;

      //      my_hists->track_vz->Fill(vz);
        int vtx_bin = hWeight_vtx->GetXaxis()->FindBin(vz);
        float num = hWeight_vtx->GetBinContent(vtx_bin);  /// DATA
        float denom = hWeight_MC_vtx->GetBinContent(vtx_bin); /// MC

// cout << "MZ I am working " << endl;
        double wvz = 1.;
        if( denom > 0.0001 ) wvz = num / denom;

        double wcen=1;
        int cent_bin = hWeight_cent->GetXaxis()->FindBin(hiBin);
        float num_cent = hWeight_cent->GetBinContent(cent_bin);  /// DATA
        float denom_cent = hWeight_MC_cent->GetBinContent(cent_bin); /// MC
        if( denom_cent > 0.0001 ) wcen = num_cent / denom_cent;
       


      if(is_data) { 

        if (my_primary->HLT_HIJet80_v1==0) continue;
        my_hists->NEvents_after_trigger->Fill(1.0);

        int event_selection = my_primary->pcollisionEventSelection;
        if(event_selection==0) continue;
        my_hists->NEvents_test->Fill(1.0);

        int noise_event_selection = my_primary->pHBHENoiseFilter; 
        if(noise_event_selection==0) continue;
        my_hists->NEvents_after_noise->Fill(1.0);


        wvz=1;
        wcen=1;
      }else{
        //        wvz=fVz->Eval(vz);
        //      wcen = fcen->Eval(hiBin);
//        wvz=1;
  //      wcen=1;
      }
      my_hists->track_vz->Fill(vz);
      my_hists->track_vz_weighted->Fill(vz, wvz*wcen);

      my_hists->Centrality->Fill(hiBin);
      my_hists->Centrality_weighted->Fill(hiBin, wvz*wcen);

///      std::cout << " wcent : " << wcen << "  wvz : " << wvz << std::endl;

      if(fabs(vz) > 15.)continue;
      
      my_hists->NEvents_after_spike->Fill(1.0);



      ///////////// -------- Dijet starts ------/////////////////
      ///////////// -------- Dijet starts ------/////////////////
      ///////////// -------- Dijet starts ------/////////////////


      TVector3 highest_jet_vec;
      highest_jet_vec.SetPtEtaPhi(0, 0, 0);
      TVector3 second_highest_jet_vec;
      second_highest_jet_vec.SetPtEtaPhi(0, 0, 0);

      const double etacut = 2.0 ;
      const double dphicut = 5.*(TMath::Pi())/6. ;
      const double leadingjetcut = 120. ;
      //const double leadingjetcut = 115. ;
      const double subleadjetcut = 50. ;

      double lead_pt=0. ;
      int lead_index=-1 ;
      double sublead_pt=0. ;
      int second_highest_idx=-1 ;
      int highest_idx=-1 ;

      if(is_dijet) {

        //search for leading jet
        for(int j4i = 0; j4i < my_primary->jtpt->size() ; j4i++) {
          double jet_pt= my_primary->jtpt->at(j4i);
          if(TMath::Abs(my_primary->jteta->at(j4i))>etacut) continue ;
          if(jet_pt<120) continue ;
          if(jet_pt >lead_pt){
            lead_pt=jet_pt;
            highest_idx=j4i;
            highest_jet_vec.SetPtEtaPhi(jet_pt, my_primary->jteta->at(j4i), my_primary->jtphi->at(j4i));

          }
        } //search for leading jet loop


        //search for subleading jet
        for(int ijet = 0 ; ijet <my_primary->jtpt->size(); ijet++){
          if(ijet==highest_idx) continue ;
          if(TMath::Abs(my_primary->jteta->at(ijet))> etacut) continue ;
          if(my_primary->jtpt->at(ijet)<50) continue ;
          if(my_primary->jtpt->at(ijet) > sublead_pt){
            sublead_pt=my_primary->jtpt->at(ijet);
            second_highest_idx=ijet;
            second_highest_jet_vec.SetPtEtaPhi(my_primary->jtpt->at(ijet), my_primary->jteta->at(ijet), my_primary->jtphi->at(ijet));

          }
        }  //end of subleading jet search

        if( highest_idx < 0 ) continue;
        if( second_highest_idx < 0 ) continue;
        if(my_primary->jtpt->at(highest_idx)< 120. ) continue;
        if( my_primary->jtpt->at(second_highest_idx)< 50) continue;


        float dphi = highest_jet_vec.DeltaPhi(second_highest_jet_vec);
        if (dphi>TMath::Pi()) dphi = 2*TMath::Pi() - dphi;


        for (int ibin=0;ibin<nCBins;ibin++){
          if (my_primary->hiBin >=CBins[ibin] && my_primary->hiBin<CBins[ibin+1]){
            my_hists->dPhi_hist[ibin]->Fill(fabs(dphi));

          }
        }

        my_hists->NEvents_before_dphi->Fill(1.0);

        if(fabs(dphi)<=(5.*TMath::Pi())/6.) continue;
        my_hists->NEvents_after_dphi->Fill(1.0);

        if(TMath::Abs(my_primary->jteta->at(highest_idx)) > 1.6 ) continue;
        if(TMath::Abs(my_primary->jteta->at(second_highest_idx)) > 1.6 ) continue;

        my_hists->NEvents_dijets->Fill(1.0);



        for (int ibin=0;ibin<nCBins;ibin++){
          if (my_primary->hiBin >=CBins[ibin] && my_primary->hiBin<CBins[ibin+1]){

            my_hists->dPhi_after_hist[ibin]->Fill(fabs(dphi));

            double Aj = (my_primary->jtpt->at(highest_idx) - my_primary->jtpt->at(second_highest_idx))/(my_primary->jtpt->at(highest_idx) + my_primary->jtpt->at(second_highest_idx));
            my_hists->Aj[ibin]->Fill(Aj); 


            my_hists->LeadingJetPt_hist[ibin]->Fill(my_primary->jtpt->at(highest_idx));
            my_hists->LeadingJetPhi_hist[ibin]->Fill(my_primary->jtphi->at(highest_idx));
            my_hists->LeadingJetEta_hist[ibin]->Fill(my_primary->jteta->at(highest_idx));

            my_hists->SubJetPt_hist[ibin]->Fill(my_primary->jtpt->at(second_highest_idx));
            my_hists->SubJetPhi_hist[ibin]->Fill(my_primary->jtphi->at(second_highest_idx));
            my_hists->SubJetEta_hist[ibin]->Fill(my_primary->jteta->at(second_highest_idx));

          }
        }

      } /// is_dijet ends


      ///////////// -------- Dijet ends ------/////////////////
      ///////////// -------- Dijet ends ------/////////////////
      ///////////// -------- Dijet ends ------/////////////////


      /// SIGNAL+BACKGROUND///

      double thiscent = my_primary->hiBin;
      double thisz = my_primary->vz->at(0);
      double cent = my_primary->hiBin;

      double thisjetpt =0 ;
      double thisjetphi =0 ;
      double thisjeteta =0 ;
      double Aj = 0;
      double pt_reweight = 1.;

      TVector3 jet_vec;
      TVector3 track_vec;
      jet_vec.SetPtEtaPhi(0, 0, 0);
      track_vec.SetPtEtaPhi(0, 0, 0);

      //// now loop over jets in the original events
      for(int j4i = 0; j4i < (int) my_primary->jtpt->size(); j4i++) {



        if(is_dijet){
          if ( (j4i != highest_idx) && (j4i != second_highest_idx) ) continue;  /// select 2 leading jet
          if(fabs(my_primary->jteta->at(j4i)) > 1.6 ) continue;
          if( my_primary->jtpt->at(highest_idx) <120 ) continue;
          if( my_primary->jtpt->at(highest_idx) >300 ) continue;
          if( my_primary->jtpt->at(second_highest_idx) <50 ) continue;
          if( my_primary->jtpt->at(second_highest_idx) >300 ) continue;
        }else{
          if(fabs(my_primary->jteta->at(j4i)) > 1.6 ) continue;
          if(my_primary->trackMax->at(j4i)/my_primary->jtpt->at(j4i) <=0.01) continue ;
          if( my_primary->jtpt->at(j4i) <120 ) continue;
          if( my_primary->jtpt->at(j4i) >300 ) continue;
        }


        int ibin = 0;      int ibin2 = 0;  int ibin3=0;

  /*      for (int centi=0;centi<nCBins;centi++){
          if (my_primary->hiBin >=CBins[centi] && my_primary->hiBin <CBins[centi+1])  ibin = centi;
        }
*/

        if(is_dijet){
          if ( is_subleadingjet &&  j4i != second_highest_idx ) continue;  /// select subleading jets
          if ( is_leadingjet && j4i != highest_idx ) continue;  /// select leading jet
        }

        for(int ibin = 0; ibin < nCBins; ibin++) {     /// centrality bin loop


        if(ptwt[ibin]->GetBinContent(ptwt[ibin]->FindBin(my_primary->jtpt->at(j4i)))){
          pt_reweight=ptwt[ibin]->GetBinContent(ptwt[ibin]->FindBin(my_primary->jtpt->at(j4i)));
        }else {
          pt_reweight = 1. ;
        }

///       std::cout << " jtpt : " << my_primary->jtpt->at(j4i) << "  reweight: " << pt_reweight << std::endl;


        for(int pti = 0; pti < nPtBins; pti++) {
          if (my_primary->jtpt->at(j4i) >=PtBins[pti] && my_primary->jtpt->at(j4i) < PtBins[pti+1])  ibin2 = pti ;

        }




        thisjetpt = my_primary->jtpt->at(j4i);
        thisjetphi = my_primary->jtphi->at(j4i);
        thisjeteta = my_primary->jteta->at(j4i);

        my_hists->all_jets_corrpT[ibin][ibin2]->Fill(my_primary->jtpt->at(j4i), wvz*wcen*pt_reweight); 
        my_hists->all_jets_phi[ibin][ibin2]->Fill(my_primary->jtphi->at(j4i), wvz*wcen); 
        my_hists->all_jets_eta[ibin][ibin2]->Fill(my_primary->jteta->at(j4i), wvz*wcen); 

        my_hists->only_leadingjets_corrpT[ibin][ibin2]->Fill(my_primary->jtpt->at(j4i),wvz*wcen*pt_reweight);
        my_hists->only_leadingjets_phi[ibin][ibin2]->Fill(my_primary->jtphi->at(j4i), wvz*wcen);
        my_hists->only_leadingjets_eta[ibin][ibin2]->Fill(my_primary->jteta->at(j4i), wvz*wcen);


        //// genparticles

        /// closure with gen particles instead of simTracks ///
        int ngen=0;
        if( !is_data ) {

          int lead_sube = 0;
          int highest_idx_sube=-1 ;

          for(int ipart = 0 ; ipart < (int) my_primary->pt->size(); ipart++){ //sim track loop 
            double gen_pt = my_primary->pt->at(ipart);
            double gen_phi = my_primary->phi->at(ipart);
            double gen_eta = my_primary->eta->at(ipart);
            int chg = my_primary->chg->at(ipart);
            int sube = my_primary->sube->at(ipart);
            if(chg==0) continue ;
            if(gen_pt<trkPtCut)continue ;
            if(TMath::Abs(gen_eta)>2.4)continue ;
            if(sube >lead_sube){
              lead_sube=sube;
              highest_idx_sube=sube;
            }
          }

 
          for(int ipart = 0 ; ipart < (int) my_primary->pt->size(); ipart++){ //sim track loop 
            double gen_pt = my_primary->pt->at(ipart);
            double gen_phi = my_primary->phi->at(ipart);
            double gen_eta = my_primary->eta->at(ipart);
            int chg = my_primary->chg->at(ipart);
            int sube = my_primary->sube->at(ipart);
            
            
           /// std::cout << " sube : " << sube << std::endl;
             
            if(chg==0) continue ;
            if(gen_pt<trkPtCut)continue ;
            if(TMath::Abs(gen_eta)>2.4)continue ;

            for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
              if (gen_pt >=TrkPtBins[trkpti] && gen_pt < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
            } /// trkpti loop

            my_hists->genparticles_pT_hist[ibin][ibin2]->Fill(gen_pt, wvz*wcen);
            my_hists->gensube_hist[ibin][ibin2]->Fill(sube, wvz*wcen);

            double deta_gen = my_primary->jteta->at(j4i) - gen_eta;
            double dphi_gen = my_primary->jtphi->at(j4i) - gen_phi;

              my_hists->hJetGenTrackSignalBackground[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen, wvz*wcen*pt_reweight);
              my_hists->hJetGenTrackSignalBackground[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen+ (2.*TMath::Pi()), wvz*wcen*pt_reweight);
              my_hists->hJetGenTrackSignalBackground[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen- (2.*TMath::Pi()), wvz*wcen*pt_reweight);

              
            if(sube==0) {   // selects pythia
              my_hists->hJetGenTrackSignalBackground_pythia[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen, wvz*wcen*pt_reweight);
              my_hists->hJetGenTrackSignalBackground_pythia[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen+ (2.*TMath::Pi()), wvz*wcen*pt_reweight);
              my_hists->hJetGenTrackSignalBackground_pythia[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen- (2.*TMath::Pi()), wvz*wcen*pt_reweight);

              my_hists->track_gen_pT_pythia[ibin][ibin2]->Fill(gen_pt, wvz*wcen);
              my_hists->track_gen_eta_pythia[ibin][ibin2]->Fill(gen_eta, wvz*wcen);
              my_hists->track_gen_phi_pythia[ibin][ibin2]->Fill(gen_phi, wvz*wcen);
      
              my_hists->GenTrkEta[ibin][ibin2][ibin3]->Fill(gen_eta,wvz*wcen);
              my_hists->GenTrkPhi[ibin][ibin2][ibin3]->Fill(gen_phi,wvz*wcen);
              my_hists->GenTrkPt[ibin][ibin2][ibin3]->Fill(gen_pt,wvz*wcen);


            }else{
              my_hists->hJetGenTrackSignalBackground_hydjet[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen, wvz*wcen*pt_reweight);
              my_hists->hJetGenTrackSignalBackground_hydjet[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen+ (2.*TMath::Pi()), wvz*wcen*pt_reweight);
              my_hists->hJetGenTrackSignalBackground_hydjet[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen- (2.*TMath::Pi()), wvz*wcen*pt_reweight);

              my_hists->track_gen_pT_hydjet[ibin][ibin2]->Fill(gen_pt, wvz*wcen);
              my_hists->track_gen_eta_hydjet[ibin][ibin2]->Fill(gen_eta, wvz*wcen);
              my_hists->track_gen_phi_hydjet[ibin][ibin2]->Fill(gen_phi, wvz*wcen);
            }



            if(sube > 0){
              ///std::cout << " sube : " << sube << "    highestSube :" << lead_sube << std::endl;
              if(sube == lead_sube){
                my_hists->hJetGenTrackSignalBackground_hydro[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen, wvz*wcen*pt_reweight);
                my_hists->hJetGenTrackSignalBackground_hydro[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen+ (2.*TMath::Pi()), wvz*wcen*pt_reweight);
                my_hists->hJetGenTrackSignalBackground_hydro[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen- (2.*TMath::Pi()), wvz*wcen*pt_reweight);
                my_hists->track_gen_pT_hydjet2[ibin][ibin2]->Fill(gen_pt, wvz*wcen);
              }else{
                my_hists->hJetGenTrackSignalBackground_minijets[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen, wvz*wcen*pt_reweight);
                my_hists->hJetGenTrackSignalBackground_minijets[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen+ (2.*TMath::Pi()), wvz*wcen*pt_reweight);
                my_hists->hJetGenTrackSignalBackground_minijets[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen- (2.*TMath::Pi()), wvz*wcen*pt_reweight);
                my_hists->track_gen_pT_hydjet3[ibin][ibin2]->Fill(gen_pt, wvz*wcen);
              }
            }

          } // genparticle loop


          for(int ipart = 0 ; ipart < (int) my_primary->pPt->size() ; ipart++){ //sim track loop 
            double gen_pt = my_primary->pPt->at(ipart);
            double gen_phi = my_primary->pPhi->at(ipart);
            double gen_eta = my_primary->pEta->at(ipart);  
            if(gen_pt<trkPtCut)continue ;

            if(TMath::Abs(gen_eta)>2.4)continue ;
            double rr = mydeltaR(my_primary->jteta->at(j4i),my_primary->jtphi->at(j4i), gen_eta,gen_phi);
            for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
              if (gen_pt >=TrkPtBins[trkpti] && gen_pt < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
            } /// trkpti loop

            double deta_gen = my_primary->jteta->at(j4i) - gen_eta;
            double dphi_gen = my_primary->jtphi->at(j4i) - gen_phi;


            my_hists->track_gen_pT_hist[ibin][ibin2]->Fill(gen_pt, wvz*wcen);
            my_hists->track_gen_eta_hist[ibin][ibin2]->Fill(gen_eta, wvz*wcen);
            my_hists->track_gen_phi_hist[ibin][ibin2]->Fill(gen_phi, wvz*wcen);


            if(rr<0.3){
              ngen += (wvz*wcen);
            }

          }/// simtrack loop 
        } /// isdata!

        ///=================================================================================================
        ///------------------------------------- GenChgParticle Ends----------------------------------------
        ///=================================================================================================


        double ntracks=0;
        float SumpTLeadJet = 0.; 
        double ntracks_away=0;
        Float_t trkweight=1.0;


          float cent=my_primary->hiBin;

        //// reconstructed tracks
        for(int tracks =0; tracks < (int) my_primary->trkPt->size(); tracks++){
          if(fabs(my_primary->trkEta->at(tracks))>2.4) continue;
          if (my_primary->highPurity->at(tracks)==1) {
            if(my_primary->trkPt->at(tracks)>trkPtCut) {

              double pt=my_primary->trkPt->at(tracks);
              double phi=my_primary->trkPhi->at(tracks);
              double eta=my_primary->trkEta->at(tracks);

              double rr = mydeltaR(my_primary->jteta->at(j4i),my_primary->jtphi->at(j4i), my_primary->trkEta->at(tracks),my_primary->trkPhi->at(tracks));

   /*           //find rmin;
              float rmin=99;
              for(int ijet = 0; ijet < (int) my_primary->jtpt->size(); ijet++) {
                if(fabs(my_primary->jteta->at(ijet))>2 || my_primary->jtpt->at(ijet)<50) continue;
                float r_reco=sqrt(pow(eta-my_primary->jteta->at(ijet),2)+pow(acos(cos(phi-my_primary->jtphi->at(ijet))),2));
                if(r_reco<rmin)rmin=r_reco;
              }
              //get efficiency correction for the track   
              float fake_pt,fake_cent,fake_accept,fake_rmin;
              fake_pt=fake_cent=fake_accept=fake_rmin=0;
              float eff_pt,eff_cent,eff_accept,eff_rmin;
              eff_pt=eff_cent=eff_accept=eff_rmin=1;

              //given the pt,centrality,eta,phi and rmin of the track find the factorized efficiencies
              for(int ipt=0;ipt<npt;ipt++){
                if(pt>=ptmin[ipt] && pt<ptmax[ipt] && cent>=cent_min[ipt] && cent<cent_max[ipt]){
                  eff_pt=p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->FindBin(pt));
                  eff_cent=p_eff_cent[ipt]->GetBinContent(p_eff_cent[ipt]->FindBin(cent));
                  eff_accept=p_eff_accept[ipt]->GetBinContent(p_eff_accept[ipt]->GetXaxis()->FindBin(phi),p_eff_accept[ipt]->GetYaxis()->FindBin(eta));
                  if(rmin<5)eff_rmin=p_eff_rmin[ipt]->GetBinContent(p_eff_rmin[ipt]->FindBin(rmin));
                }
              }
              for(int ipt=0;ipt<npt;ipt++){
                if(pt>=ptmin[ipt] && pt<ptmax[ipt] && cent>=cent_min[ipt] && cent<cent_max[ipt]){
                  fake_pt=p_fake_pt[ipt]->GetBinContent(p_fake_pt[ipt]->FindBin(pt));
                  fake_cent=p_fake_cent[ipt]->GetBinContent(p_fake_cent[ipt]->FindBin(cent));
                  fake_accept=p_fake_accept[ipt]->GetBinContent(p_fake_accept[ipt]->GetXaxis()->FindBin(phi),p_fake_accept[ipt]->GetYaxis()->FindBin(eta));
                  if(rmin<5) fake_rmin=p_fake_rmin[ipt]->GetBinContent(p_fake_rmin[ipt]->FindBin(rmin));
                }
              }
              float temp_eff=1;
              temp_eff=eff_accept*eff_cent*eff_pt*eff_rmin;
              float temp_fake=0;
              temp_fake=fake_accept+fake_cent+fake_pt+fake_rmin;
              if(temp_eff==0){ //the if statements are temporary next corrections won't have these
                if(pt>100)temp_eff=0.8;
                else temp_eff=1;
              }

              Float_t trkweight = (1-temp_fake)/temp_eff;
 */
///              std::cout << " weigght : " << trkweight << std::endl;
               
              for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
                if (my_primary->trkPt->at(tracks) >=TrkPtBins[trkpti] && my_primary->trkPt->at(tracks) < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
              } /// trkpti loop

              Float_t trkweight = (1-my_primary->fake->at(tracks))/my_primary->eff->at(tracks);

              my_hists->TrkPt[ibin][ibin2][ibin3]->Fill(my_primary->trkPt->at(tracks),wvz*wcen);
              my_hists->TrkEta[ibin][ibin2][ibin3]->Fill(my_primary->trkEta->at(tracks),wvz*wcen);
              my_hists->TrkPhi[ibin][ibin2][ibin3]->Fill(my_primary->trkPhi->at(tracks),wvz*wcen);

              my_hists->TrkPt_weighted[ibin][ibin2][ibin3]->Fill(my_primary->trkPt->at(tracks),trkweight*wvz*wcen);
              my_hists->TrkEta_weighted[ibin][ibin2][ibin3]->Fill(my_primary->trkEta->at(tracks),trkweight*wvz*wcen);
              my_hists->TrkPhi_weighted[ibin][ibin2][ibin3]->Fill(my_primary->trkPhi->at(tracks),trkweight*wvz*wcen);

              my_hists->track_cand_pT_hist[ibin][ibin2]->Fill(my_primary->trkPt->at(tracks), wvz*wcen);
              my_hists->track_cand_eta_hist[ibin][ibin2]->Fill(my_primary->trkEta->at(tracks), wvz*wcen);
              my_hists->track_cand_phi_hist[ibin][ibin2]->Fill(my_primary->trkPhi->at(tracks), wvz*wcen);

              my_hists->track_cand_pT_hist_weighted[ibin][ibin2]->Fill(my_primary->trkPt->at(tracks), trkweight*wvz*wcen);
              my_hists->track_cand_eta_hist_weighted[ibin][ibin2]->Fill(my_primary->trkEta->at(tracks), trkweight*wvz*wcen);
              my_hists->track_cand_phi_hist_weighted[ibin][ibin2]->Fill(my_primary->trkPhi->at(tracks), trkweight*wvz*wcen);

              double deta = my_primary->jteta->at(j4i) - my_primary->trkEta->at(tracks);
              double dphi = my_primary->jtphi->at(j4i) - my_primary->trkPhi->at(tracks);

              my_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight*wvz*wcen*pt_reweight);
              my_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3]->Fill(deta,dphi+ (2.*TMath::Pi()), trkweight*wvz*wcen*pt_reweight);
              my_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3]->Fill(deta,dphi- (2.*TMath::Pi()), trkweight*wvz*wcen*pt_reweight);

              my_hists->hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen*pt_reweight);
              my_hists->hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi+ (2.*TMath::Pi()), wvz*wcen*pt_reweight);
              my_hists->hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi- (2.*TMath::Pi()), wvz*wcen*pt_reweight);


              jet_vec.SetPtEtaPhi(my_primary->jtpt->at(j4i), my_primary->jteta->at(j4i), my_primary->jtphi->at(j4i));
              track_vec.SetPtEtaPhi(my_primary->trkPt->at(tracks), my_primary->trkEta->at(tracks), my_primary->trkPhi->at(tracks));
              double dphi_vector= jet_vec.DeltaPhi(track_vec);
              my_hists->dPhi_vector[ibin]->Fill(dphi_vector);

              int ybin1 = my_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3]->GetYaxis()->FindBin(-TMath::Pi());
              int ybin2 = my_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3]->GetYaxis()->FindBin(TMath::Pi());
              int nbinsx = my_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3]->GetNbinsX();


              float my_int = 0.;
              for(int bxi = 0; bxi <= my_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3]->GetNbinsX()+1; bxi++ ) {
                float x_cent = my_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3]->GetXaxis()->GetBinCenter(bxi);
                for(int byi = 0; byi <= my_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3]->GetNbinsY()+1; byi++ ) {
                  float y_cent = my_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3]->GetYaxis()->GetBinCenter(byi);
                  if( fabs(y_cent) > TMath::Pi() / 2. ) continue;
                  float cont = my_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3]->GetBinContent(bxi, byi);
                  my_int += cont;
                }
              }
              double integral_2D = my_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3]->Integral(0, nbinsx, ybin1, ybin2);
              float n_jets_int = my_hists->all_jets_corrpT[ibin][ibin2]->Integral();
              double integral_2D_corr = integral_2D / n_jets_int;
              if( evi == 97 ) {
                /// std::cout << "dphi first measure: " << dphi << ", second: " << dphi_vector << ", track num: " << tracks  << ", integral_2D_corr: " << integral_2D_corr << ", my_int: " << my_int << "\n";
              }
              if(fabs(dphi_vector) < (TMath::Pi())/2){
                //                 std::cout << " near_side_dphi : " << dphi << std::endl;
                ntracks += (trkweight*wvz*wcen);
                if( ibin == 0 && ibin2 == 0 && ibin3 == 0 && evi == 97 ) {
                  double mean_1D = my_hists->Ntracks[ibin][ibin2][ibin3]->GetMean();
                  ///std::cout << "track num: " << tracks << ", pT: " << my_primary->trkPt->at(tracks) << ", eta: " << my_primary->trkEta->at(tracks) << ", phi: " << my_primary->trkPhi->at(tracks) << ", dphi: " << fabs(dphi_vector) << ", ibin: " << ibin << ", ibin2: " << ibin2 << ", ibin3: " << ibin3 << ", integral_2D: " << integral_2D << ", n_jets_int: " << n_jets_int << ", integral_2D_corr: " << integral_2D_corr << ", my_int: " << my_int << ", ntracks: " << ntracks << ", mean_1D: " << mean_1D << "\n";
                }
              }else{
                ///               std::cout << " away_side_dphi : " << dphi << std::endl;
                ntracks_away += (trkweight*wvz*wcen);
              }





              if (rr<0.3){
                SumpTLeadJet  += my_primary->trkPt->at(tracks)*trkweight;
              }
            }   // track pt cut
            } ///track algo and purity
          } /// track loop


          my_hists->Ntracks[ibin][ibin2][ibin3]->Fill(ntracks);
          my_hists->Ntracks_awayside[ibin][ibin2][ibin3]->Fill(ntracks_away);
          my_hists->NumChargedParticles[ibin][ibin2]->Fill(ntracks);
          my_hists->Ngen[ibin][ibin2][ibin3]->Fill(ngen);

          float my_int = 0.;
          for(int bxi = 0; bxi <= my_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3]->GetNbinsX()+1; bxi++ ) {
            float x_cent = my_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3]->GetXaxis()->GetBinCenter(bxi);
            for(int byi = 0; byi <= my_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3]->GetNbinsY()+1; byi++ ) {
              float y_cent = my_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3]->GetYaxis()->GetBinCenter(byi);
              if( fabs(y_cent) > TMath::Pi() / 2. ) continue;
              float cont = my_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3]->GetBinContent(bxi, byi);
              my_int += cont;
            }
          }
          int ybin1 = my_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3]->GetYaxis()->FindBin(-TMath::Pi());
          int ybin2 = my_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3]->GetYaxis()->FindBin(TMath::Pi());
          int nbinsx = my_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3]->GetNbinsX();
          double integral_2D = my_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3]->Integral(0, nbinsx, ybin1, ybin2);
          float n_jets_int = my_hists->all_jets_corrpT[ibin][ibin2]->Integral();
          double integral_2D_corr = integral_2D / n_jets_int;
          double mean_1D = my_hists->Ntracks[ibin][ibin2][ibin3]->GetMean();
          float my_int_corr = my_int / n_jets_int;
          ///std::cout << "evi: " << evi << ", ibin: " << ibin << ", ibin2: " << ibin2 << ", ibin3: " << ibin3 << ", integral_2D: " << integral_2D << ", n_jets_int: " << n_jets_int << ", integral_2D_corr: " << integral_2D_corr << ", my 2D int: " << my_int_corr << ", mean_1D: " << mean_1D << "\n";

          //// reconstructed tracks with eta reflection
          if( fabs(my_primary->jteta->at(j4i)) > 0.3 ){

            for(int tracks =0; tracks < (int) my_primary->trkPt->size(); tracks++){
              if(fabs(my_primary->trkEta->at(tracks))>2.4) continue;
              if (my_primary->highPurity->at(tracks)==1) {
                if(my_primary->trkPt->at(tracks)>trkPtCut) {

                  double rr = mydeltaR(-(my_primary->jteta->at(j4i)),my_primary->jtphi->at(j4i), my_primary->trkEta->at(tracks),my_primary->trkPhi->at(tracks));

                  Float_t trkweight_er = (1-my_primary->fake->at(tracks))/my_primary->eff->at(tracks);

                  if(rr<0.3){
                    my_hists->EtaRef_bkg_pt_weighted[ibin][ibin2]->Fill(my_primary->trkPt->at(tracks), trkweight_er*wvz*wcen);
                    my_hists->EtaRef_bkg_pt[ibin][ibin2]->Fill(my_primary->trkPt->at(tracks),wvz*wcen);
                  }

                } /// trkpt cut
              }  /// purity and algo
            } /// trk loop
          } // jet eta loop



          ///   Numerator for signal
          const int NRBIN=6;
          std::vector<double> SumPtRBin(NRBIN);
          for(int vi = 0; vi < (int) SumPtRBin.size(); vi++) {
            SumPtRBin.at(vi) = 0.;
          }
          const float RBIN[NRBIN+1]={0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3};
          double ptParticle=0;

          for(int tracks =0; tracks < (int) my_primary->trkPt->size(); tracks++){
            if(fabs(my_primary->trkEta->at(tracks))>2.4) continue;
            if (my_primary->highPurity->at(tracks)==1) {
              if(my_primary->trkPt->at(tracks)>trkPtCut) {

                double rr = mydeltaR(my_primary->jteta->at(j4i),my_primary->jtphi->at(j4i), my_primary->trkEta->at(tracks),my_primary->trkPhi->at(tracks));

                Float_t trkweight_num = (1-my_primary->fake->at(tracks))/my_primary->eff->at(tracks);

                ptParticle = my_primary->trkPt->at(tracks);


                for(int ir=0;ir<NRBIN;ir++){
                  if(RBIN[ir]<rr && rr <= RBIN[ir+1]){
                    SumPtRBin[ir]+=ptParticle*trkweight_num;
                    continue;
                  }
                }   

              } /// trkpt cut
            }  /// purity and algo
          } /// trk loop




          ///=================================================================================================
          ///------------------------ Filling Shapes histograms  ---------------------------------------------
          ///=================================================================================================


          std::vector<double> DifferentialPt(NRBIN);

          for(int ir=0;ir<NRBIN;ir++){
            if(my_primary->jtpt->at(j4i)>0.0001){
              DifferentialPt[ir] =SumPtRBin[ir]/(my_primary->jtpt->at(j4i));
            }
          }

          for(int ir=0;ir<NRBIN;ir++){
            float fill_val = 0.05*(ir+1)-0.025;
            if(my_primary->jtpt->at(j4i)>0.0001){
              my_hists->JetShapeDiffParticles_1D[ibin][ibin2]->Fill(fill_val,DifferentialPt[ir]*wvz*wcen);
            }
          }
        } //// cent loop
        } /// jet loop in the original event



      ///////===== EMix starts ===========/////////
      ///////===== EMix starts ===========/////////
      ///////===== EMix starts ===========/////////
/*

      if( thisjetpt < 1. ) continue;    /// failed event selection

      ///====  Here: Load MB tree ===////
      int n_MBevs_found = 0;
      int mixentry = evi;

      while(1) {
        if( evi % 100 == 0 ) {
          //  std::cout << "Mixing for evi " << evi << ", mixentry " << mixentry << std::endl;
        }
        if( n_MBevs_found > 30 ) break;
        mixentry++;
        if( mixentry >= n_evt ) mixentry=0;
        if( mixentry == evi ) break;
        my_primary->GetEntry(mixentry);
        if( !is_data ) {
          double pthat_mix_new = my_primary->pthat;
          if (pthat_mix_new > pthatmax) continue;
        }

        /// event matching
//        if( fabs(thisz - my_primary->vz->at(0)) > 1.0 ) continue;
//        if( fabs(thisz - my_primary->vz->at(0)) > 5.0 ) continue;
  ///      if( fabs(thiscent-my_primary->hiBin) / thiscent > 0.25 ) continue;

  /*      bool is_centmatched = false;
        for (int cbin=0; cbin<nCentralityBins; cbin++){
          if (my_primary->hiBin >=CentralityBins[cbin] && my_primary->hiBin<CentralityBins[cbin+1]){
            if (thiscent >=CentralityBins[cbin] && thiscent<CentralityBins[cbin+1]){
              is_centmatched = true;
              break;
            }
          }
        }
        if(!is_centmatched) continue;

*/

/*        n_MBevs_found++;

        int ibin = 0;      int ibin2 = 0;  int ibin3=0;
        for (int centi=0;centi<nCBins;centi++){
          if (my_primary->hiBin >=CBins[centi] && my_primary->hiBin <CBins[centi+1])  ibin = centi;
        } 

        for(int pti = 0; pti < nPtBins; pti++) { 
          if (thisjetpt >=PtBins[pti] && thisjetpt < PtBins[pti+1])  ibin2 = pti ;
        } /// pti loop

        my_hists->emix_jetpt[ibin][ibin2]->Fill(1, wvz*wcen);

        float trkweight_bkg = 1.;

        if( !is_data ) {
          /// closure with gen particles instead of simTracks ///

          for(int ipart = 0 ; ipart < (int) my_primary->pt->size(); ipart++){ //sim track loop 
            double gen_pt = my_primary->pt->at(ipart);
            double gen_phi = my_primary->phi->at(ipart);
            double gen_eta = my_primary->eta->at(ipart);
            int chg = my_primary->chg->at(ipart);
            int sube = my_primary->sube->at(ipart);

            if(chg==0) continue ;
            if(gen_pt<trkPtCut)continue ;

            if(TMath::Abs(gen_eta)>2.4)continue ;
            my_hists->genparticles_bkg_pT_hist[ibin][ibin2]->Fill(gen_pt, wvz*wcen);

            double deta_gen = thisjeteta - gen_eta;
            double dphi_gen = thisjetphi - gen_phi;

            if(sube==0) {   // selects pythia
              my_hists->hJetGenTrackBackground_pythia[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen, wvz*wcen);
              my_hists->hJetGenTrackBackground_pythia[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen+ (2.*TMath::Pi()), wvz*wcen);
              my_hists->hJetGenTrackBackground_pythia[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen- (2.*TMath::Pi()), wvz*wcen);

            }else{    // selects hydjet part
              my_hists->hJetGenTrackBackground_hydjet[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen, wvz*wcen);
              my_hists->hJetGenTrackBackground_hydjet[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen+ (2.*TMath::Pi()), wvz*wcen);
              my_hists->hJetGenTrackBackground_hydjet[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen- (2.*TMath::Pi()), wvz*wcen);

            }
          } // genparticle loop


         for(int ipart = 0 ; ipart < (int) my_primary->pPt->size(); ipart++){ //sim track loop 
            double gen_pt = my_primary->pPt->at(ipart);
            double gen_phi = my_primary->pPhi->at(ipart);
            double gen_eta = my_primary->pEta->at(ipart);
            if(gen_pt<trkPtCut)continue ;

            if(TMath::Abs(gen_eta)>2.4)continue ;
            double rr = mydeltaR(thisjeteta, thisjetphi, gen_eta,gen_phi);

            for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
              if (gen_pt >=TrkPtBins[trkpti] && gen_pt < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
            } /// trkpti loop

            double deta_bkg_gen = thisjeteta - gen_eta;
            double dphi_bkg_gen = thisjetphi - gen_phi;

            my_hists->GenTrkEta_bkg[ibin][ibin2][ibin3]->Fill(gen_eta,wvz*wcen);
            my_hists->GenTrkPhi_bkg[ibin][ibin2][ibin3]->Fill(gen_phi,wvz*wcen);

            my_hists->hJetGenTrackBackground[ibin][ibin2][ibin3]->Fill(deta_bkg_gen,dphi_bkg_gen,wvz*wcen);
            my_hists->hJetGenTrackBackground[ibin][ibin2][ibin3]->Fill(deta_bkg_gen,dphi_bkg_gen+(2*TMath::Pi()),wvz*wcen);
            my_hists->hJetGenTrackBackground[ibin][ibin2][ibin3]->Fill(deta_bkg_gen,dphi_bkg_gen-(2*TMath::Pi()),wvz*wcen);

            my_hists->track_gen_bkg_pT_hist[ibin][ibin2]->Fill(gen_pt, wvz*wcen);
            my_hists->track_gen_bkg_eta_hist[ibin][ibin2]->Fill(gen_eta, wvz*wcen);
            my_hists->rr_gen_leadingjet[ibin][ibin3]->Fill(rr, wvz*wcen);


          } // simtrk loop
        } /// isdata!




        ////==== place reco level ====///////
        for(int tracks =0; tracks < (int) my_primary->trkPt->size(); tracks++){
          if (my_primary->highPurity->at(tracks)==1) {
            if( fabs(my_primary->trkEta->at(tracks)) > 2.4 ) continue;
            if(my_primary->trkPt->at(tracks)>trkPtCut) {

              double rr = mydeltaR(thisjeteta, thisjetphi, my_primary->trkEta->at(tracks), my_primary->trkPhi->at(tracks));

              double pt=my_primary->trkPt->at(tracks);
              double phi=my_primary->trkPhi->at(tracks);
              double eta=my_primary->trkEta->at(tracks);

              //find rmin;
              float rmin=99;
              for(int ijet = 0; ijet < (int) my_primary->jtpt->size(); ijet++) {
                if(fabs(thisjeteta)>2 || thisjetpt<50) continue;
                float r_reco=sqrt(pow((eta-thisjeteta),2)+pow(acos(cos(phi-thisjetphi)),2));
                if(r_reco<rmin)rmin=r_reco;
              }
              //get efficiency correction for the track   
              float fake_pt,fake_cent,fake_accept,fake_rmin;
              fake_pt=fake_cent=fake_accept=fake_rmin=0;
              float eff_pt,eff_cent,eff_accept,eff_rmin;
              eff_pt=eff_cent=eff_accept=eff_rmin=1;

              //given the pt,centrality,eta,phi and rmin of the track find the factorized efficiencies
              for(int ipt=0;ipt<npt;ipt++){
                if(pt>=ptmin[ipt] && pt<ptmax[ipt] && cent>=cent_min[ipt] && cent<cent_max[ipt]){
                  eff_pt=p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->FindBin(pt));
                  eff_cent=p_eff_cent[ipt]->GetBinContent(p_eff_cent[ipt]->FindBin(cent));
                  eff_accept=p_eff_accept[ipt]->GetBinContent(p_eff_accept[ipt]->GetXaxis()->FindBin(phi),p_eff_accept[ipt]->GetYaxis()->FindBin(eta));
                  if(rmin<5)eff_rmin=p_eff_rmin[ipt]->GetBinContent(p_eff_rmin[ipt]->FindBin(rmin));
                }
              }
              for(int ipt=0;ipt<npt;ipt++){
                if(pt>=ptmin[ipt] && pt<ptmax[ipt] && cent>=cent_min[ipt] && cent<cent_max[ipt]){
                  fake_pt=p_fake_pt[ipt]->GetBinContent(p_fake_pt[ipt]->FindBin(pt));
                  fake_cent=p_fake_cent[ipt]->GetBinContent(p_fake_cent[ipt]->FindBin(cent));
                  fake_accept=p_fake_accept[ipt]->GetBinContent(p_fake_accept[ipt]->GetXaxis()->FindBin(phi),p_fake_accept[ipt]->GetYaxis()->FindBin(eta));
                  if(rmin<5) fake_rmin=p_fake_rmin[ipt]->GetBinContent(p_fake_rmin[ipt]->FindBin(rmin));
                }
              }
              float temp_eff=1;
              temp_eff=eff_accept*eff_cent*eff_pt*eff_rmin;
              float temp_fake=0;
              temp_fake=fake_accept+fake_cent+fake_pt+fake_rmin;
              if(temp_eff==0){ //the if statements are temporary next corrections won't have these
                if(pt>100)temp_eff=0.8;
                else temp_eff=1;
              }

              Float_t trkweight_bkg = (1-temp_fake)/temp_eff;





              for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) { 
                if (my_primary->trkPt->at(tracks) >=TrkPtBins[trkpti] && my_primary->trkPt->at(tracks) < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
              } /// trkpti loop

///               Float_t trkweight_bkg = (1-my_primary->fake->at(tracks))/my_primary->eff->at(tracks);

              double deta_bkg = thisjeteta - my_primary->trkEta->at(tracks);
              double dphi_bkg = thisjetphi - my_primary->trkPhi->at(tracks);

              my_hists->TrkEta_bkg[ibin][ibin2][ibin3]->Fill(my_primary->trkEta->at(tracks),trkweight_bkg*wvz*wcen);
              my_hists->TrkPhi_bkg[ibin][ibin2][ibin3]->Fill(my_primary->trkPhi->at(tracks),trkweight_bkg*wvz*wcen);
              my_hists->hJetTrackBackground[ibin][ibin2][ibin3]->Fill(deta_bkg,dphi_bkg,trkweight_bkg*wvz*wcen);
              my_hists->hJetTrackBackground[ibin][ibin2][ibin3]->Fill(deta_bkg,dphi_bkg+(2*TMath::Pi()),trkweight_bkg*wvz*wcen);
              my_hists->hJetTrackBackground[ibin][ibin2][ibin3]->Fill(deta_bkg,dphi_bkg-(2*TMath::Pi()),trkweight_bkg*wvz*wcen);

              my_hists->hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta_bkg,dphi_bkg,wvz*wcen);
              my_hists->hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta_bkg,dphi_bkg+(2*TMath::Pi()),wvz*wcen);
              my_hists->hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta_bkg,dphi_bkg-(2*TMath::Pi()),wvz*wcen);

              my_hists->track_bkg_pT_hist[ibin][ibin2]->Fill(my_primary->trkPt->at(tracks), wvz*wcen);
              my_hists->track_bkg_pT_hist_weighted[ibin][ibin2]->Fill(my_primary->trkPt->at(tracks), trkweight_bkg*wvz*wcen);

              my_hists->track_bkg_eta_hist[ibin][ibin2]->Fill(my_primary->trkEta->at(tracks), wvz*wcen);
              my_hists->track_bkg_eta_hist_weighted[ibin][ibin2]->Fill(my_primary->trkEta->at(tracks), trkweight_bkg*wvz*wcen);

              my_hists->rr_reco_leadingjet[ibin][ibin3]->Fill(rr, wvz*wcen);
              my_hists->rr_reco_leadingjet_weighted[ibin][ibin3]->Fill(rr, trkweight_bkg*wvz*wcen);


              if(rr<0.3)
              {
              }

            } /// trkpt cut
          } /// high purity and trkAlgo cut	    
        } /// trk loop


      }  /// looping over Mixed events until you find 20


      ////==== end reco level =====/////////
      my_primary->GetEntry(evi);

      my_hists->NumberOfMatches->Fill(n_MBevs_found);

*/

 
      if( evi % 500 == 0 ) my_hists->Write();

    }  ///event loop

      ///==========================   Event Loop ends ==================================================
      ///==========================   Event Loop ends ==================================================
      ///==========================   Event Loop ends ==================================================


   std::cout << "done with file " << fi << std::endl;

  }  /// file loop

}




    ///========================= background subtruction ============================================================
    ///========================= background subtruction ============================================================
    ///========================= background subtruction ============================================================


    void GetBkgShape(TVector3 highest_jet_vec, TVector3 second_highest_jet_vec, TVector3& bkg_dir, TVector3& bkg_dir2) { 


      //TVector3 bkg_dir = highest_jet_vec.Cross(second_highest_jet_vec);

      bkg_dir = highest_jet_vec;
      bkg_dir2 = second_highest_jet_vec;

      float bkg_dir_Eta = -bkg_dir.Eta();
      float bkg_dir_Phi = bkg_dir.Phi();
      float bkg_dir_Pt = bkg_dir.Pt();

    float bkg_dir_Eta2 = -bkg_dir2.Eta();
    float bkg_dir_Phi2 = bkg_dir2.Phi();
    float bkg_dir_Pt2 = bkg_dir2.Pt();

    bkg_dir.SetPtEtaPhi(bkg_dir_Pt, bkg_dir_Eta, bkg_dir_Phi);
    bkg_dir2.SetPtEtaPhi(bkg_dir_Pt2, bkg_dir_Eta2, bkg_dir_Phi2);

  }

  void ReadFileList(std::vector<TString> &my_file_names, TString file_of_names, bool debug)
  {
    ifstream file_stream(file_of_names);
    std::string line;
    my_file_names.clear();
    if( debug ) std::cout << "Open file " << file_of_names << " to extract files to run over" << std::endl;
    if( file_stream.is_open() ) {
      if( debug ) std::cout << "Opened " << file_of_names << " for reading" << std::endl;
      int line_num = 0;
      while( !file_stream.eof() ) {
	getline(file_stream, line);
	if( debug ) std::cout << line_num << ": " << line << std::endl;
	TString tstring_line(line);
	if( tstring_line.CompareTo("", TString::kExact) != 0 ) my_file_names.push_back(tstring_line);
	line_num++;
      }
    } else {
      std::cout << "Error, could not open " << file_of_names << " for reading" << std::endl;
      assert(0);
    }
  }




  //void GetFilesOfFileNames(std::vector<TString> &files_of_file_names, std::vector<float> &Xsections) ////getting the data/simulation based on whatever dataset you pick (arg)
  void GetFilesOfFileNames(std::vector<TString> &files_of_file_names, std::vector<float> &Xsections, std::vector<double> &assumed_n_evt) ////getting the data/simulation based on whatever dataset you pick (arg)


  {
    files_of_file_names.clear();   Xsections.clear();

    if( dataset_type_code == e_HydJet30 ) {
      files_of_file_names.push_back("Hydjet_Pt30.txt");
      Xsections.push_back(1.079e-02);
      assumed_n_evt.push_back(109918);

    }else if( dataset_type_code == e_HydJet50 ) {
      files_of_file_names.push_back("Hydjet_Pt50.txt");
      Xsections.push_back(1.021e-03 - 9.913E-05);
      assumed_n_evt.push_back(112166);

    }else if( dataset_type_code == e_HydJet80 ) {
      files_of_file_names.push_back("Hydjet_Pt80.txt");
      Xsections.push_back(9.913E-05 - 3.069E-05);
      assumed_n_evt.push_back(146162);

    }else if( dataset_type_code == e_HydJet100 ) {
      files_of_file_names.push_back("Hydjet_Pt100.txt");
      Xsections.push_back(3.069E-05 - 1.128E-05);
      assumed_n_evt.push_back(593463);

    }else if( dataset_type_code == e_HydJet120 ) {
      files_of_file_names.push_back("Hydjet_Pt120.txt");
      Xsections.push_back(1.128E-05 - 1.470E-06);
      assumed_n_evt.push_back(141263);

    }else if( dataset_type_code == e_HydJet170 ) {
      files_of_file_names.push_back("Hydjet_Pt170.txt");
      Xsections.push_back(1.470E-06 - 5.310E-07);
      assumed_n_evt.push_back(62944);

    }else if( dataset_type_code == e_HydJet200 ) {
      files_of_file_names.push_back("Hydjet_Pt200.txt");
      Xsections.push_back(5.310E-07 - 1.192E-07);
      assumed_n_evt.push_back(54952);

    }else if( dataset_type_code == e_HydJet250 ) {
      files_of_file_names.push_back("Hydjet_Pt250.txt");
      Xsections.push_back(1.192E-07 - 3.176E-08);
      assumed_n_evt.push_back(37856);

    }else if( dataset_type_code == e_HydJet300 ) {
      files_of_file_names.push_back("Hydjet_Pt300.txt");
      Xsections.push_back(3.176E-08);
      assumed_n_evt.push_back(53009);

    }else if( dataset_type_code == e_Data2011 ) {
      files_of_file_names.push_back("Data2011.txt"); 
    }else {
      std::cout << "I don't understand the fileset" << std::endl;
      assert(0);
    }
  }



  float GetDatasetWeight(double n_evt_raw, double Xsection)
  {
    assert( n_evt_raw > 0.5 );

    double wt  = Xsection / n_evt_raw;

    return wt;
  }





