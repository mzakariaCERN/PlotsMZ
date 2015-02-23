#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2D.h"

#include "class_def/JetAna.h"
#include "class_def/Tracks.h"
#include "class_def/HLT.h"
#include "class_def/HiTree.h"
#include "class_def/Skim.h"
#include "class_def/GenParticles.h"
#include "class_def/pfcand.h"  // needed to access particle flow candidates, make sure the two files are executable and included in the compuling script

#include "TH2F.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TCut.h>

#include "assert.h"
#include <fstream>
#include "TMath.h"
#include "TVector3.h"
#include <vector>

using namespace std;


class hist_class {
  public:
    hist_class(TString the_desc);
    void Delete();
    void Write();
    TString desc;

    //  TH1F* NEvents;




};


hist_class::hist_class(TString the_desc)
{

  desc = the_desc;
  //NEvents = new TH1F((TString) (desc + "_Nevents"), "", 100, 0., 100.);
}

void hist_class::Delete()
{
  //delete NEvents;

}

void hist_class::Write()
{

  // TString out_name = (TString) ("root_output/" + desc + "_Data_JetID_test.root");
  TString out_name = (TString) ("root_output/" + desc + "_MC_Pthat80.root");

  TFile *out_file = new TFile(out_name, "RECREATE");

  //NEvents->Write();

  out_file->Close();
}




int main(int argc, char *argv[])
{
  assert( argc == 2);
  int output_file_num = atoi(argv[1]);


  //////////###### centrality Bins ###########///////////////

  std::cout << "I am working\n";
  TTree *inp_tree;
  TTree *inp_tree3;
  TTree *inp_tree4;
  TTree *inp_tree5;
  TTree *inp_tree6;
  TTree *inp_tree7;
  TTree *inp_tree2;

  //// PYTHIA+HYDJET files:
  ///TFile *my_file = TFile::Open("/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet120_FOREST_Track8_Jet24_FixedPtHatJES_v0/0.root");
  //TFile *my_file = TFile::Open("/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet100_FOREST_Track8_Jet24_FixedPtHatJES_v0/0.root");
 // TFile *my_file = TFile::Open("/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet80_FOREST_Track8_Jet24_FixedPtHatJES_v0/0.root");
//  TFile *my_file = TFile::Open("/afs/cern.ch/user/m/mozakari/CMSSW_5_3_20/src/FromPelin/ntuple_maker_pbpb/06B1A8BB-EDCA-E311-87A8-00266CF9C018.root");
//         TFile *my_file = TFile::Open("/data/mzakaria/PythiaHydjet_OfficialProduction_20150211/HiForest_PYTHIA_HYDJET_pthat370_Track9_Jet30_matchEqR_merged_forest_0.root");

  // TFile *my_file = TFile::Open("/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet50_FOREST_Track8_Jet24_FixedPtHatJES_v0/0.root");
  ///TFile *my_file = TFile::Open("/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet30_FOREST_Track8_Jet24_FixedPtHatJES_v0/0.root");
  TFile *my_file = TFile::Open("/data/mzakaria/PythiaHydjet_OfficialProduction_20150211/HiForest_PYTHIA_HYDJET_pthat80_Track9_Jet30_matchEqR_merged_forest_0.root");

  /// PbPb data:
  //TFile *my_file = TFile::Open("/mnt/hadoop/cms/store/user/dgulhan/PbPbForest_MatchEqR_Calo_HIHighPt_HIRun2011-14Mar2014-v4_merged_nodublicate/PbPbForest_MatchEqR_Calo_HIHighPt_HIRun2011-14Mar2014-v4.root");

  /// MINBIAS DATA
  ///   TFile *my_file = TFile::Open("/mnt/hadoop/cms/store/user/velicanu/HIMinBias2011_GR_R_53_LV6_CMSSW_5_3_16_Forest_Track8_Jet21/0.root");

  /// MINBIAS DATA
    //TFile *my_file = TFile::Open("/mnt/hadoop/cms/store/user/velicanu/HI_53X_MB_Hydjet_STARTHI53_LV1_5_3_16_trk8_jet26_2Apr2014_1200EST_FOREST/0.root");



  if(my_file->IsZombie()) { 
    std::cout << "Is zombie" << std::endl;
  }    

  inp_tree = (TTree*)  my_file->Get("akVs3CaloJetAnalyzer/t");
  inp_tree3 = (TTree*) my_file->Get("hiEvtAnalyzer/HiTree");
  inp_tree4 = (TTree*) my_file->Get("skimanalysis/HltTree");
  inp_tree5 = (TTree*) my_file->Get("anaTrack/trackTree");
  inp_tree6 = (TTree*) my_file->Get("hltanalysis/HltTree");
  inp_tree7 = (TTree*) my_file->Get("HiGenParticleAna/hi");
  inp_tree2 = (TTree*) my_file->Get("pfcandAnalyzer/pfTree");

  JetAna *my_ct = new JetAna(inp_tree);
  HiTree *my_ct3 = new HiTree(inp_tree3);
  Skim *my_ct4   = new Skim(inp_tree4);
  Tracks *my_ct5 = new Tracks(inp_tree5);
  HLT *my_ct6 = new HLT(inp_tree6);
  GenParticles *my_ct7 = new GenParticles(inp_tree7);
  pfcand *my_ct2 = new pfcand(inp_tree2);


  std::cout << "Got CT" << std::endl;

  int n_evt = my_ct->fChain->GetEntriesFast();



  ///    TString output_file_base = "root_output/PyhtiaHydjet30_mini_ntuple_v1";
  // TString output_file_base = "root_output/PyhtiaHydjet50_mini_ntuple_v1";
  TString output_file_base = "root_output/PyhtiaHydjet80_mini_ntuple_v1";
  ///TString output_file_base = "root_output/PyhtiaHydjet100_mini_ntuple_v1";
  //TString output_file_base = "root_output/PyhtiaHydjet120_mini_ntuple_v1";
  //    TString output_file_base = "root_output/PyhtiaHydjet120_mini_ntuple_v3";

  TString output_file_extension = "_p";   output_file_extension += output_file_num;   output_file_extension += ".root";
  TFile *output_file = new TFile((TString) (output_file_base + output_file_extension), "RECREATE");
  output_file->cd();
  TTree *mixing_tree = new TTree("mixing_tree", "");

  Int_t mult = -999;
  Int_t nTrk = -999;
  Int_t nParticle = -999;
  Int_t nPFpart = -999;  // this variable used to loop over PF candidates

  std::vector<Int_t> *sube= new std::vector<Int_t>();  sube->clear();
  std::vector<Float_t> *pt = new std::vector<Float_t>();  pt->clear();
  std::vector<Float_t> *phi = new std::vector<Float_t>();  phi->clear();
  std::vector<Float_t> *eta = new std::vector<Float_t>();  eta->clear();
  std::vector<Int_t> *chg= new std::vector<Int_t>();  chg->clear();
  std::vector<Float_t> *pPt = new std::vector<Float_t>();  pPt->clear();
  std::vector<Float_t> *mpPt = new std::vector<Float_t>();  mpPt->clear();
  std::vector<Float_t> *pPhi = new std::vector<Float_t>();  pPhi->clear();
  std::vector<Float_t> *pEta = new std::vector<Float_t>();  pEta->clear();
  std::vector<Float_t> *pNRec = new std::vector<Float_t>();  pNRec->clear();
  std::vector<Float_t> *mtrkAlgo = new std::vector<Float_t>();  mtrkAlgo->clear();
  std::vector<Int_t> *mtrkQual = new std::vector<Int_t>();  mtrkQual->clear();
  std::vector<Float_t> *trkEta = new std::vector<Float_t>();   trkEta->clear();
  std::vector<Float_t> *trkPhi = new std::vector<Float_t>();   trkPhi->clear();
  std::vector<Float_t> *trkPt = new std::vector<Float_t>();   trkPt->clear();
  std::vector<Float_t> *trkAlgo = new std::vector<Float_t>();   trkAlgo->clear();
  std::vector<Float_t> *highPurity = new std::vector<Float_t>();   highPurity->clear();
  std::vector<Float_t> *vz = new std::vector<Float_t>();   vz->clear();
  //std::vector<Float_t> *trkweight = new std::vector<Float_t>();   trkweight->clear();
  std::vector<Float_t> *jteta = new std::vector<Float_t>();   jteta->clear();
  std::vector<Float_t> *jtphi = new std::vector<Float_t>();   jtphi->clear();
  std::vector<Float_t> *jtpt = new std::vector<Float_t>();   jtpt->clear();
  std::vector<Float_t> *rawpt = new std::vector<Float_t>();   rawpt->clear(); //From HT, rawpt is jtpt with jff-dep corrections
  std::vector<Float_t> *refeta = new std::vector<Float_t>();  refeta->clear();
  std::vector<Float_t> *refphi = new std::vector<Float_t>();  refphi->clear();
  std::vector<Float_t> *refpt = new std::vector<Float_t>();   refpt->clear();


  std::vector<Float_t> *refdrjt = new std::vector<Float_t>();   refdrjt->clear();
  std::vector<Float_t> *refparton_pt = new std::vector<Float_t>();   refparton_pt->clear();
  std::vector<Int_t> *refparton_flavor = new std::vector<Int_t>();  refparton_flavor->clear();


  std::vector<Float_t> *trackMax = new std::vector<Float_t>();   trackMax->clear();
  std::vector<Float_t> *trkDxy1 = new std::vector<Float_t>();   trkDxy1->clear();
  std::vector<Float_t> *trkDxyError1 = new std::vector<Float_t>();   trkDxyError1->clear();
  std::vector<Float_t> *trkDz1 = new std::vector<Float_t>();   trkDz1->clear();
  std::vector<Float_t> *trkDzError1 = new std::vector<Float_t>();   trkDzError1->clear();
  std::vector<Float_t> *trkPtError = new std::vector<Float_t>();   trkPtError->clear();
 
  std::vector<Float_t> *eff = new std::vector<Float_t>();   eff->clear();
  std::vector<Float_t> *fake = new std::vector<Float_t>();   fake->clear();
  std::vector<Float_t> *trkfake = new std::vector<Float_t>();   trkfake->clear();


  std::vector<Float_t> *mtrkDxy1 = new std::vector<Float_t>();   mtrkDxy1->clear();
  std::vector<Float_t> *mtrkDxyError1 = new std::vector<Float_t>();   mtrkDxyError1->clear();
  std::vector<Float_t> *mtrkDz1 = new std::vector<Float_t>();   mtrkDz1->clear();
  std::vector<Float_t> *mtrkDzError1 = new std::vector<Float_t>();   mtrkDzError1->clear();
  std::vector<Float_t> *mtrkPtError = new std::vector<Float_t>();   mtrkPtError->clear();
  // std::vector<Float_t> *meff = new std::vector<Float_t>();   meff->clear();

  std::vector<Float_t> *geneta = new std::vector<Float_t>();  geneta->clear();
  std::vector<Float_t> *genphi = new std::vector<Float_t>();  genphi->clear();
  std::vector<Float_t> *genpt = new std::vector<Float_t>();   genpt->clear();

  std::vector<Int_t> *pfId = new std::vector<Int_t>(); pfId->clear();
  std::vector<Float_t> *pfPt = new std::vector<Float_t>(); pfPt->clear();
  std::vector<Float_t> *pfVsPt = new std::vector<Float_t>(); pfVsPt->clear();
  std::vector<Float_t> *pfEta = new std::vector<Float_t>(); pfEta->clear();
  std::vector<Float_t> *pfPhi = new std::vector<Float_t>(); pfPhi->clear();
  std::vector<Float_t> *sumpt = new std::vector<Float_t>(); sumpt->clear();

  Int_t pHBHENoiseFilter = -999;
  Int_t pcollisionEventSelection = -999;
  Int_t HLT_HIJet80_v1 = -999;
  Int_t HLT_HIMinBiasHfOrBSC_v1 = -999;
  Int_t hiBin = -999;
  Float_t pthat = -999;

  /// higenparticles
  mixing_tree->Branch("mult", &mult, "mult/I");
  mixing_tree->Branch("pt", "vector<Float_t>", &pt);
  mixing_tree->Branch("phi", "vector<Float_t>", &phi);
  mixing_tree->Branch("eta", "vector<Float_t>", &eta);
  mixing_tree->Branch("chg", "vector<Int_t>", &chg);
  mixing_tree->Branch("sube", "vector<Int_t>", &sube);

  mixing_tree->Branch("nTrk", &nTrk, "nTrk/I");
  mixing_tree->Branch("nParticle", &nParticle, "nParticle/I");
  mixing_tree->Branch("pPt", "vector<Float_t>", &pPt);
  mixing_tree->Branch("mpPt", "vector<Float_t>", &mpPt);
  mixing_tree->Branch("pPhi", "vector<Float_t>", &pPhi);
  mixing_tree->Branch("pEta", "vector<Float_t>", &pEta);
  mixing_tree->Branch("pNRec", "vector<Float_t>", &pNRec);
  mixing_tree->Branch("mtrkAlgo", "vector<Float_t>", &mtrkAlgo);
  mixing_tree->Branch("mtrkQual", "vector<Int_t>", &mtrkQual);
  mixing_tree->Branch("trkEta", "vector<Float_t>", &trkEta);
  mixing_tree->Branch("trkPhi", "vector<Float_t>", &trkPhi);
  mixing_tree->Branch("trkPt", "vector<Float_t>", &trkPt);
  mixing_tree->Branch("trkAlgo", "vector<Float_t>", &trkAlgo);
  mixing_tree->Branch("highPurity", "vector<Float_t>", &highPurity);
  mixing_tree->Branch("vz", "vector<Float_t>", &vz);
  // mixing_tree->Branch("trkweight", "vector<Float_t>", &trkweight);
  mixing_tree->Branch("pHBHENoiseFilter", &pHBHENoiseFilter, "pHBHENoiseFilter/I");
  mixing_tree->Branch("pcollisionEventSelection", &pcollisionEventSelection, "pcollisionEventSelection/I");
  mixing_tree->Branch("HLT_HIJet80_v1", &HLT_HIJet80_v1, "HLT_HIJet80_v1/I");
  mixing_tree->Branch("HLT_HIMinBiasHfOrBSC_v1", &HLT_HIMinBiasHfOrBSC_v1, "HLT_HIMinBiasHfOrBSC_v1/I");

  mixing_tree->Branch("hiBin", &hiBin, "hiBin/I");
  mixing_tree->Branch("jteta", "vector<Float_t>", &jteta);
  mixing_tree->Branch("jtphi", "vector<Float_t>", &jtphi);
  mixing_tree->Branch("jtpt", "vector<Float_t>", &jtpt);
  mixing_tree->Branch("rawpt", "vector<Float_t>", &rawpt);

  mixing_tree->Branch("pthat", &pthat, "pthat/F");
  mixing_tree->Branch("refeta", "vector<Float_t>", &refeta);
  mixing_tree->Branch("refphi", "vector<Float_t>", &refphi);
  mixing_tree->Branch("refpt", "vector<Float_t>", &refpt);
  mixing_tree->Branch("trackMax", "vector<Float_t>", &trackMax);
  mixing_tree->Branch("refdrjt", "vector<Float_t>", &refdrjt);
  mixing_tree->Branch("refparton_pt", "vector<Float_t>", &refparton_pt);
  mixing_tree->Branch("refparton_flavor", "vector<Int_t>", &refparton_flavor);
  mixing_tree->Branch("trkDxy1", "vector<Float_t>", &trkDxy1);
  mixing_tree->Branch("trkDxyError1", "vector<Float_t>", &trkDxyError1);
  mixing_tree->Branch("trkDz1", "vector<Float_t>", &trkDz1);
  mixing_tree->Branch("trkDzError1", "vector<Float_t>", &trkDzError1);
  mixing_tree->Branch("trkPtError", "vector<Float_t>", &trkPtError);
  mixing_tree->Branch("eff", "vector<Float_t>", &eff);
  mixing_tree->Branch("fake", "vector<Float_t>", &fake);
  mixing_tree->Branch("trkfake", "vector<Float_t>", &trkfake);
  mixing_tree->Branch("mtrkDxy1", "vector<Float_t>", &mtrkDxy1);
  mixing_tree->Branch("mtrkDxyError1", "vector<Float_t>", &mtrkDxyError1);
  mixing_tree->Branch("mtrkDz1", "vector<Float_t>", &mtrkDz1);
  mixing_tree->Branch("mtrkDzError1", "vector<Float_t>", &mtrkDzError1);
  mixing_tree->Branch("mtrkPtError", "vector<Float_t>", &mtrkPtError);
  //mixing_tree->Branch("meff", "vector<Float_t>", &meff);
  mixing_tree->Branch("geneta", "vector<Float_t>", &geneta);
  mixing_tree->Branch("genphi", "vector<Float_t>", &genphi);
  mixing_tree->Branch("genpt", "vector<Float_t>", &genpt);

  mixing_tree->Branch("nPFpart", &nPFpart, "nPFpart/I");
  mixing_tree->Branch("pfId", "vector<Int_t>", &pfId);
  mixing_tree->Branch("pfPt", "vector<Float_t>", &pfPt);
  mixing_tree->Branch("pfVsPt", "vector<Float_t>", &pfVsPt);
  mixing_tree->Branch("pfEta", "vector<Float_t>", &pfEta);
  mixing_tree->Branch("pfPhi", "vector<Float_t>", &pfPhi);
  mixing_tree->Branch("sumpt", "vector<Float_t>", &sumpt);
  


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

  for(int ipt=0; ipt<npt;ipt++){
    f_eff[ipt]= new TFile(Form("TrackCorrectionTables_FixedByMIT/akVs3Calo_20140920/eff/eff_pt%d_%d_cent%d_%d.root",(int)(ptmin[ipt]*100),(int)(ptmax[ipt]*100),(int)(0.5*cent_min[ipt]),(int)(0.5*cent_max[ipt])));

    p_eff_cent[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_cent");
    p_eff_pt[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_pt");
    p_eff_accept[ipt]=(TProfile2D*)f_eff[ipt]->Get("p_eff_acceptance");
    p_eff_rmin[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_rmin");

    f_fake[ipt]= new TFile(Form("TrackCorrectionTables_FixedByMIT/akVs3Calo_20140920/fake/fake_pt%d_%d_cent%d_%d.root",(int)(ptmin[ipt]*100),(int)(ptmax[ipt]*100),(int)(0.5*cent_min[ipt]),(int)(0.5*cent_max[ipt])));

    p_fake_cent[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_cent");
    p_fake_pt[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_pt");
    p_fake_accept[ipt]=(TProfile2D*)f_fake[ipt]->Get("p_fake_acceptance");
    p_fake_rmin[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_rmin");
  }



  output_file->cd();

  int ev_min = output_file_num*10000;
  int ev_max = ev_min + 9999;

  cout << "ev_min: " << ev_min << ", Entries: " << n_evt << std::endl;

  assert( ev_min < n_evt );

  if( ev_max >= n_evt ) ev_max = n_evt;
  std::cout << "Will run from event number " << ev_min << " to " << ev_max << "\n";
  for(int evi = ev_min; evi < ev_max; evi++) {

    if( evi % 1000 == 0 )  std::cout << "evi: " << evi <<  " of " << n_evt << "\n";


    my_ct->fChain->GetEntry(evi);

    my_ct2->fChain->GetEntry(evi);  //used to include PF
    my_ct3->fChain->GetEntry(evi);
    my_ct4->fChain->GetEntry(evi);
    my_ct5->fChain->GetEntry(evi);
    my_ct6->fChain->GetEntry(evi);
    my_ct7->fChain->GetEntry(evi);

    //std::cout << "evi: " << evi << ", number of tracks: " << my_ct5->nTrk << std::endl;
    nTrk = my_ct5->nTrk;
    mult = my_ct7->mult;
    nParticle = my_ct5->nParticle; 
    nPFpart = my_ct2->nPFpart;


    if( evi % 1000 == 0 ) std::cout << "Filled successfully" << std::endl;



    //// genparticles
    for(int part = 0; part < my_ct7->mult; part++) {
      if( fabs(my_ct7->eta[part])>2.4) continue;
      if( my_ct7->pt[part] < 1. ) continue;
      int charge = my_ct7->chg[part];
      if(charge==0) continue ;

      pt->push_back(my_ct7->pt[part]);
      phi->push_back(my_ct7->phi[part]);
      eta->push_back(my_ct7->eta[part]);
      chg->push_back(my_ct7->chg[part]);
      sube->push_back(my_ct7->sube[part]);

    }


    for(int j4i = 0; j4i < my_ct->ngen ; j4i++) {
      geneta->push_back(my_ct->geneta[j4i]);
      genphi->push_back(my_ct->genphi[j4i]); 
      genpt->push_back(my_ct->genpt[j4i]);
    }   


    for(int pfi = 0; pfi< my_ct2->nPFpart ; pfi++) {


      pfId->push_back(my_ct2->pfId[pfi]);
      pfPt->push_back(my_ct2->pfPt[pfi]);
      pfVsPt->push_back(my_ct2->pfVsPt[pfi]);
      pfEta->push_back(my_ct2->pfEta[pfi]);
      pfPhi->push_back(my_ct2->pfPhi[pfi]);
      sumpt->push_back(my_ct2->sumpt[pfi]);

    } /// particle flow candidate loop


    for(int j4i = 0; j4i < my_ct->nref ; j4i++) {
      refeta->push_back(my_ct->refeta[j4i]);
      refphi->push_back(my_ct->refphi[j4i]);
      refpt->push_back(my_ct->refpt[j4i]);

      refdrjt->push_back(my_ct->refdrjt[j4i]);
      refparton_pt->push_back(my_ct->refparton_pt[j4i]);
      refparton_flavor->push_back(my_ct->refparton_flavor[j4i]);
    }



    for(int j4i = 0; j4i < my_ct->nref ; j4i++) {

      jteta->push_back(my_ct->jteta[j4i]);
      jtphi->push_back(my_ct->jtphi[j4i]);
      jtpt->push_back(my_ct->jtpt[j4i]);
      rawpt->push_back(my_ct->rawpt[j4i]);

      trackMax->push_back(my_ct->trackMax[j4i]);

    } /// jet loop


    ////====== New Reco Tracks =====///
    float cent=my_ct3->hiBin;

    //loop over tracks
    for(int itrk=0;itrk<my_ct5->nParticle;itrk++){

      float eta=my_ct5->pEta[itrk];
      if(fabs(eta)>2.4) continue; //acceptance of the tracker   

      float pt=my_ct5->pPt[itrk];
      //if(pt<0.5) continue; //acceptance of the tracker
      if(pt<1.0) continue; //acceptance of the tracker


      pPt->push_back(my_ct5->pPt[itrk]);
      mpPt->push_back(my_ct5->mtrkPt[itrk]);
      pPhi->push_back(my_ct5->pPhi[itrk]);
      pEta->push_back(my_ct5->pEta[itrk]);

      mtrkDxy1 -> push_back(my_ct5->mtrkDxy1[itrk]);
      mtrkDxyError1 -> push_back(my_ct5->mtrkDxyError1[itrk]);
      mtrkDz1 -> push_back(my_ct5->mtrkDz1[itrk]);
      mtrkDzError1 -> push_back(my_ct5->mtrkDzError1[itrk]);
      mtrkPtError -> push_back(my_ct5->mtrkPtError[itrk]);
      pNRec->push_back(my_ct5->pNRec[itrk]);
      mtrkAlgo->push_back(my_ct5->mtrkAlgo[itrk]);
      mtrkQual->push_back(my_ct5->mtrkQual[itrk]);

    } /// simtrk loop

    //// reco track loop
    for(int itrk=0;itrk<my_ct5->nTrk;itrk++){

      if(my_ct5->highPurity[itrk]==0 || my_ct5->trkPtError[itrk]/my_ct5->trkPt[itrk]>=0.1 || TMath::Abs(my_ct5->trkDz1[itrk]/my_ct5->trkDzError1[itrk])>=3.0 ||TMath::Abs(my_ct5->trkDxy1[itrk]/my_ct5->trkDxyError1[itrk])>=3.0) continue ;


      float eta=my_ct5->trkEta[itrk];
      if(fabs(eta)>2.4) continue; //acceptance of the tracker   

      float pt=my_ct5->trkPt[itrk];
      //if(pt<0.5) continue; //acceptance of the tracker
      if(pt<1.0) continue; //acceptance of the tracker

      //acceptance of the tracker
      float phi=my_ct5->trkPhi[itrk];

      //float trkfake=my_ct5->trkFake[itrk];
      float rmin=99;

      //find rmin; 

      for(int ijet=0;ijet<my_ct->nref;ijet++){
	if(fabs(my_ct->jteta[ijet])>2 || my_ct->jtpt[ijet]<50) continue;
	float r_reco=sqrt(pow(eta-my_ct->jteta[ijet],2)+pow(acos(cos(phi-my_ct->jtphi[ijet])),2));
	if(r_reco<rmin)rmin=r_reco;
      }

      //get efficiency correction for the track   
      float eff_accept=1;
      float eff_pt=1;
      float eff_cent=1;
      float eff_rmin=1;

      float fake_pt,fake_cent,fake_accept,fake_rmin;
      fake_pt=fake_cent=fake_accept=fake_rmin=0;


      //given the pt,centrality,eta,phi and rmin of the track find the factorized efficiencies

      for(int ipt=0;ipt<npt;ipt++){
	//if(pt>=ptmin[ipt] && pt<ptmax[ipt] && (0.5*cent)>=cent_min[ipt] && (0.5*cent)<cent_max[ipt]){
	if(pt>=ptmin[ipt] && pt<ptmax[ipt] && cent>=cent_min[ipt] && cent<cent_max[ipt]){
	  eff_pt=p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->FindBin(pt));
	  eff_cent=p_eff_cent[ipt]->GetBinContent(p_eff_cent[ipt]->FindBin(cent));
	  eff_accept=p_eff_accept[ipt]->GetBinContent(p_eff_accept[ipt]->GetXaxis()->FindBin(phi),p_eff_accept[ipt]->GetYaxis()->FindBin(eta));
	  if(rmin<5)eff_rmin=p_eff_rmin[ipt]->GetBinContent(p_eff_rmin[ipt]->FindBin(rmin));

	  fake_pt=p_fake_pt[ipt]->GetBinContent(p_fake_pt[ipt]->FindBin(pt));
	  fake_cent=p_fake_cent[ipt]->GetBinContent(p_fake_cent[ipt]->FindBin(cent));
	  fake_accept=p_fake_accept[ipt]->GetBinContent(p_fake_accept[ipt]->GetXaxis()->FindBin(phi),p_fake_accept[ipt]->GetYaxis()->FindBin(eta));
	  if(rmin<5)fake_rmin=p_fake_rmin[ipt]->GetBinContent(p_fake_rmin[ipt]->FindBin(rmin));
	}
      }
      float temp_eff=1;
      temp_eff=eff_accept*eff_cent*eff_pt*eff_rmin;
      float temp_fake=0;
      if(pt<100)temp_fake=fake_accept+fake_cent+fake_pt+fake_rmin;
      if(temp_eff==0){
	cout<<"zero efficiency"<<" eta="<<eta<<" pt="<<pt<<" phi="<<phi<<" cent="<<cent<<endl;
	if(pt>100)temp_eff=0.8;
	else temp_eff=1;
      }

      trkEta->push_back(my_ct5->trkEta[itrk]);
      trkPhi->push_back(my_ct5->trkPhi[itrk]);
      trkPt->push_back(my_ct5->trkPt[itrk]);
      trkAlgo->push_back(my_ct5->trkAlgo[itrk]);
      highPurity->push_back(my_ct5->highPurity[itrk]);

      trkDxy1 -> push_back(my_ct5->trkDxy1[itrk]);
      trkDxyError1 -> push_back(my_ct5->trkDxyError1[itrk]);
      trkDz1 -> push_back(my_ct5->trkDz1[itrk]);
      trkDzError1 -> push_back(my_ct5->trkDzError1[itrk]);
      trkPtError -> push_back(my_ct5->trkPtError[itrk]);


      eff->push_back(temp_eff);
      fake->push_back(temp_fake);
      trkfake -> push_back(my_ct5->trkFake[itrk]);


      }




      //// ======  End of Tracking eff corrections ====== ////

      pHBHENoiseFilter = my_ct4->pHBHENoiseFilter;
      pcollisionEventSelection = my_ct4->pcollisionEventSelection;
      HLT_HIJet80_v1 = my_ct6->HLT_HIJet80_v7;
      HLT_HIMinBiasHfOrBSC_v1 = my_ct6->HLT_HIMinBiasHfOrBSC_v1;
      hiBin = my_ct3->hiBin;
      vz->push_back(my_ct3->vz);
      pthat = my_ct->pthat;

      ///// Fill it
      mixing_tree->Fill();


      ///// Reset
      pt->clear();
      phi->clear();
      eta->clear();
      chg->clear();
      sube->clear();

      pPt->clear();
      mpPt->clear();
      pPhi->clear();
      pEta->clear();
      pNRec->clear();
      mtrkAlgo->clear();
      mtrkQual->clear();

      trkEta->clear();
      trkPhi->clear();
      trkPt->clear();
      trkAlgo->clear();
      highPurity->clear();

      vz->clear();

      //   trkweight->clear();
      jteta->clear();
      jtphi->clear();
      jtpt->clear();
      rawpt->clear();

      refeta->clear();
      refphi->clear();
      refpt->clear();
      refdrjt->clear();
      refparton_pt->clear();
      refparton_flavor->clear();
      trackMax->clear();
      trkDxy1->clear();
      trkDxyError1->clear();
      trkDz1->clear();
      trkDzError1->clear();
      trkPtError->clear();
      eff->clear();
      fake->clear();
      trkfake->clear();
      mtrkDxy1->clear();
      mtrkDxyError1->clear();
      mtrkDz1->clear();
      mtrkDzError1->clear();
      mtrkPtError->clear();
      //meff->clear();
      geneta->clear();
      genphi->clear();
      genpt->clear();

	
 	pfId->clear();
    pfPt->clear();
    pfVsPt->clear();
    pfEta->clear();
    pfPhi->clear();
    sumpt->clear();


    }  ///event loop


    mixing_tree->Write();
    output_file->Write();
    output_file->Close();

  }




