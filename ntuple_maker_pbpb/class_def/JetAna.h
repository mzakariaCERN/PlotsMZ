//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Mar 13 20:15:30 2014 by ROOT version 5.32/00
// from TTree t/akVs3CalopatJets Jet Analysis Tree
// found on file: /tmp/HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0.root
//////////////////////////////////////////////////////////

#ifndef JetAna_h
#define JetAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class JetAna {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           evt;
   Float_t         b;
   Int_t           nref;
   Float_t         rawpt[134];   //[nref]
   Float_t         jtpt[134];   //[nref]
   Float_t         jteta[134];   //[nref]
   Float_t         jty[134];   //[nref]
   Float_t         jtphi[134];   //[nref]
   Float_t         jtpu[134];   //[nref]
   Float_t         jtm[134];   //[nref]
   Float_t         discr_fr01[134];   //[nref]
   Float_t         trackMax[134];   //[nref]
   Float_t         trackSum[134];   //[nref]
   Int_t           trackN[134];   //[nref]
   Float_t         trackHardSum[134];   //[nref]
   Int_t           trackHardN[134];   //[nref]
   Float_t         chargedMax[134];   //[nref]
   Float_t         chargedSum[134];   //[nref]
   Int_t           chargedN[134];   //[nref]
   Float_t         chargedHardSum[134];   //[nref]
   Int_t           chargedHardN[134];   //[nref]
   Float_t         photonMax[134];   //[nref]
   Float_t         photonSum[134];   //[nref]
   Int_t           photonN[134];   //[nref]
   Float_t         photonHardSum[134];   //[nref]
   Int_t           photonHardN[134];   //[nref]
   Float_t         neutralMax[134];   //[nref]
   Float_t         neutralSum[134];   //[nref]
   Int_t           neutralN[134];   //[nref]
   Float_t         hcalSum[134];   //[nref]
   Float_t         ecalSum[134];   //[nref]
   Float_t         eMax[134];   //[nref]
   Float_t         eSum[134];   //[nref]
   Int_t           eN[134];   //[nref]
   Float_t         muMax[134];   //[nref]
   Float_t         muSum[134];   //[nref]
   Int_t           muN[134];   //[nref]
   Int_t           beamId1;
   Int_t           beamId2;
   Float_t         pthat;
   Float_t         refpt[134];   //[nref]
   Float_t         refeta[134];   //[nref]
   Float_t         refy[134];   //[nref]
   Float_t         refphi[134];   //[nref]
   Float_t         refdphijt[134];   //[nref]
   Float_t         refdrjt[134];   //[nref]
   Float_t         refparton_pt[134];   //[nref]
   Int_t           refparton_flavor[134];   //[nref]
   Int_t           refparton_flavorForB[134];   //[nref]
   Float_t         genChargedSum[134];   //[nref]
   Float_t         genHardSum[134];   //[nref]
   Float_t         signalChargedSum[134];   //[nref]
   Float_t         signalHardSum[134];   //[nref]
   Int_t           subid[134];   //[nref]
   Int_t           ngen;
   Int_t           genmatchindex[10];   //[ngen]
   Float_t         genpt[10];   //[ngen]
   Float_t         geneta[10];   //[ngen]
   Float_t         geny[10];   //[ngen]
   Float_t         genphi[10];   //[ngen]
   Float_t         gendphijt[10];   //[ngen]
   Float_t         gendrjt[10];   //[ngen]
   Int_t           gensubid[10];   //[ngen]

   // List of branches
   TBranch        *b_evt;   //!
   TBranch        *b_b;   //!
   TBranch        *b_nref;   //!
   TBranch        *b_rawpt;   //!
   TBranch        *b_jtpt;   //!
   TBranch        *b_jteta;   //!
   TBranch        *b_jty;   //!
   TBranch        *b_jtphi;   //!
   TBranch        *b_jtpu;   //!
   TBranch        *b_jtm;   //!
   TBranch        *b_discr_fr01;   //!
   TBranch        *b_trackMax;   //!
   TBranch        *b_trackSum;   //!
   TBranch        *b_trackN;   //!
   TBranch        *b_trackHardSum;   //!
   TBranch        *b_trackHardN;   //!
   TBranch        *b_chargedMax;   //!
   TBranch        *b_chargedSum;   //!
   TBranch        *b_chargedN;   //!
   TBranch        *b_chargedHardSum;   //!
   TBranch        *b_chargedHardN;   //!
   TBranch        *b_photonMax;   //!
   TBranch        *b_photonSum;   //!
   TBranch        *b_photonN;   //!
   TBranch        *b_photonHardSum;   //!
   TBranch        *b_photonHardN;   //!
   TBranch        *b_neutralMax;   //!
   TBranch        *b_neutralSum;   //!
   TBranch        *b_neutralN;   //!
   TBranch        *b_hcalSum;   //!
   TBranch        *b_ecalSum;   //!
   TBranch        *b_eMax;   //!
   TBranch        *b_eSum;   //!
   TBranch        *b_eN;   //!
   TBranch        *b_muMax;   //!
   TBranch        *b_muSum;   //!
   TBranch        *b_muN;   //!
   TBranch        *b_beamId1;   //!
   TBranch        *b_beamId2;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_refpt;   //!
   TBranch        *b_refeta;   //!
   TBranch        *b_refy;   //!
   TBranch        *b_refphi;   //!
   TBranch        *b_refdphijt;   //!
   TBranch        *b_refdrjt;   //!
   TBranch        *b_refparton_pt;   //!
   TBranch        *b_refparton_flavor;   //!
   TBranch        *b_refparton_flavorForB;   //!
   TBranch        *b_genChargedSum;   //!
   TBranch        *b_genHardSum;   //!
   TBranch        *b_signalChargedSum;   //!
   TBranch        *b_signalHardSum;   //!
   TBranch        *b_subid;   //!
   TBranch        *b_ngen;   //!
   TBranch        *b_genmatchindex;   //!
   TBranch        *b_genpt;   //!
   TBranch        *b_geneta;   //!
   TBranch        *b_geny;   //!
   TBranch        *b_genphi;   //!
   TBranch        *b_gendphijt;   //!
   TBranch        *b_gendrjt;   //!
   TBranch        *b_gensubid;   //!

   JetAna(TTree *tree=0);
   virtual ~JetAna();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef JetAna_cxx
JetAna::JetAna(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
  /*    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/tmp/HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/tmp/HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/tmp/HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0.root:/akVs3CaloJetAnalyzer");
      dir->GetObject("t",tree);
*/
   }
   Init(tree);
}

JetAna::~JetAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t JetAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t JetAna::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void JetAna::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("b", &b, &b_b);
   fChain->SetBranchAddress("nref", &nref, &b_nref);
   fChain->SetBranchAddress("rawpt", rawpt, &b_rawpt);
   fChain->SetBranchAddress("jtpt", jtpt, &b_jtpt);
   fChain->SetBranchAddress("jteta", jteta, &b_jteta);
   fChain->SetBranchAddress("jty", jty, &b_jty);
   fChain->SetBranchAddress("jtphi", jtphi, &b_jtphi);
   fChain->SetBranchAddress("jtpu", jtpu, &b_jtpu);
   fChain->SetBranchAddress("jtm", jtm, &b_jtm);
   fChain->SetBranchAddress("discr_fr01", discr_fr01, &b_discr_fr01);
   fChain->SetBranchAddress("trackMax", trackMax, &b_trackMax);
   fChain->SetBranchAddress("trackSum", trackSum, &b_trackSum);
   fChain->SetBranchAddress("trackN", trackN, &b_trackN);
   fChain->SetBranchAddress("trackHardSum", trackHardSum, &b_trackHardSum);
   fChain->SetBranchAddress("trackHardN", trackHardN, &b_trackHardN);
   fChain->SetBranchAddress("chargedMax", chargedMax, &b_chargedMax);
   fChain->SetBranchAddress("chargedSum", chargedSum, &b_chargedSum);
   fChain->SetBranchAddress("chargedN", chargedN, &b_chargedN);
   fChain->SetBranchAddress("chargedHardSum", chargedHardSum, &b_chargedHardSum);
   fChain->SetBranchAddress("chargedHardN", chargedHardN, &b_chargedHardN);
   fChain->SetBranchAddress("photonMax", photonMax, &b_photonMax);
   fChain->SetBranchAddress("photonSum", photonSum, &b_photonSum);
   fChain->SetBranchAddress("photonN", photonN, &b_photonN);
   fChain->SetBranchAddress("photonHardSum", photonHardSum, &b_photonHardSum);
   fChain->SetBranchAddress("photonHardN", photonHardN, &b_photonHardN);
   fChain->SetBranchAddress("neutralMax", neutralMax, &b_neutralMax);
   fChain->SetBranchAddress("neutralSum", neutralSum, &b_neutralSum);
   fChain->SetBranchAddress("neutralN", neutralN, &b_neutralN);
   fChain->SetBranchAddress("hcalSum", hcalSum, &b_hcalSum);
   fChain->SetBranchAddress("ecalSum", ecalSum, &b_ecalSum);
   fChain->SetBranchAddress("eMax", eMax, &b_eMax);
   fChain->SetBranchAddress("eSum", eSum, &b_eSum);
   fChain->SetBranchAddress("eN", eN, &b_eN);
   fChain->SetBranchAddress("muMax", muMax, &b_muMax);
   fChain->SetBranchAddress("muSum", muSum, &b_muSum);
   fChain->SetBranchAddress("muN", muN, &b_muN);
   fChain->SetBranchAddress("beamId1", &beamId1, &b_beamId1);
   fChain->SetBranchAddress("beamId2", &beamId2, &b_beamId2);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("refpt", refpt, &b_refpt);
   fChain->SetBranchAddress("refeta", refeta, &b_refeta);
   fChain->SetBranchAddress("refy", refy, &b_refy);
   fChain->SetBranchAddress("refphi", refphi, &b_refphi);
   fChain->SetBranchAddress("refdphijt", refdphijt, &b_refdphijt);
   fChain->SetBranchAddress("refdrjt", refdrjt, &b_refdrjt);
   fChain->SetBranchAddress("refparton_pt", refparton_pt, &b_refparton_pt);
   fChain->SetBranchAddress("refparton_flavor", refparton_flavor, &b_refparton_flavor);
   fChain->SetBranchAddress("refparton_flavorForB", refparton_flavorForB, &b_refparton_flavorForB);
   fChain->SetBranchAddress("genChargedSum", genChargedSum, &b_genChargedSum);
   fChain->SetBranchAddress("genHardSum", genHardSum, &b_genHardSum);
   fChain->SetBranchAddress("signalChargedSum", signalChargedSum, &b_signalChargedSum);
   fChain->SetBranchAddress("signalHardSum", signalHardSum, &b_signalHardSum);
   fChain->SetBranchAddress("subid", subid, &b_subid);
   fChain->SetBranchAddress("ngen", &ngen, &b_ngen);
   fChain->SetBranchAddress("genmatchindex", genmatchindex, &b_genmatchindex);
   fChain->SetBranchAddress("genpt", genpt, &b_genpt);
   fChain->SetBranchAddress("geneta", geneta, &b_geneta);
   fChain->SetBranchAddress("geny", geny, &b_geny);
   fChain->SetBranchAddress("genphi", genphi, &b_genphi);
   fChain->SetBranchAddress("gendphijt", gendphijt, &b_gendphijt);
   fChain->SetBranchAddress("gendrjt", gendrjt, &b_gendrjt);
   fChain->SetBranchAddress("gensubid", gensubid, &b_gensubid);
   Notify();
}

Bool_t JetAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void JetAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t JetAna::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef JetAna_cxx
