//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar  4 16:33:17 2014 by ROOT version 5.14/00e
// from TTree HltTree/
// found on file: /mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet80/Track8_Jet19_STARTHI53_LV1/merged_300k_2/HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0.root
//////////////////////////////////////////////////////////

#ifndef Skim_h
#define Skim_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class Skim {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Int_t           photonStep;
   Int_t           extrapatstep;
   Int_t           ana_step;
   Int_t           phltJetHI;
   Int_t           pcollisionEventSelection;
   Int_t           pHBHENoiseFilter;
   Int_t           phiEcalRecHitSpikeFilter;
   Int_t           hltAna;

   // List of branches
   TBranch        *b_photonStep;   //!
   TBranch        *b_extrapatstep;   //!
   TBranch        *b_ana_step;   //!
   TBranch        *b_phltJetHI;   //!
   TBranch        *b_pcollisionEventSelection;   //!
   TBranch        *b_pHBHENoiseFilter;   //!
   TBranch        *b_phiEcalRecHitSpikeFilter;   //!
   TBranch        *b_hltAna;   //!

   Skim(TTree *tree=0);
   virtual ~Skim();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Skim_cxx
Skim::Skim(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
  /*    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet80/Track8_Jet19_STARTHI53_LV1/merged_300k_2/HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0.root");
      if (!f) {
         f = new TFile("/mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet80/Track8_Jet19_STARTHI53_LV1/merged_300k_2/HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0.root");
      }
      tree = (TTree*)gDirectory->Get("HltTree");
*/
   }
   Init(tree);
}

Skim::~Skim()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Skim::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Skim::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Skim::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normaly not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("photonStep", &photonStep, &b_photonStep);
   fChain->SetBranchAddress("extrapatstep", &extrapatstep, &b_extrapatstep);
   fChain->SetBranchAddress("ana_step", &ana_step, &b_ana_step);
   fChain->SetBranchAddress("phltJetHI", &phltJetHI, &b_phltJetHI);
   fChain->SetBranchAddress("pcollisionEventSelection", &pcollisionEventSelection, &b_pcollisionEventSelection);
   fChain->SetBranchAddress("pHBHENoiseFilter", &pHBHENoiseFilter, &b_pHBHENoiseFilter);
   fChain->SetBranchAddress("phiEcalRecHitSpikeFilter", &phiEcalRecHitSpikeFilter, &b_phiEcalRecHitSpikeFilter);
   fChain->SetBranchAddress("hltAna", &hltAna, &b_hltAna);
   Notify();
}

Bool_t Skim::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Skim::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Skim::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Skim_cxx
