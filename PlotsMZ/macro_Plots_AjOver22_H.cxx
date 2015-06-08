
#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TTree.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TLatex.h"
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

int macro_Plots_AjOver22_H(int i=0)
{
  TFile  *finbg, *fout_inc ;
  TFile *fin = TFile::Open("Data2011_AjUnder22_Only13.root");  // The input file
  fout_inc = new TFile("PbPb_Inclusive_Correlations.root","RECREATE"); // The output file

  const int nCBins = 4;   //Centrality
  const int nPtBins = 1;  // Hardest Jet range
  const int nTrkPtBins = 5; //Associated tracks range

  float me00[nCBins][nPtBins][nTrkPtBins]; // used to save the mixing at (Deta,Dphi) = (0,0) (the output of fitting)
  float SideBandRatio[nCBins][nPtBins][nTrkPtBins]; // used to save the mixing at (Deta,Dphi) = (0,0) (the output of fitting)
  int NumberOfDijets = 0;
  int NEtaBins = 0;
  float PhiBinWidth = 0.0;
  float EtaBinWidth = 0.0;
  //	int me00bin;

  float PtBins[nPtBins+1] = {100, 300};
  int NDijets[nCBins];
  TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};
  TLatex *centtex, *pttex; 

  float CBins[nCBins+1] = {0, 20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%","Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};
  //MZS
  float MaxAj = 0; 
		
  float AjBins[1] = {0.11};
  //TString AjBin_str[1] = {"Aj < 0.22"};
  TString AjBin_str[1] = {"Aj > 0.22"};
  //TString AjBin_labels[1] = {"Aj < 0.22"};
  TString AjBin_labels[1] = {"Aj > 0.22"};
  //MZE
  float TrkPtBins[nTrkPtBins+1] = {1, 2, 3, 4, 8, 999};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt999" };
  TString TrkPtBin_labels[nTrkPtBins] = {"1<pT<2","2<pT<3","3<pT<4","4<pT<8","pT>8"};

  float JetLabel[2] = {1, 2};
  TString JetLabel_str[2] = {"Leading Jet","Leading Jet"};
  //TString JetLabel_str[2] = {"    Leading","    Leading"};
  TString JetLabel_labels[2] = {"Leading Jet","Leading Jet"};
  //TString JetLabel_labels[2] = {"   Leading","   Leading"};


  float left_level, right_level;

	
  TCanvas *cc1 = new TCanvas("cc1", "CC1_me_projeta_HT", 0, 0, 1500, 1500);	
  TCanvas *cc2 = new TCanvas("cc2", "CC2_me_projphi_HT", 0, 0, 1500, 1500);	
  TCanvas *cc3 = new TCanvas("cc3", "CC3_SignalNoScale_HT", 0, 0, 1500, 1500);	
  TCanvas *cc4 = new TCanvas("cc4", "CC4_JetTrackOverNormalizedME_etaProj", 0, 0, 1500, 1500);	
  TCanvas *cc5 = new TCanvas("cc5", "CC5_ME_HT", 0, 0, 1500, 1500);	
  TCanvas *cc6 = new TCanvas("cc6", "CC6_JetTrackOverNormalizedME_etaProj_LBoverRB", 0, 0, 1500, 1500);	
  TCanvas *cc7 = new TCanvas("cc7", "CC7_EtaBandSymmetry", 0, 0, 1500, 1500);	
  TCanvas *cc8 = new TCanvas("cc8", "CC8_EtaBandSymmetry_HT", 0, 0, 1500, 1500);	
  TCanvas *cc9 = new TCanvas("cc9", "CC9_EtaBandSymmetry_HT", 0, 0, 1500, 1500);	

	
  cc1->Divide(4,4,0.0001,0.0001);   // 4 Raw by 4 Columns and the separation is 0.1 and 0.1
  cc2->Divide(4,4,0.0001,0.0001);
  cc3->Divide(4,4,0.0001,0.0001);
  cc4->Divide(4,4,0.0001,0.0001);
  cc5->Divide(4,4,0.0001,0.0001);
  cc6->Divide(4,4,0.0001,0.0001);
  cc7->Divide(4,4,0.0001,0.0001);
  cc8->Divide(4,4,0.0001,0.0001);
  cc9->Divide(4,4,0.0001,0.0001);
 
	

  TF1 *fa = new TF1("fa","[0]",-5,5);
  TF1 *flat_fit = new TF1("flat_fit", "[0] + x - x " );	
 
  cout << "the default parameter = " <<  i <<  endl;
  int Logical = 1;
  int lbin, rbin;

  TString  datalabel;

  gStyle->SetOptStat(0);   
  gStyle->SetPadBottomMargin(0.05);
  gStyle->SetPadTopMargin   (0.05);
  gStyle->SetPadLeftMargin  (0.05);
  gStyle->SetPadRightMargin (0.05);
  gStyle->SetTitleSize(0.18,"xyz");

  TH2D* hJetTrackSignalBackgroundLeading[nCBins][nPtBins][nTrkPtBins];  // The signal
  TH2D* hJetTrackMELeading[nCBins][nPtBins][nTrkPtBins];                // The mixed events
  TH2D *yield_inc[nCBins][nPtBins][nTrkPtBins];      	       // used to save a clone of hJetTrackSignalBackgroundLeading[][][] and then become the Histo of Signal over ME
  TH2D *yield_inc_copy[nCBins][nPtBins][nTrkPtBins];      	       // used to save a clone of hJetTrackSignalBackgroundLeading[][][] and then become the Histo of Signal over ME
  //TH2D *yield_inc_copy2[nCBins][nPtBins][nTrkPtBins];      	       // used to save a clone of yield_inc for separat studies
  TH2D* Correlation[nCBins][nPtBins][nTrkPtBins];


  TH1D *yield_inc_SBL[nCBins][nPtBins][nTrkPtBins];              // The Y projection of the yield (phi) at the left range
  TH1D *yield_inc_SBR[nCBins][nPtBins][nTrkPtBins];		// The Y projection of the yield (phi) at the right range 
  TH1D *yield_inc_Signal[nCBins][nPtBins][nTrkPtBins];		// The Y projection of the yield (phi) at the center (the signal) 
  TH1D *yield_inc_SBL_ForRatio[nCBins][nPtBins][nTrkPtBins];		// The Y projection of the yield (phi) at the center (the signal) 

  TH1D *yield_corrected_inc_SBL_corrected[nCBins][nPtBins][nTrkPtBins];
  TF1 *bg_fit[nCBins][nPtBins][nTrkPtBins]; 
  //The Y projection of signal/ME, Left band

	

  char *histname = new char[10];
  int xx = 0;
  double A[nCBins][nPtBins][nTrkPtBins] ;
  double B[nCBins][nPtBins][nTrkPtBins] ;
  double C[nCBins][nPtBins][nTrkPtBins] ;

  sprintf(histname, "h_x_%d",xx);
  //	TH1D *yield_corrected_inc_SBL_corrected[nCBins][nPtBins][nTrkPtBins]=  new TH1D("yield_corrected_inc_SBL_corrected"[nCBins][nPtBins][nTrkPtBins], "Mass;mass (kg)", 100, 0, 5);
  TH1D *yield_corrected_inc_SBR_corrected[nCBins][nPtBins][nTrkPtBins];             //The Y projection of signal/ME, Right band

  TH1D *yield_inc_Signal_ForBinCheck[nCBins][nPtBins][nTrkPtBins];		// copies from yield_inc_Signal to make studies on the side bands 

  TH1D *me_projeta[nCBins][nPtBins][nTrkPtBins];
  TH1D *yield_inc_proj[nCBins][nPtBins][nTrkPtBins];
  TH1D *yield_inc_proj_copy[nCBins][nPtBins][nTrkPtBins];
  TH1D *yield_inc_proj_signal[nCBins][nPtBins][nTrkPtBins];
  TH1D *test_me_inc_proj[nCBins][nPtBins][nTrkPtBins];
  TH1D *ProEta;
  TH1D *me_projphi[nCBins][nPtBins][nTrkPtBins];
  TH1D *EtaProj[nCBins][nPtBins][nTrkPtBins];
  TH1D *PhiProj_L[nCBins][nPtBins][nTrkPtBins];
  TH1D *PhiProj_R[nCBins][nPtBins][nTrkPtBins];
  TH1D *PhiProj_RPlusL[nCBins][nPtBins][nTrkPtBins];
	
  TH1D *EtaProjDummy[nCBins][nPtBins][nTrkPtBins];
  TH1D *Signal[nCBins][nPtBins][nTrkPtBins];
  TH1D *SignalMinusRL[nCBins][nPtBins][nTrkPtBins];
	
  TH1D *Aj[nCBins];	

  TString desc = "hists";

  for (int ibin=0;ibin<nCBins;ibin++){
    Aj[ibin] = (TH1D*)fin->Get((TString)(desc + "_Aj"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]))->Clone((TString) ("Aj_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]));    
    for (int ibin2=0;ibin2<nPtBins;ibin2++){ 
      for (int ibin3=0;ibin3<nTrkPtBins-1;ibin3++){
	cout<<ibin<<" "<<ibin2<<" "<<ibin3<<endl;
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	
	hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

      } // end of ibin3
    } //end of ibin2

    cout << "after 2 loops" << endl;
    cout << " Aj based Number of Dijets = " << Aj[ibin]->GetEntries() << " For Centrality bin = " << ibin+1 << endl;
    int AjBinMax = Aj[ibin]->GetMaximumBin();
    float AjMax =  Aj[ibin]->GetXaxis()->GetBinUpEdge(AjBinMax);
    cout << " Aj for this sampel = " << AjMax << endl;

    // Normalizing

    for(int ibin2=0; ibin2<nPtBins;ibin2++){
      for (int ibin3=0; ibin3< 4; ibin3++){

	lbin = hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->FindBin(-1.5);
	rbin = hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->FindBin(1.5);
	//		printf("the line number is: %d \n \n", __LINE__);

	cout << "hJetTrackSignalBackground " << ibin << " " << ibin2 << " " << ibin3 << endl;
	//	hJetTrackSignalBackground[ibin][ibin2][ibin3]->Draw("LEGO2");

	cout << "Fitting the central part of eta" << endl;

	fout_inc->cd();
	//cout << "Canvas number = " << 4*(ibin3+1)-ibin << endl;
	//me_proj_canvas->cd(4*(ibin3+1)-ibin);
				
	yield_inc[ibin][ibin2][ibin3] =  (TH2D*) hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Clone((TString)("Yield_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

	yield_inc_copy[ibin][ibin2][ibin3] =  (TH2D*) hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Clone((TString)("Yield_copy"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));


	me_projeta[ibin][ibin2][ibin3] =  hJetTrackMELeading[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_eta%d%d%d",ibin,ibin2,ibin3),1,100); //what is the 1,100 for?
	me_projphi[ibin][ibin2][ibin3] =  hJetTrackMELeading[ibin][ibin2][ibin3]->ProjectionY(Form("me_proj_phi%d%d%d",ibin,ibin2,ibin3),1,100); //what is the 1,100 for?
		

	cout << "Trying to get Number of eta bins = " << endl;
	cout << me_projeta[ibin][ibin2][ibin3]->GetNbinsX() << endl;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	cc1->cd(4*(ibin3+1)-ibin);
	//me_projeta[ibin][ibin2][ibin3]->Draw() ;
	me_projeta[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("Delta eta");
	NEtaBins =  me_projeta[ibin][ibin2][ibin3]->GetNbinsX() ;
	me_projeta[ibin][ibin2][ibin3] -> Fit("pol0","","",-.2,.2);
	me00[ibin][ibin2][ibin3] = me_projeta[ibin][ibin2][ibin3]->GetFunction("pol0")->GetParameter(0);		

	me_projeta[ibin][ibin2][ibin3]->SetAxisRange(-2, 2,"X");
	me_projeta[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");
	//me_projeta[ibin][ibin2][ibin3]->DrawCopy() ;
	me_projeta[ibin][ibin2][ibin3]->Draw() ;

	if(ibin<3){
	  TLatex *centtex = new TLatex(0.15,0.85,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	  TLatex *pttex = new TLatex(0.15,0.80,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	  TLatex *Ajtex = new TLatex(0.76,0.85,AjBin_labels[0]);
	  Ajtex->SetNDC();
	  Ajtex->Draw();
	  TLatex *Jettex = new TLatex(0.70, 0.80, JetLabel_labels[0]);
	  Jettex->SetNDC();
	  Jettex->Draw();

	}
	if(ibin==3){
	  TLatex *centtex = new TLatex(0.15,0.85,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	  TLatex *pttex = new TLatex(0.15,0.80,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	  TLatex *Ajtex = new TLatex(0.76,0.85,AjBin_labels[0]);
	  Ajtex->SetNDC();
	  Ajtex->Draw();
	  TLatex *Jettex = new TLatex(0.70, 0.80, JetLabel_labels[0]);
	  Jettex->SetNDC();
	  Jettex->Draw();
				

	}
		
	cc2->cd(4*(ibin3+1) - ibin);
	me_projeta[ibin][ibin2][ibin3] =  hJetTrackMELeading[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_eta%d%d%d",ibin,ibin2,ibin3),1,100); //what is the 1,100 for?

	me_projeta[ibin][ibin2][ibin3]->Scale (100./me00[ibin][ibin2][ibin3]);

	hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(100./me00[ibin][ibin2][ibin3]);

	//me_projeta[ibin][ibin2][ibin3]->Draw();
	me_projeta[ibin][ibin2][ibin3]->Draw();

	//##//
	{
	  TLatex *centtex = new TLatex(0.15,0.85,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	  TLatex *pttex = new TLatex(0.15,0.80,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	  TLatex *Ajtex = new TLatex(0.76,0.85,AjBin_labels[0]);
	  Ajtex->SetNDC();
	  Ajtex->Draw();
	  TLatex *Jettex = new TLatex(0.70, 0.80, JetLabel_labels[0]);
	  Jettex->SetNDC();
	  Jettex->Draw();

	}


	//##//
	//##//							
	cc3->cd(4*(ibin3+1)-ibin);

	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->CenterTitle();
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->CenterTitle();
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitleOffset(1.6);
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitleOffset(1.6);
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Draw("LEGO2");
	
	if(ibin<3){
	  TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	  TLatex *pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	  TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
	  Ajtex->SetNDC();
	  Ajtex->Draw();
	  TLatex *Jettex = new TLatex(0.70, 0.90, JetLabel_labels[0]);
	  Jettex->SetNDC();
	  Jettex->Draw();


	}
	if(ibin==3){
	  TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	  TLatex *pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	  TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
	  Ajtex->SetNDC();
	  Ajtex->Draw();
	  TLatex *Jettex = new TLatex(0.70, 0.90, JetLabel_labels[0]);
	  Jettex->SetNDC();
	  Jettex->Draw();


	}

//return 3;
	cc4->cd(4*(ibin3+1)-ibin);

	PhiBinWidth =  hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->GetBinWidth(1) ;
	EtaBinWidth =  hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->GetBinWidth(1) ;
				

	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() * PhiBinWidth*EtaBinWidth));
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-3.,3.);	
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> CenterTitle();
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> CenterTitle();
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(1.6);
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(1.6);

	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Draw("LEGO2");	

	{
	  TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	  TLatex *pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	  TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
	  Ajtex->SetNDC();
	  Ajtex->Draw();
	  TLatex *Jettex = new TLatex(0.70, 0.85, JetLabel_labels[0]);
	  Jettex->SetNDC();
	  Jettex->Draw();


	}
	
//return 4;
	cc5->cd(4*(ibin3+1)-ibin);
		
	hJetTrackMELeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
	hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
	hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis() -> CenterTitle();
	hJetTrackMELeading[ibin][ibin2][ibin3]->GetYaxis() -> CenterTitle();
	hJetTrackMELeading[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(1.6);
	hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(1.6);


	hJetTrackMELeading[ibin][ibin2][ibin3]->Draw("LEGO2");

	{
	  TLatex *centtex = new TLatex(0.15,0.85,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	  TLatex *pttex = new TLatex(0.15,0.80,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	  TLatex *Ajtex = new TLatex(0.76,0.85,AjBin_labels[0]);
	  Ajtex->SetNDC();
	  Ajtex->Draw();
	  TLatex *Jettex = new TLatex(0.70, 0.80, JetLabel_labels[0]);
	  Jettex->SetNDC();
	  Jettex->Draw();

	}

//return 5;
	cc6->cd(4*(ibin3+1)-ibin);
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Divide(hJetTrackMELeading[ibin][ibin2][ibin3]);
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Draw("LEGO2");	

	{
	  TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	  TLatex *pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	  TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
	  Ajtex->SetNDC();
	  Ajtex->Draw();
	  TLatex *Jettex = new TLatex(0.70, 0.85, JetLabel_labels[0]);
	  Jettex->SetNDC();
	  Jettex->Draw();

	}
//return 6;

	cc7->cd(4*(ibin3+1)-ibin);
	
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	

	yield_inc_proj[ibin][ibin2][ibin3] = (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionX((Form("ResultEta%d%d%d",ibin,ibin2,ibin3)),1, 25);	
	yield_inc_proj[ibin][ibin2][ibin3]->Rebin(2);
	yield_inc_proj[ibin][ibin2][ibin3]->Scale(EtaBinWidth/2.);
			
	yield_inc_proj[ibin][ibin2][ibin3]->Fit("flat_fit","0","",-3.,-1.5);
	left_level = flat_fit->GetParameter(0);
	yield_inc_proj[ibin][ibin2][ibin3]->Fit("flat_fit","0","",1.5,3.);
	right_level =  flat_fit->GetParameter(0);
	bg_fit[ibin][ibin2][ibin3] = new TF1("bg_fit%d%d%d","[0]+x-x",-3.,3.);
	double level = (left_level+right_level)/2.;
	bg_fit[ibin][ibin2][ibin3]->SetParameter(0,(left_level+right_level)/2.);
	yield_inc_proj[ibin][ibin2][ibin3]->SetAxisRange(-3., 3., "X");

					
	yield_inc_proj[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");
	if(ibin==0&&ibin3==0)	yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(80.);
	if(ibin==1&&ibin3==0)	yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(50.);
	if(ibin==2&&ibin3==0)	yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(20.);
	yield_inc_proj[ibin][ibin2][ibin3]->Draw();
	bg_fit[ibin][ibin2][ibin3]->Draw("same");

						
	{
	  TLatex *centtex = new TLatex(0.15,0.85,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	  TLatex *pttex = new TLatex(0.15,0.8,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	  TLatex *Ajtex = new TLatex(0.75,0.85,AjBin_labels[0]);
	  Ajtex->SetNDC();
	  Ajtex->Draw();
	  TLatex *Jettex = new TLatex(0.70, 0.80, JetLabel_labels[0]);
	  Jettex->SetNDC();
	  Jettex->Draw();

	}

								
//return 7;
	cc8->cd((4*(ibin3+1)-ibin));

	Signal[ibin][ibin2][ibin3] = (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionY(Form("Phi_Signal%d%d%d",ibin,ibin2,ibin3),36, 65);		
	Signal[ibin][ibin2][ibin3]->SetMarkerStyle(20);
	Signal[ibin][ibin2][ibin3]->SetMarkerColor(kRed);
	Signal[ibin][ibin2][ibin3]->Rebin(2);
	Signal[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#phi");
	Signal[ibin][ibin2][ibin3]->Scale(EtaBinWidth/2.);

	Signal[ibin][ibin2][ibin3]->Draw();


	PhiProj_L[ibin][ibin2][ibin3] = (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionY(Form("ResultPhi_L%d%d%d",ibin,ibin2,ibin3),21, 35);
	//PhiProj_L[ibin][ibin2][ibin3]->Draw("same");

	PhiProj_R[ibin][ibin2][ibin3] =  (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionY(Form("ResultPhi_R%d%d%d",ibin,ibin2,ibin3),66, 80);	
	//PhiProj_R[ibin][ibin2][ibin3]->Draw("same");

	PhiProj_R[ibin][ibin2][ibin3]->Add(PhiProj_L[ibin][ibin2][ibin3], 1);
	PhiProj_R[ibin][ibin2][ibin3]->SetMarkerStyle(20);
	PhiProj_R[ibin][ibin2][ibin3]->SetMarkerColor(kBlue);
	PhiProj_R[ibin][ibin2][ibin3]->Rebin(2);
	PhiProj_R[ibin][ibin2][ibin3]->Scale(EtaBinWidth/2.);
	PhiProj_R[ibin][ibin2][ibin3]->Draw("same");


	{
	  TLatex *centtex = new TLatex(0.15,0.85,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	  TLatex *pttex = new TLatex(0.15,0.80,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	  TLatex *Ajtex = new TLatex(0.75,0.85,AjBin_labels[0]);
	  Ajtex->SetNDC();
	  Ajtex->Draw();
	  TLatex *Jettex = new TLatex(0.70, 0.80, JetLabel_labels[0]);
	  Jettex->SetNDC();
	  Jettex->Draw();

	}
//return 8;			

	cc9->cd(4*(ibin3+1)-ibin);
	SignalMinusRL[ibin][ibin2][ibin3] = (TH1D*)Signal[ibin][ibin2][ibin3]->Clone((TString)("signal_MinusRL"+CBin_strs[ibin]+"_"+TrkPtBin_strs[ibin3]));
	SignalMinusRL[ibin][ibin2][ibin3]->Add(	PhiProj_R[ibin][ibin2][ibin3], -1);
	SignalMinusRL[ibin][ibin2][ibin3]->SetAxisRange(-1.5, 1.5,"X");
	if(ibin3==0)	SignalMinusRL[ibin][ibin2][ibin3]->SetMinimum(-2.);
	if(ibin3==0)	SignalMinusRL[ibin][ibin2][ibin3]->SetMaximum(8.);
	SignalMinusRL[ibin][ibin2][ibin3]->Draw();
				
	{
	  TLatex *centtex = new TLatex(0.15,0.85,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	  TLatex *pttex = new TLatex(0.15,0.8,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	  TLatex *Ajtex = new TLatex(0.75,0.85,AjBin_labels[0]);
	  Ajtex->SetNDC();
	  Ajtex->Draw();
	  TLatex *Jettex = new TLatex(0.70, 0.80, JetLabel_labels[0]);
	  Jettex->SetNDC();
	  Jettex->Draw();

	}



      }
    }
  }

  cc1->SaveAs("cc1_ME_EtaProj_HT.png");
  cc2->SaveAs("cc2_MEOverME00_HT.png");
  cc3->SaveAs("cc3_SignalNoScale_HT.png");
  cc4->SaveAs("cc4_S_HT.png");
  cc5->SaveAs("cc5_ME_HT.png");
  cc6->SaveAs("cc6_SOverME.png");
  cc7->SaveAs("cc7_RightBandLeftBandComparison.png");
  cc8->SaveAs("cc8_PhiProjectionofSignalandBands.png");
  cc9->SaveAs("cc9_resultOnPhi.png");
  return i;
}                                           

