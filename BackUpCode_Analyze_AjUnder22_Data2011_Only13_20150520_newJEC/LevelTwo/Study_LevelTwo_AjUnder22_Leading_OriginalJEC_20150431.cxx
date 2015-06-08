
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
	#include "HIN_14_016_functions.h" // used for rebinning
	using namespace std;

	int Study_LevelTwo_AjUnder22_Leading_OriginalJEC_20150431(double i=0, double j = 0)
	{
	double eta_limit, phi_limit;
	eta_limit = i;
	phi_limit = j;
	cout << "eta_limit = " << eta_limit << " phi_limit = " << phi_limit << endl;	
	TFile  *finbg, *fout_inc ;

	//TFile *fin = TFile::Open("Data2011_trkPtCut1_All_AJ11.root");  // The input file
	//TFile *fin = TFile::Open("Data2011_All_AJ22.root");  // The input file
	//TFile *fin = TFile::Open("Data2011_All_AJOver22.root");  // The input file
	//TFile *fin = TFile::Open("Data2011_AjUnder22.root");  // The input file
	//TFile *fin = TFile::Open("Data2011_AjUnder22_Only13.root");  // The input file
	TFile *fin = TFile::Open("/home/mzakaria3/Documents/Research/MyGitProject/PlotsMZ/Data2011_AjUnder22_Only13_20150416.root");  // The input file
	TFile *fin_LevelOne = TFile::Open("~/Documents/Research/MyGitProject/PlotsMZ/BackUpCode_Analyze_AjUnder22_Data2011_Only13_20150520_newJEC/LevelOne/LevelOne_PbPb_AjUnder22_Correlations_20150520.root");  // The input file
	//TFile *fin = TFile::Open("FromH/PbPb_Leading_Correlations.root");  // The input file
	//TFile *fin = TFile::Open("Data2011_All_AJ33.root");  // The input file
	//TFile *fin = TFile::Open("/home/mzakaria3/Documents/Research/MyGitProject/PlotsMZ/Data2011_trkPtCut1_All_AJ11.root");  // The input file
	//	TFile *fin = TFile::Open("/home/mzakaria3/Documents/Research/MyGitProject/PlotsMZ/Data2011_trkPtCut1_All_Aj11_20150318.root");  // The input file
	//TFile *fin = TFile::Open("/home/mzakaria3/Documents/Research/MyGitProject/PlotsMZ/Data2011_trkPtCut1_All_Aj11_20150319.root");  // The input file

	fout_inc = new TFile("LevelTwo_PbPb_Under22_Correlations_20150520.root","RECREATE"); // The output file


	const int nCBins = 4;   //Centrality
	const int nPtBins = 1;  // Hardest Jet range
	const int nTrkPtBins = 5; //Associated tracks range

	float me00[nCBins][nPtBins][nTrkPtBins]; // used to save the mixing at (Deta,Dphi) = (0,0) (the output of fitting)
	float SideBandRatio[nCBins][nPtBins][nTrkPtBins]; // used to save the mixing at (Deta,Dphi) = (0,0) (the output of fitting)
	int NumberOfDijets = 0;
	int NEtaBins = 0;
	float PhiBinWidth = 0.0;
	float EtaBinWidth = 0.0;
	float left_level = 0, right_level = 0, left_level_error = 0, right_level_error = 0;
	//	int me00bin;

	float PtBins[nPtBins+1] = {100, 300};
	int NDijets[nCBins] ;
	TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};
	TLatex *centtex, *pttex; 

	float CBins[nCBins+1] = {0, 20, 60, 100, 200};
	TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
	TString CBin_labels[nCBins] = {"Cent. 0-10%","Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};
	//MZS
	float MaxAj = 0; 
		
	float AjBins[1] = {0.11};
	//TString AjBin_str[1] = {"Aj < 0.22"};
	TString AjBin_str[1] = {"Aj < 0.22"};
	//TString AjBin_labels[1] = {"Aj < 0.22"};
	TString AjBin_labels[3] = {"Aj < 0.22", "No Aj Cut", "Aj > 0.22"};
	//MZE
	float TrkPtBins[nTrkPtBins+1] = {1, 2, 3, 4, 8, 999};
	TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt999" };
	TString TrkPtBin_labels[nTrkPtBins] = {"1<p^{assoc.}_{T}<2 GeV/c","2<p^{assoc.}_{T}<3 GeV/c","3<p^{assoc.}_{T}<4 GeV/c","4<p^{assoc.}_{T}<8 GeV/c","p^{assoc.}_{T}>8GeV/c"};
	TString COM_label[1] = {"#sqrt{S_{NN}}=2.76 TeV"};
	TString CMS_label[1] = {"CMS Preliminary"};
	TString Pb_label[1] = {"PbPb 169 #mu b^{-1}"};
	TString eta_label[1] = {"|#eta| < 1.6"};
	TString pT1_label[1] = {"p_{T,jet}>120 GeV/c"};

	float JetLabel[2] = {1, 2};
	TString JetLabel_str[2] = {"Leading Jet","Leading Jet"};
	//TString JetLabel_str[2] = {"    Leading","    Leading"};
	TString JetLabel_labels[2] = {"Leading Jet","Leading Jet"};
	//TString JetLabel_labels[2] = {"   Leading","   Leading"};




	TCanvas *dd1 = new TCanvas("dd1", "DD1_DPhiProjSandB", 0, 0, 1500, 1500);	
	TCanvas *dd2 = new TCanvas("dd2", "DD2_DPhiProjS_Minu_B", 0, 0, 1500, 1500);	
	TCanvas *dd3 = new TCanvas("dd3", "DD3_DPhiProjS_Minu_B_Rebinned", 0, 0, 1500, 1500);	
	TCanvas *cc8 = new TCanvas("cc8", "CC7_EtaBandSymmetry", 0, 0, 1500, 1500);	
	TCanvas *cc9 = new TCanvas("cc9", "CC8_EtaBandSymmetry_HT", 0, 0, 1500, 1500);	
	TCanvas *cc10 = new TCanvas("cc10", "CC10_EtaBandSymmetry_HT", 0, 0, 1500, 1500);	
	TCanvas *cc11 = new TCanvas("cc11", "CC11_EtaBandSymmetry_HT", 0, 0, 1500, 1500);	
	TCanvas *cc12 = new TCanvas("cc12", "CC12_CC6j_LBoverRB_magnified", 0, 0, 1500, 1500);	
	TCanvas *cc13 = new TCanvas("cc13", "CC13_SignalLBoverRB_magnified", 0, 0, 1500, 1500);	
	TCanvas *cc14 = new TCanvas("cc14", "CC14_ME", 0, 0, 1500, 1500);	
	TCanvas *cc15 = new TCanvas("cc15", "CC15_ME_EtaProj", 0, 0, 1500, 1500);	
	TCanvas *cc16 = new TCanvas("cc16", "CC16_ME_EtaProj_ScaledTo1", 0, 0, 1500, 1500);	
	TCanvas *cc17 = new TCanvas("cc17", "CC17_MEoverM0", 0, 0, 1500, 1500);	
	TCanvas *cc18 = new TCanvas("cc18", "CC18_MEoverM0_ByJets", 0, 0, 1500, 1500);	
	TCanvas *cc19 = new TCanvas("cc19", "CC19_MEoverM0_ByJets_EtaProj_PhiNarrow", 0, 0, 1500, 1500);	
	TCanvas *cc20 = new TCanvas("cc20", "CC20_MEoverM0_ByJets_PhiProh_EtaSide", 0, 0, 1500, 1500);	
	TCanvas *cc21 = new TCanvas("cc21", "CC21_Yield", 0, 0, 1500, 1500);	
	TCanvas *cc22 = new TCanvas("cc22", "CC22_S_Scaled_NotPerEtaNotPerPhi", 0, 0, 1500, 1500);	
	TCanvas *cc23 = new TCanvas("cc23", "CC23_S_Scaled_overMEoverM0_NotScaled_NotPerEtaNotPerPhi", 0, 0, 1500, 1500);	
	TCanvas *cc24 = new TCanvas("cc24", "CC24_S_Scaled_overMEoverM0_NotScaled_PerEtaPerPhi_BkgSubtracted", 0, 0, 1500, 1500);	
	TCanvas *cc25 = new TCanvas("cc25", "CC25_S_Scaled_overMEoverM0_NotScaled_PerEtaPerPhi_etaSideBand", 0, 0, 1500, 1500);	
	TCanvas *cc26 = new TCanvas("cc26", "CC26_S_Scaled_overMEoverM0_NotScaled_PerEtaPerPhi_MinusEtaSideBand_ProjectionY", 0, 0, 1500, 1500);	
	TCanvas *cc27 = new TCanvas("cc27", "CC27_Rebinning", 0, 0, 1500, 1500);
	//cc27->SetRightMargin(1.001);
	//cc27->SetBorderMode(1);

	TCanvas *cc28 = new TCanvas("cc28", "CC28_SubtractBeforeFilling", 0, 0, 1500, 1500);	


	//	c1->Divide(4,4,0.0,0.0);
	//	c6->Divide(4,4,0.0,0.0);
	//	c8->Divide(4,4,0.0,0.0);
	//	c9->Divide(4,4,0.0,0.0);
	//	c10->Divide(4,4,0.0,0.0);


	dd1->Divide(4,4,0.0001,0.0001);
	dd2->Divide(4,4,0.0001,0.0001);
	dd3->Divide(4,4,0.0001,0.0001);
	cc9->Divide(4,4,0.0001,0.0001);
	cc10->Divide(4,4,0.0001,0.0001);
	//cc11->Divide(4,4,0.0001,0.0001);

	cc11->Divide(4,4,0.000,0.000);
	cc12->Divide(4,4,0.0001,0.0001);
	cc13->Divide(4,4,0.0001,0.0001);
	cc14->Divide(4,4,0.0001,0.0001);
	cc15->Divide(4,4,0.0001,0.0001);
	cc16->Divide(4,4,0.0001,0.0001);
	cc17->Divide(4,4,0.0001,0.0001);
	cc18->Divide(4,4,0.0001,0.0001);
	cc19->Divide(4,4,0.0001,0.0001);
	cc20->Divide(4,4,0.0001,0.0001);
	cc21->Divide(4,4,0.0001,0.0001);
	cc22->Divide(4,4,0.0001,0.0001);
	cc23->Divide(4,4,0.0001,0.0001);
	cc24->Divide(4,4,0.0001,0.0001);
	cc25->Divide(4,4,0.0001,0.0001);
	cc26->Divide(4,4,0.0000,0.0000);
	//cc26_4->SetRightMargin(0.01);
	cc27->Divide(4,4,0.0000,0.0000);
if (1) {
    cc26->cd(4);
    gPad->SetRightMargin(0.01);
    cc26->cd(8);
    gPad->SetRightMargin(0.01);
    cc26->cd(12);
    gPad->SetRightMargin(0.01);
    cc26->cd(16);
    gPad->SetRightMargin(0.01);
}	

if (1) {
    cc26->cd(1);
    gPad->SetLeftMargin(0.15);
    cc26->cd(5);
    gPad->SetLeftMargin(0.15);
    cc26->cd(9);
    gPad->SetLeftMargin(0.15);
    cc26->cd(13);
    gPad->SetLeftMargin(0.15);
}	
if (1) {
    cc27->cd(4);
    gPad->SetRightMargin(0.01);
    cc27->cd(8);
    gPad->SetRightMargin(0.01);
    cc27->cd(12);
    gPad->SetRightMargin(0.01);
    cc27->cd(16);
    gPad->SetRightMargin(0.01);
}	

if (1) {
    cc27->cd(1);
    gPad->SetLeftMargin(0.15);
    cc27->cd(5);
    gPad->SetLeftMargin(0.15);
    cc27->cd(9);
    gPad->SetLeftMargin(0.15);
    cc27->cd(13);
    gPad->SetLeftMargin(0.15);
}	

	
	cc28->Divide(4,4,0.0001,0.0001);


	TF1 *fa = new TF1("fa","[0]",-5,5);
	TF1 *flat_fit = new TF1("flat_fit", "[0] + x - x " );	
	//TCanvas *me_proj_canvas = new TCanvas("me_proj_canvas","me_proj_canvas",0,0,1500,1600);
	//me_proj_canvas->Divide(4,4,0.0000,0.0000);

	//TCanvas *result_proj_canvas = new TCanvas("result_proj_canvas","",0,0,1500,1600);
	//result_proj_canvas->Divide(4,4,0.0000,0.0000);

	//TCanvas *me_proj_sideBands = new TCanvas("me_proj_sideBands","",0,0,1500,1600);
	//me_proj_sideBands->Divide(4,4,0.0000,0.0000);

	//	TCanvas *result_proj_OnlyLeading = new TCanvas("result_proj_OnlyLeading","result_proj_OnlyLeading",0,0,1500,1600);
	//	result_proj_OnlyLeading->Divide(4,4,0.0000,0.0000);

	cout << "the default parameter = " <<  i <<  endl;
	int Logical = 1;
	int lbin, rbin;
	double llimiteta, rlimiteta;
	//if(!Logical)
	//  printf("Not logical value at line number %d in file %s\n", __LINE__, __FILE__);  // used to print the line (debugging)

	TString  datalabel;

	gStyle->SetOptStat(0);   
	//gStyle->SetPadBottomMargin(0.001);
	//gStyle->SetPadTopMargin   (0.001);
	//gStyle->SetPadLeftMargin  (0.001);
	//gStyle->SetPadRightMargin (0.001);
	gStyle->SetTitleSize(0.06,"xyz");
	gStyle->SetLabelSize(0.06);
				//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(0.8);
				//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(0.8);
				//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);

	TH2D *Scaled_Normalized_Yield_Inv[nCBins][nPtBins][nTrkPtBins];	

	TH2D* hJetBackgroundLeading[nCBins][nPtBins][nTrkPtBins];  // The signal
	TH2D* hJetTrackSignalBackgroundLeading[nCBins][nPtBins][nTrkPtBins];  // The signal
	TH2D* hJetTrackSignalBackgroundLeading_cc26[nCBins][nPtBins][nTrkPtBins];  // The signal
	TH2D* hJetTrackMELeading[nCBins][nPtBins][nTrkPtBins];                // The mixed events
	TH2D* hJetTrackMELeading_cc26[nCBins][nPtBins][nTrkPtBins];                // The mixed events
	TH2D* hJetTrackMELeading_cc5[nCBins][nPtBins][nTrkPtBins];                // The mixed events
	TH2D *yield_inc[nCBins][nPtBins][nTrkPtBins];      	       // used to save a clone of hJetTrackSignalBackgroundLeading[][][] and then become the Histo of Signal over ME
	TH2D *yield_inc_copy[nCBins][nPtBins][nTrkPtBins];      	       // used to save a clone of hJetTrackSignalBackgroundLeading[][][] and then become the Histo of Signal over ME
	//TH2D *yield_inc_copy2[nCBins][nPtBins][nTrkPtBins];      	       // used to save a clone of yield_inc for separat studies
	TH2D* Correlation[nCBins][nPtBins][nTrkPtBins];
	TH2D* result[nCBins][nPtBins][nTrkPtBins];
	TH2D* result_cc26[nCBins][nPtBins][nTrkPtBins];
	TH2D* background[nCBins][nPtBins][nTrkPtBins];
	TH2D* background_cc26[nCBins][nPtBins][nTrkPtBins];
	TH1D* background_left[nCBins][nPtBins][nTrkPtBins];
	TH1D* background_left_cc26[nCBins][nPtBins][nTrkPtBins];
	TH1D* background_proj[nCBins][nPtBins][nTrkPtBins];
	TH1D* background_proj_cc26[nCBins][nPtBins][nTrkPtBins];


	TH1D *yield_inc_SBL[nCBins][nPtBins][nTrkPtBins];              // The Y projection of the yield (phi) at the left range
	TH1D *yield_inc_SBR[nCBins][nPtBins][nTrkPtBins];		// The Y projection of the yield (phi) at the right range 
	TH1D *yield_inc_Signal[nCBins][nPtBins][nTrkPtBins];		// The Y projection of the yield (phi) at the center (the signal) 
	TH1D *yield_inc_SBL_ForRatio[nCBins][nPtBins][nTrkPtBins];		// The Y projection of the yield (phi) at the center (the signal) 

	TH1D *yield_corrected_inc_SBL_corrected[nCBins][nPtBins][nTrkPtBins];

	 //TH1D *check_old_phi_rebin[12][5][4][2];
	 //TH1D *check_old_phi_rebin[5][4][2];
	 TH1D *check_old_phi[nCBins][nPtBins][nTrkPtBins];
	 TH1D *check_old_phi_rebin[nCBins][nPtBins][nTrkPtBins];
	 TH1D *check_new_phi[nCBins][nPtBins][nTrkPtBins];

	TF1 *bg_fit[nCBins][nPtBins][nTrkPtBins]; 
	     //The Y projection of signal/ME, Left band



	char *histname = new char[10];
	int xx = 0;
	double temp1, err1;
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
	TH1D *Signal_cc26[nCBins][nPtBins][nTrkPtBins];
	TH1D *SignalMinusRL[nCBins][nPtBins][nTrkPtBins];
	TH1D *SignalMinusRL_cc10[nCBins][nPtBins][nTrkPtBins];

	TH1D *Aj[nCBins];	

	TString desc = "hists";

	for (int ibin=0;ibin<nCBins;ibin++){
		Aj[ibin] = (TH1D*)fin->Get((TString)(desc + "_Aj"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]))->Clone((TString) ("Aj_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]));    
		for (int ibin2=0;ibin2<nPtBins;ibin2++){ 
			for (int ibin3=0;ibin3<nTrkPtBins-1;ibin3++){
				cout<<ibin<<" "<<ibin2<<" "<<ibin3<<endl;
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
				
			//	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

				hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
				
				//hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

				Correlation[ibin][ibin2][ibin3] = (TH2D*)  fin->Get((TString) (desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	Scaled_Normalized_Yield_Inv[ibin][ibin2][ibin3] = (TH2D*)fin_LevelOne->Get((TString)("Scaled_Normalized_Yield_Inv_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));


			} // end of ibin3
		} //end of ibin2

		cout << "after 2 loops" << endl;
		cout << " Aj based Number of Dijets = " << Aj[ibin]->GetEntries() << " For Centrality bin = " << ibin+1 << endl;
		int AjBinMax = Aj[ibin]->GetMaximumBin();
		float AjMax =  Aj[ibin]->GetXaxis()->GetBinUpEdge(AjBinMax);
		cout << " Aj for this sampel = " << AjMax << endl;
		
		//TFile *f = TFile::Open("/home/mzakaria3/Documents/Research/Data2011_trkPtCut1_part1.root");
		//TFile *f = TFile::Open("Data2011_trkPtCut1_All_AJ11.root");
		//cout << "took the file" << endl;

		//hists_hJetTrackMECent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Draw("LEGO2");
		//cout << "Plot of the Mixed Events" << endl;

		//hists_hJetTrackMECent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionX("ProEta")->Draw();
		//cout << "Eta Projection" <<  endl;

		//hists_hJetTrackMECent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("me_projphi")->Draw();
		//cout << "Phi Projection" << endl;

		//cout << "Number of Eta bins:" << ProEta->GetNbinsX() << endl;
		//float EtaBins = ProEta->GetNbinsX() ;
		//cout << "Number of Phi bins:" << me_projphi->GetNbinsX() << endl;

		//Normalizing

		for(int ibin2=0; ibin2<nPtBins;ibin2++){
			for (int ibin3=0; ibin3< 4; ibin3++){

		 lbin = Correlation[ibin][ibin2][ibin3]->GetYaxis()->FindBin(-1.5);
		rbin = Correlation[ibin][ibin2][ibin3]->GetYaxis()->FindBin(1.5);
		//		printf("the line number is: %d \n \n", __LINE__);

				cout << "hJetTrackSignalBackground " << ibin << " " << ibin2 << " " << ibin3 << endl;
			//	hJetTrackSignalBackground[ibin][ibin2][ibin3]->Draw("LEGO2");

				cout << "Fitting the central part of eta" << endl;

			//	fout_inc->cd();
				//cout << "Canvas number = " << 4*(ibin3+1)-ibin << endl;
				//me_proj_canvas->cd(4*(ibin3+1)-ibin);
				
				yield_inc[ibin][ibin2][ibin3] =  (TH2D*) hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Clone((TString)("Yield_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

			   yield_inc_copy[ibin][ibin2][ibin3] =  (TH2D*) hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Clone((TString)("Yield_copy"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

				
				//ProEta->Fit("pol0","","",-0.2,0.2);
				//float s = pol0->GetParameter(0);  // need to make this an array
				me_projeta[ibin][ibin2][ibin3] =  hJetTrackMELeading[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_eta%d%d%d",ibin,ibin2,ibin3),1,100); 
				me_projphi[ibin][ibin2][ibin3] =  hJetTrackMELeading[ibin][ibin2][ibin3]->ProjectionY(Form("me_proj_phi%d%d%d",ibin,ibin2,ibin3),1,100); //what is the 1,100 for?
				//me_projeta[ibin][ibin2][ibin3]->Scale (1./100);
				EtaProj[ibin][ibin2][ibin3] = (TH1D*)Correlation[ibin][ibin2][ibin3]->ProjectionX((TString)("Eta_Proj_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),lbin,rbin);
	  
				EtaProjDummy[ibin][ibin2][ibin3] = (TH1D*) EtaProj[ibin][ibin2][ibin3]->Clone((TString)("Eta_Proj_Dummy_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

			

				cout << "Trying to get Number of eta bins = " << endl;
				cout << me_projeta[ibin][ibin2][ibin3]->GetNbinsX() << endl;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			dd1->cd((4*(ibin3+1)-ibin));
			 
			Scaled_Normalized_Yield_Inv[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
			Scaled_Normalized_Yield_Inv[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
			Scaled_Normalized_Yield_Inv[ibin][ibin2][ibin3]->GetXaxis() -> CenterTitle();
			Scaled_Normalized_Yield_Inv[ibin][ibin2][ibin3]->GetYaxis() -> CenterTitle();
			Scaled_Normalized_Yield_Inv[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(1.6);
			Scaled_Normalized_Yield_Inv[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(1.6);
			Scaled_Normalized_Yield_Inv[ibin][ibin2][ibin3]->SetAxisRange(-3., 3.,"X");
			Scaled_Normalized_Yield_Inv[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-3.,3.);	
			Scaled_Normalized_Yield_Inv[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
			Scaled_Normalized_Yield_Inv[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
			Scaled_Normalized_Yield_Inv[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
			Scaled_Normalized_Yield_Inv[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.06);

			Scaled_Normalized_Yield_Inv[ibin][ibin2][ibin3]->Draw("LEGO2");


			//PhiBinWidth =  me_projphi[ibin][ibin2][ibin3]->GetBinWidth(1) ;
			EtaBinWidth =  Scaled_Normalized_Yield_Inv[ibin][ibin2][ibin3]->GetXaxis()->GetBinWidth(1) ;
			PhiBinWidth =  Scaled_Normalized_Yield_Inv[ibin][ibin2][ibin3]->GetYaxis()->GetBinWidth(1) ;


			Signal[ibin][ibin2][ibin3] =    (TH1D*)Scaled_Normalized_Yield_Inv[ibin][ibin2][ibin3]->ProjectionY(Form("Phi_Signal_dd1_%d%d%d",ibin,ibin2,ibin3),36, 65);	
			std::cout << "Signal1: Number of bins (size) = " << Signal[ibin][ibin2][ibin3]->GetSize() << std::endl; 
			
			Signal[ibin][ibin2][ibin3]->SetMarkerStyle(20);
			Signal[ibin][ibin2][ibin3]->SetMarkerColor(kRed);
			//	Signal[ibin][ibin2][ibin3]->Rebin(4);
			Signal[ibin][ibin2][ibin3]->Rebin(2);
			std::cout << "Signal2: Number of bins (size) = " << Signal[ibin][ibin2][ibin3]->GetSize() << std::endl; 				
			Signal[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#phi");
			Signal[ibin][ibin2][ibin3]->Scale(EtaBinWidth/2.);
			//	Signal[ibin][ibin2][ibin3]->Scale(EtaBinWidth/4.);
			//	Signal[ibin][ibin2][ibin3]->Scale(1. / 2.);
			//	Signal[ibin][ibin2][ibin3]->Scale(1. / 29.);
				//Signal[ibin][ibin2][ibin3]->Scale(EtaBinWidth/4.);
				//Signal[ibin][ibin2][ibin3]->Scale(1./(4. * 29));
			Signal[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleOffset(0.8);
			Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleOffset(0.8);
			Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
			Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
			Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleOffset(0.8);
			Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
			Signal[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
			Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
			Signal[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.05);
			Signal[ibin][ibin2][ibin3]->SetMaximum(TMath::Max(Signal[ibin][ibin2][ibin3]->GetMaximum()+5.0, Signal[ibin][ibin2][ibin3]->GetMaximum())*1.05);

			Signal[ibin][ibin2][ibin3]->Draw();
				
			std::cout << "Signal: Number of bins (size) = " << Signal[ibin][ibin2][ibin3]->GetSize() << std::endl; 

			//PhiProj_L[ibin][ibin2][ibin3] = (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionY(Form("ResultPhi_L_cc9_%d%d%d",ibin,ibin2,ibin3),10, 35);
			PhiProj_L[ibin][ibin2][ibin3] = (TH1D*)Scaled_Normalized_Yield_Inv[ibin][ibin2][ibin3]->ProjectionY(Form("ResultPhi_L_dd1_%d%d%d",ibin,ibin2,ibin3),21, 35);
			//PhiProj_L[ibin][ibin2][ibin3]->Draw("same");
			//PhiProj_L[ibin][ibin2][ibin3]->Scale(1./14.);	
			std::cout << "PhiProj_L1: Number of bins (size) = " << PhiProj_L[ibin][ibin2][ibin3]->GetSize() << std::endl; 
			
			PhiProj_R[ibin][ibin2][ibin3] = (TH1D*)Scaled_Normalized_Yield_Inv[ibin][ibin2][ibin3]->ProjectionY(Form("ResultPhi_R_dd1_%d%d%d",ibin,ibin2,ibin3),66, 80);	
			//PhiProj_R[ibin][ibin2][ibin3] = (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionY(Form("ResultPhi_R_cc9_%d%d%d",ibin,ibin2,ibin3),66, 90);
			//PhiProj_R[ibin][ibin2][ibin3]->Scale(1./14.);				
			//PhiProj_R[ibin][ibin2][ibin3]->Draw("same");
			std::cout << "PhiProj_R1: Number of bins (size) = " << PhiProj_R[ibin][ibin2][ibin3]->GetSize() << std::endl; 


			PhiProj_R[ibin][ibin2][ibin3]->Add(PhiProj_L[ibin][ibin2][ibin3], 1);
			PhiProj_R[ibin][ibin2][ibin3]->SetMarkerStyle(20);
			PhiProj_R[ibin][ibin2][ibin3]->SetMarkerColor(kBlue);
			std::cout << "PhiProj_R2: Number of bins (size) = " << PhiProj_R[ibin][ibin2][ibin3]->GetSize() << std::endl; 
				
			PhiProj_R[ibin][ibin2][ibin3]->Rebin(2);
			std::cout << "PhiProj_R3: Number of bins (size) = " << PhiProj_R[ibin][ibin2][ibin3]->GetSize() << std::endl; 

			//PhiProj_R[ibin][ibin2][ibin3]->Rebin(2);
			PhiProj_R[ibin][ibin2][ibin3]->Scale(EtaBinWidth/2.);	
			//PhiProj_R[ibin][ibin2][ibin3]->Scale(1./2.);	
			//	PhiProj_R[ibin][ibin2][ibin3]->Scale(EtaBinWidth/4.);	
			//PhiProj_R[ibin][ibin2][ibin3]->Scale(1./ (4.));	
			//	PhiProj_R[ibin][ibin2][ibin3]->Scale(1./ (29.));	
		
			//	if(ibin == 0 && ibin3 == 0)
			//	PhiProj_R[ibin][ibin2][ibin3]->SetMaximum(120);

			PhiProj_R[ibin][ibin2][ibin3]->Draw("same");
			std::cout << "PhiProj_R: Number of bins (size) = " << PhiProj_R[ibin][ibin2][ibin3]->GetSize() << std::endl; 
			//yield_inc_proj[ibin][ibin2][ibin3]->Draw();	
			std::cout<<"BBB Bin Content:"<<Signal[ibin][ibin2][ibin3]->GetBinContent(1)<<std::endl;	


			Signal[ibin][ibin2][ibin3]->Draw("same");

	{
					TLatex *centtex = new TLatex(0.15,0.85,CBin_labels[ibin]);
					centtex->SetNDC();
					centtex->Draw();
					TLatex *pttex = new TLatex(0.15,0.78,TrkPtBin_labels[ibin3]);
					pttex->SetNDC();
					pttex->Draw();
					TLatex *Ajtex = new TLatex(0.75,0.85,AjBin_labels[0]);
					Ajtex->SetNDC();
					Ajtex->Draw();
					TLatex *Jettex = new TLatex(0.70, 0.80, JetLabel_labels[0]);
					Jettex->SetNDC();
					Jettex->Draw();
					TLatex *COMtex = new TLatex(0.15, 0.93, COM_label[0]);
					COMtex->SetNDC();
				//	if(ibin == 3)
					COMtex->Draw();
					TLatex *CMStex = new TLatex(0.15, 0.70, CMS_label[0]);
					CMStex->SetNDC();
				//	if(ibin == 2)
					CMStex->Draw();
					TLatex *Pbtex = new TLatex(0.65, 0.93, Pb_label[0]);
					Pbtex->SetNDC();
					Pbtex->Draw();
				}

					
						
			dd1->SaveAs("dd1_PhiProjectionofSignalandBands_Leading_AjUnder22.png");	
			
			dd2->cd((4*(ibin3+1)-ibin));
		
			Signal[ibin][ibin2][ibin3]->Draw();
			
			Signal[ibin][ibin2][ibin3]->Add(	PhiProj_R[ibin][ibin2][ibin3], -1);
			Signal[ibin][ibin2][ibin3]->Draw();
			Signal[ibin][ibin2][ibin3]->SetAxisRange(-1.5, 1.5,"X");
	//		Signal[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-1.5, 1.5);
			Signal[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("1/N_{jet} dN/(d#Delta#phi dp_{T})");
			Signal[ibin][ibin2][ibin3]->GetYaxis()->SetTitleSize(0.06);
			Signal[ibin][ibin2][ibin3]->GetYaxis()->CenterTitle();
			Signal[ibin][ibin2][ibin3]->GetYaxis()->SetTitleOffset(0.6);
	
			
			{
					TLatex *centtex = new TLatex(0.15,0.85,CBin_labels[ibin]);
					centtex->SetNDC();
					centtex->Draw();
					TLatex *pttex = new TLatex(0.15,0.78,TrkPtBin_labels[ibin3]);
					pttex->SetNDC();
					pttex->Draw();
					TLatex *Ajtex = new TLatex(0.75,0.85,AjBin_labels[0]);
					Ajtex->SetNDC();
					Ajtex->Draw();
					TLatex *Jettex = new TLatex(0.70, 0.80, JetLabel_labels[0]);
					Jettex->SetNDC();
					Jettex->Draw();
					TLatex *COMtex = new TLatex(0.15, 0.93, COM_label[0]);
					COMtex->SetNDC();
				//	if(ibin == 3)
					COMtex->Draw();
					TLatex *CMStex = new TLatex(0.15, 0.70, CMS_label[0]);
					CMStex->SetNDC();
				//	if(ibin == 2)
					CMStex->Draw();
					TLatex *Pbtex = new TLatex(0.65, 0.93, Pb_label[0]);
					Pbtex->SetNDC();
					Pbtex->Draw();

			}

			dd2->SaveAs("dd2_resultOnPhi_Leading_AjUnder.png");
			

			dd3->cd((4*(ibin3+1)-ibin));
			Signal[ibin][ibin2][ibin3]->Draw();
		        double etalim = 1.5;
			double dx_phi = 0.;
			llimiteta = Scaled_Normalized_Yield_Inv[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-etalim + .001);
			rlimiteta = Scaled_Normalized_Yield_Inv[ibin][ibin2][ibin3]->GetXaxis()->FindBin( etalim - .001);

			cout << "llimiteta = " << llimiteta << " rlimiteta = " << rlimiteta << endl; 

			check_old_phi[ibin][ibin2][ibin3] =  (TH1D*)Scaled_Normalized_Yield_Inv[ibin][ibin2][ibin3]->ProjectionY(Form("Phi_Signal_dd3%d%d%d",ibin,ibin2,ibin3),llimiteta,rlimiteta);	
		//	check_old_phi[ibin][ibin2][ibin3]->Draw();
//	return 3;
			//check_old_phi[ibin][ibin2][ibin3]->Scale(1./30);
			//check_old_phi[ibin][ibin2][ibin3]->Scale(EtaBinWidth);
			dx_phi = check_old_phi[ibin][ibin2][ibin3]->GetBinWidth(1);
			//check_old_phi[ibin][ibin2][ibin3]->Scale(dx_phi);
			
			check_old_phi[ibin][ibin2][ibin3]->Draw();
//return 3;
			Signal[ibin][ibin2][ibin3]->Scale(1./PhiBinWidth);
			check_old_phi_rebin[ibin][ibin2][ibin3] = (TH1D*)Rebin_dPhi(check_old_phi[ibin][ibin2][ibin3]);
		//	check_old_phi_rebin[ibin][ibin2][ibin3] = (TH1D*)Rebin_dPhi(Signal[ibin][ibin2][ibin3]);
		//	check_old_phi_rebin[ibin][ibin2][ibin3] = (TH1D*)Rebin_dEta(check_old_phi[ibin][ibin2][ibin3]);
		//	check_old_phi_rebin[ibin][ibin2][ibin3]->Rebin(2);
		//	check_old_phi_rebin[ibin][ibin2][ibin3]->Scale(1./2);
			dx_phi = check_old_phi[ibin][ibin2][ibin3]->GetBinWidth(1);
			cout << "dx_phi = " << dx_phi << endl;	
		//	check_old_phi[ibin][ibin2][ibin3]->Scale(1. / dx_phi);							
			check_old_phi[ibin][ibin2][ibin3]->Draw();
//return 3;
			check_new_phi[ibin][ibin2][ibin3] =  (TH1D*)Scaled_Normalized_Yield_Inv[ibin][ibin2][ibin3]->ProjectionY(Form("Phi_Signal2_dd3%d%d%d",ibin,ibin2,ibin3),llimiteta,rlimiteta);
			//check_new_phi[ibin][ibin2][ibin3]->Scale(1/dx_phi);
			//check_new_phi[ibin][ibin2][ibin3]->Scale(dx_phi);
	
			check_old_phi_rebin[ibin][ibin2][ibin3]->SetMarkerStyle(20);
			check_old_phi_rebin[ibin][ibin2][ibin3]->SetMarkerColor(kRed);
			check_old_phi_rebin[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#phi");
			check_old_phi_rebin[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#frac{1}{N_{jet}}#frac{d^{2}N}{ d#Delta#phi dp_{T}}");
			check_old_phi_rebin[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(0.8);
			check_old_phi_rebin[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(0.8);
			check_old_phi_rebin[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
			check_old_phi_rebin[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
			check_old_phi_rebin[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(0.8);
			check_old_phi_rebin[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
			check_old_phi_rebin[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);   			
			check_old_phi_rebin[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
			check_old_phi_rebin[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.05);   	//	result[ibin][ibin2][ibin3]->Draw("LEGO2");
			check_old_phi_rebin[ibin][ibin2][ibin3]->GetXaxis()-> SetRangeUser(-1.5,1.5);	
			check_old_phi_rebin[ibin][ibin2][ibin3]->GetYaxis()-> SetRangeUser(-0.5, 7);
//			check_old_phi_rebin[ibin][ibin2][ibin3]->SetMaximum(TMath::Max(check_old_phi_rebin[ibin][ibin2][ibin3]->GetMaximum()+5.0, check_old_phi_rebin[ibin][ibin2][ibin3]->GetMaximum())*1.05);
			
			//check_old_phi_rebin[ibin][ibin2][ibin3]->Rebin(3);
			//check_old_phi_rebin[ibin][ibin2][ibin3]->Scale(1/3.);
			//check_old_phi_rebin[ibin][ibin2][ibin3]->Draw();
			//	check_old_phi_rebin[ibin][ibin2][ibin3]->Draw("Y+");
			//	check_old_phi_rebin[ibin][ibin2][ibin3]->Draw("Y-", "SAME");
		
			check_old_phi_rebin[ibin][ibin2][ibin3]->Draw();
			//check_new_phi[ibin][ibin2][ibin3]->Draw();
			
		
				
				if(ibin==3){
					//TLatex *centtex2 = new TLatex(0.15,0.9,CBin_labels[ibin]);
					TLatex *centtex2 = new TLatex(0.20,0.90,CBin_labels[ibin]);
					centtex2->SetNDC();
					centtex2->Draw();
					//TLatex *pttex2 = new TLatex(0.15,0.83,TrkPtBin_labels[ibin3]);
					TLatex *pttex2 = new TLatex(0.20,0.83,TrkPtBin_labels[ibin3]);
					pttex2->SetNDC();
					pttex2->Draw();
					TLatex *Ajtex2 = new TLatex(0.85,0.90,AjBin_labels[0]);
					Ajtex2->SetNDC();
					Ajtex2->Draw();
					TLatex *Jettex2 = new TLatex(0.80, 0.85, JetLabel_labels[0]);
					Jettex2->SetNDC();
					Jettex2->Draw();


				}
	
				//	if(ibin<3){
				//	TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
				//	centtex->SetNDC();
				//	centtex->Draw();
				//	TLatex *pttex = new TLatex(0.05,0.83,TrkPtBin_labels[ibin3]);
				//	pttex->SetNDC();
				//	pttex->Draw();
				//	TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
				//	Ajtex->SetNDC();
				//	Ajtex2->Draw();
				//	TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
				//	Jettex->SetNDC();
				//	Jettex2->Draw();
				//}
	
		
			dd3->SaveAs("dd3_resultOnPhi_Leading_AjUnder22.png");
	
					cc9->cd((4*(ibin3+1)-ibin));
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


			//hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
			hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));

			//hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(100./me00[ibin][ibin2][ibin3]);
				cout << "Trying to get Number of eta bins = " << endl;
				cout << me_projeta[ibin][ibin2][ibin3]->GetNbinsX() << endl;

				PhiBinWidth =  me_projphi[ibin][ibin2][ibin3]->GetBinWidth(1) ;
				EtaBinWidth =  hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->GetBinWidth(1) ;
				PhiBinWidth =  hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->GetBinWidth(1) ;
				cout << "Eta Bin Width = " << 	 hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->GetBinWidth(1) << endl;
				cout << "####################ibin = " << ibin << endl;

		me_projeta[ibin][ibin2][ibin3] =  hJetTrackMELeading[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_eta_cc9_%d%d%d",ibin,ibin2,ibin3),1,100);	
				///me_projeta[ibin][ibin2][ibin3]->Scale(1. / 100.);	
				me_projeta[ibin][ibin2][ibin3]->Scale(1. / 99.);
				me_projeta[ibin][ibin2][ibin3]->Draw();	
				me_projeta[ibin][ibin2][ibin3]->Fit("pol0","","",-.2,.2);
				me_projeta[ibin][ibin2][ibin3]->Draw();	
				me00[ibin][ibin2][ibin3] = me_projeta[ibin][ibin2][ibin3]->GetFunction("pol0")->GetParameter(0);	

	
			hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1./me00[ibin][ibin2][ibin3]);
			//hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() * PhiBinWidth*EtaBinWidth));					
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() * PhiBinWidth*EtaBinWidth));	
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Divide(hJetTrackMELeading[ibin][ibin2][ibin3]);

				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> CenterTitle();
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> CenterTitle();
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(1.6);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(1.6);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->SetAxisRange(-3., 3.,"X");

				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-3.,3.);	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.06);

			//	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

				
				
			//	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Divide(hJetTrackMELeading[ibin][ibin2][ibin3]);

				 //hists_hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Rebin2D(5,2);
				

				//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Draw("LEGO2");

				//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionX("Check%d%d%d",1, 50)->Draw();
				//yield_inc_proj[ibin][ibin2][ibin3] = (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionX((Form("ResultEta%d%d%d",ibin,ibin2,ibin3)),1, 50);

				//TH1D *Signal = (TH1D*)hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("Signal",36, 65);
				//Signal[ibin][ibin2][ibin3] =    (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionY(Form("Phi_Signal_cc9_%d%d%d",ibin,ibin2,ibin3),36, 66);	
				Signal[ibin][ibin2][ibin3] =    (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionY(Form("Phi_Signal_cc9_%d%d%d",ibin,ibin2,ibin3),36, 65);	
				std::cout << "Signal1: Number of bins (size) = " << Signal[ibin][ibin2][ibin3]->GetSize() << std::endl; 
			
				Signal[ibin][ibin2][ibin3]->SetMarkerStyle(20);
				Signal[ibin][ibin2][ibin3]->SetMarkerColor(kRed);
			//	Signal[ibin][ibin2][ibin3]->Rebin(4);
				Signal[ibin][ibin2][ibin3]->Rebin(2);
				std::cout << "Signal2: Number of bins (size) = " << Signal[ibin][ibin2][ibin3]->GetSize() << std::endl; 				
				Signal[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#phi");
				Signal[ibin][ibin2][ibin3]->Scale(EtaBinWidth/2.);
			//	Signal[ibin][ibin2][ibin3]->Scale(EtaBinWidth/4.);
			//	Signal[ibin][ibin2][ibin3]->Scale(1. / 2.);
			//	Signal[ibin][ibin2][ibin3]->Scale(1. / 29.);
				//Signal[ibin][ibin2][ibin3]->Scale(EtaBinWidth/4.);
				//Signal[ibin][ibin2][ibin3]->Scale(1./(4. * 29));
				Signal[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleOffset(0.8);
				Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleOffset(0.8);
				Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleOffset(0.8);
				Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				Signal[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
				Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				Signal[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.05);

			//	if(ibin == 0 && ibin3 == 0)   // top right
				Signal[ibin][ibin2][ibin3]->SetMaximum(TMath::Max(Signal[ibin][ibin2][ibin3]->GetMaximum()+5.0, Signal[ibin][ibin2][ibin3]->GetMaximum())*1.05);
		


				Signal[ibin][ibin2][ibin3]->Draw();
				std::cout << "Signal: Number of bins (size) = " << Signal[ibin][ibin2][ibin3]->GetSize() << std::endl; 

				//PhiProj_L[ibin][ibin2][ibin3] = (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionY(Form("ResultPhi_L_cc9_%d%d%d",ibin,ibin2,ibin3),10, 35);
				PhiProj_L[ibin][ibin2][ibin3] = (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionY(Form("ResultPhi_L_cc9_%d%d%d",ibin,ibin2,ibin3),21, 35);
				//PhiProj_L[ibin][ibin2][ibin3]->Draw("same");
				//PhiProj_L[ibin][ibin2][ibin3]->Scale(1./14.);	
				std::cout << "PhiProj_L1: Number of bins (size) = " << PhiProj_L[ibin][ibin2][ibin3]->GetSize() << std::endl; 
			

				PhiProj_R[ibin][ibin2][ibin3] = (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionY(Form("ResultPhi_R_cc9_%d%d%d",ibin,ibin2,ibin3),66, 80);	
				//PhiProj_R[ibin][ibin2][ibin3] = (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionY(Form("ResultPhi_R_cc9_%d%d%d",ibin,ibin2,ibin3),66, 90);
				//PhiProj_R[ibin][ibin2][ibin3]->Scale(1./14.);				
				//PhiProj_R[ibin][ibin2][ibin3]->Draw("same");
				std::cout << "PhiProj_R1: Number of bins (size) = " << PhiProj_R[ibin][ibin2][ibin3]->GetSize() << std::endl; 


				PhiProj_R[ibin][ibin2][ibin3]->Add(PhiProj_L[ibin][ibin2][ibin3], 1);
				PhiProj_R[ibin][ibin2][ibin3]->SetMarkerStyle(20);
				PhiProj_R[ibin][ibin2][ibin3]->SetMarkerColor(kBlue);
				std::cout << "PhiProj_R2: Number of bins (size) = " << PhiProj_R[ibin][ibin2][ibin3]->GetSize() << std::endl; 
				
				PhiProj_R[ibin][ibin2][ibin3]->Rebin(2);
				std::cout << "PhiProj_R3: Number of bins (size) = " << PhiProj_R[ibin][ibin2][ibin3]->GetSize() << std::endl; 

				//PhiProj_R[ibin][ibin2][ibin3]->Rebin(2);
				PhiProj_R[ibin][ibin2][ibin3]->Scale(EtaBinWidth/2.);	
				//PhiProj_R[ibin][ibin2][ibin3]->Scale(1./2.);	
		//		PhiProj_R[ibin][ibin2][ibin3]->Scale(EtaBinWidth/4.);	
				//PhiProj_R[ibin][ibin2][ibin3]->Scale(1./ (4.));	
			//	PhiProj_R[ibin][ibin2][ibin3]->Scale(1./ (29.));	
		
				//	if(ibin == 0 && ibin3 == 0)
				//	PhiProj_R[ibin][ibin2][ibin3]->SetMaximum(120);

				PhiProj_R[ibin][ibin2][ibin3]->Draw("same");
	std::cout << "PhiProj_R: Number of bins (size) = " << PhiProj_R[ibin][ibin2][ibin3]->GetSize() << std::endl; 
				//yield_inc_proj[ibin][ibin2][ibin3]->Draw();	
	std::cout<<"BBB Bin Content:"<<Signal[ibin][ibin2][ibin3]->GetBinContent(1)<<std::endl;	


//			Signal[ibin][ibin2][ibin3]->Add(PhiProj_R[ibin][ibin2][ibin3], -1);
//cc9->Modified();
//cc9->Update();
//		std::cout<<"AAA Bin Content:"<<Signal[ibin][ibin2][ibin3]->GetBinContent(1)<<std::endl;	

	Signal[ibin][ibin2][ibin3]->Draw("same");

					{
					TLatex *centtex = new TLatex(0.15,0.85,CBin_labels[ibin]);
					centtex->SetNDC();
					centtex->Draw();
					TLatex *pttex = new TLatex(0.15,0.78,TrkPtBin_labels[ibin3]);
					pttex->SetNDC();
					pttex->Draw();
					TLatex *Ajtex_dd3 = new TLatex(0.75,0.85,AjBin_labels[0]);
					Ajtex_dd3->SetNDC();
					Ajtex_dd3->Draw();
					TLatex *Jettex_dd3 = new TLatex(0.70, 0.80, JetLabel_labels[0]);
					Jettex_dd3->SetNDC();
					Jettex_dd3->Draw();
					TLatex *COMtex = new TLatex(0.15, 0.93, COM_label[0]);
					COMtex->SetNDC();
				//	if(ibin == 3)
					COMtex->Draw();
					TLatex *CMStex = new TLatex(0.15, 0.70, CMS_label[0]);
					CMStex->SetNDC();
				//	if(ibin == 2)
					CMStex->Draw();
					TLatex *Pbtex = new TLatex(0.65, 0.93, Pb_label[0]);
					Pbtex->SetNDC();
					Pbtex->Draw();
				}
		//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Draw("LEGO2");
		//	cc9->SaveAs("cc9_PhiProjectionofSignalandBands_Leading_Under22.png");			
//return 9;
			cc10->cd(4*(ibin3+1)-ibin);
			SignalMinusRL_cc10[ibin][ibin2][ibin3] = (TH1D*)Signal[ibin][ibin2][ibin3]->Clone((TString) ("SignalMinusRL_cc10_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
		//	SignalMinusRL_cc10[ibin][ibin2][ibin3]->SetMarkerColor(kRed);
		//		SignalMinusRL_cc10[ibin][ibin2][ibin3]->Draw();

			PhiProj_R[ibin][ibin2][ibin3] = (TH1D*)PhiProj_R[ibin][ibin2][ibin3]->Clone((TString) ("Phiroj_R_cc10_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

		//	PhiProj_R[ibin][ibin2][ibin3]->Draw("same");
		//	PhiProj_R[ibin][ibin2][ibin3]->Draw();
		//	SignalMinusRL_cc10[ibin][ibin2][ibin3]->Add(	PhiProj_R[ibin][ibin2][ibin3], -1);
		//	SignalMinusRL_cc10[ibin][ibin2][ibin3]->Draw();
			
	//	std::cout<<"BBB Bin Content:"<<Signal[ibin][ibin2][ibin3]->GetBinContent(1)<<std::endl;	
		SignalMinusRL_cc10[ibin][ibin2][ibin3]->Add(	PhiProj_R[ibin][ibin2][ibin3], -1);
	//	std::cout<<"AAA Bin Content:"<<Signal[ibin][ibin2][ibin3]->GetBinContent(1)<<std::endl;	
		SignalMinusRL_cc10[ibin][ibin2][ibin3]->Draw();
		SignalMinusRL_cc10[ibin][ibin2][ibin3]->SetAxisRange(-1.5, 1.5,"X");
	//		SignalMinusRL_cc10[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-1.5, 1.5);
		SignalMinusRL_cc10[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("1/N_{jet} dN/(d#Delta#phi dp_{T})");
		SignalMinusRL_cc10[ibin][ibin2][ibin3]->GetYaxis()->SetTitleSize(0.06);
		SignalMinusRL_cc10[ibin][ibin2][ibin3]->GetYaxis()->CenterTitle();
		SignalMinusRL_cc10[ibin][ibin2][ibin3]->GetYaxis()->SetTitleOffset(0.6);
			
		//	SignalMinusRL[ibin][ibin2][ibin3]->Scale(1/(3.0 *  Aj[ibin]->GetEntries() * PhiBinWidth));

//			if(ibin3==0)	SignalMinusRL[ibin][ibin2][ibin3]->SetMinimum(-2.);
//			if(ibin3==0)	SignalMinusRL[ibin][ibin2][ibin3]->SetMaximum(8.);


		//	SignalMinusRL_cc10[ibin][ibin2][ibin3]->SetMarkerStyle(20);
	//		SignalMinusRL_cc10[ibin][ibin2][ibin3]->SetMarkerColor(kBlue);
//			SignalMinusRL[ibin][ibin2][ibin3]->SetMaximum(TMath::Max(SignalMinusRL[ibin][ibin2][ibin3]->GetMaximum()+5.0, SignalMinusRL[ibin][ibin2][ibin3]->GetMaximum()*1.05));
			
			//SignalMinusRL_cc10[ibin][ibin2][ibin3]->Draw();
				
			{
					TLatex *centtex = new TLatex(0.15,0.85,CBin_labels[ibin]);
					centtex->SetNDC();
					centtex->Draw();
					TLatex *pttex = new TLatex(0.15,0.78,TrkPtBin_labels[ibin3]);
					pttex->SetNDC();
					pttex->Draw();
					TLatex *Ajtex_cc10 = new TLatex(0.75,0.85,AjBin_labels[0]);
					Ajtex_cc10->SetNDC();
					Ajtex_cc10->Draw();
					TLatex *Jettex_cc10 = new TLatex(0.70, 0.80, JetLabel_labels[0]);
					Jettex_cc10->SetNDC();
					Jettex_cc10->Draw();
					TLatex *COMtex = new TLatex(0.15, 0.93, COM_label[0]);
					COMtex->SetNDC();
				//	if(ibin == 3)
					COMtex->Draw();
					TLatex *CMStex = new TLatex(0.15, 0.70, CMS_label[0]);
					CMStex->SetNDC();
				//	if(ibin == 2)
					CMStex->Draw();
					TLatex *Pbtex = new TLatex(0.65, 0.93, Pb_label[0]);
					Pbtex->SetNDC();
					Pbtex->Draw();

				}
			cc10->SaveAs("cc10_resultOnPhi_Leading_AjUnder22.png");


				cc11->cd(4*(ibin3+1)-ibin);
			SignalMinusRL[ibin][ibin2][ibin3] = (TH1D*)Signal[ibin][ibin2][ibin3]->Clone((TString)("signal_MinusRL"+CBin_strs[ibin]+"_"+TrkPtBin_strs[ibin3]));
			SignalMinusRL[ibin][ibin2][ibin3]->Add(	PhiProj_R[ibin][ibin2][ibin3], -1);
			SignalMinusRL[ibin][ibin2][ibin3]->SetAxisRange(-1.5, 1.5,"X");
			SignalMinusRL[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("1/N_{jet} dN/(d#Delta#phi dp_{T})");
		//	SignalMinusRL[ibin][ibin2][ibin3]->GetYaxis()->SetTitleSize(0.06);
			SignalMinusRL[ibin][ibin2][ibin3]->GetYaxis()->CenterTitle();
			SignalMinusRL[ibin][ibin2][ibin3]->GetYaxis()->SetTitleOffset(0.6);
			
		//	SignalMinusRL[ibin][ibin2][ibin3]->Scale(1/(3.0 *  Aj[ibin]->GetEntries() * PhiBinWidth));

			if(ibin3==0)	SignalMinusRL[ibin][ibin2][ibin3]->SetMinimum(-2.);
			if(ibin3==0)	SignalMinusRL[ibin][ibin2][ibin3]->SetMaximum(8.);


			SignalMinusRL[ibin][ibin2][ibin3]->SetMarkerStyle(20);
			SignalMinusRL[ibin][ibin2][ibin3]->SetMarkerColor(kBlue);
			SignalMinusRL[ibin][ibin2][ibin3]->Draw();
				
			{
					TLatex *centtex = new TLatex(0.15,0.85,CBin_labels[ibin]);
					if(ibin != 3)
					{
					TLatex *centtex = new TLatex(0.05,0.85,CBin_labels[ibin]);
					}
					centtex->SetNDC();
					centtex->Draw();
					TLatex *pttex = new TLatex(0.15,0.78,TrkPtBin_labels[ibin3]);
					if(ibin !=3)
					{
					TLatex *pttex = new TLatex(0.05,0.78,TrkPtBin_labels[ibin3]);
					}
					pttex->SetNDC();
					pttex->Draw();
					TLatex *Ajtex = new TLatex(0.85,0.85,AjBin_labels[0]);
					Ajtex->SetNDC();
					Ajtex->Draw();
					TLatex *Jettex = new TLatex(0.80, 0.80, JetLabel_labels[0]);
					Jettex->SetNDC();
					Jettex->Draw();
					TLatex *COMtex = new TLatex(0.15, 0.93, COM_label[0]);
					COMtex->SetNDC();
				if(ibin == 3)
					COMtex->Draw();
					TLatex *CMStex = new TLatex(0.05, 0.92, CMS_label[0]);
					CMStex->SetNDC();
				if(ibin == 2)
					CMStex->Draw();
					TLatex *Pbtex = new TLatex(0.75, 0.93, Pb_label[0]);
					Pbtex->SetNDC();
				if(ibin == 3)
					Pbtex->Draw();
					TLatex *etatex = new TLatex(0.85, 0.91, eta_label[0]);
					etatex->SetNDC();
				if(ibin == 2)
					etatex->Draw();
					
					TLatex *pt1tex = new TLatex(0.75, 0.9, pT1_label[0]);
					pt1tex->SetNDC();
				if(ibin == 1)
					pt1tex->Draw();
				}

			//SignalMinusRL[ibin][ibin2][ibin3] = (TH1D*)Signal[ibin][ibin2][ibin3]->Clone((TString)("signal_MinusRL"+CBin_strs[ibin]+"_"+TrkPtBin_strs[ibin3]));
			//SignalMinusRL[ibin][ibin2][ibin3]->Add(	PhiProj_R[ibin][ibin2][ibin3], -1);
			//SignalMinusRL[ibin][ibin2][ibin3]->SetAxisRange(-1.5, 1.5,"X");



				

		cc11->SaveAs("cc11_resultOnPhi_Leading_AjUnder22.png");	
	//##///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				cc12->cd(4*(ibin3+1)-ibin);
	
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
			//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Draw("LEGO2");   //testing that this is the S/ME from cc6

		//	yield_inc[ibin][ibin2][ibin3] =  (TH2D*) hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Clone((TString)("Yield_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));	

			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

			
			hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

			hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(100./me00[ibin][ibin2][ibin3]);	
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() * PhiBinWidth*EtaBinWidth));	
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Divide(hJetTrackMELeading[ibin][ibin2][ibin3]);

				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> CenterTitle();
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> CenterTitle();
				//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(0.8);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(0.8);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.05);


//			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Divide(hJetTrackMELeading[ibin][ibin2][ibin3]);
			yield_inc_proj[ibin][ibin2][ibin3] = (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionX((Form("ResultEta%d%d%d",ibin,ibin2,ibin3)),1, 50);	
			//yield_inc_proj[ibin][ibin2][ibin3] = (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionX((Form("ResultEta%d%d%d",ibin,ibin2,ibin3)),1, 25);	
					yield_inc_proj[ibin][ibin2][ibin3]->Rebin(4);
					yield_inc_proj[ibin][ibin2][ibin3]->Scale(EtaBinWidth/4.);   //
						
					//yield_inc_proj[ibin][ibin2][ibin3]->Draw();
						yield_inc_proj[ibin][ibin2][ibin3]->Fit("flat_fit","0","",-3.,-1.5);
					left_level = flat_fit->GetParameter(0);
							yield_inc_proj[ibin][ibin2][ibin3]->Fit("flat_fit","0","",1.5,3.);
						right_level =  flat_fit->GetParameter(0);
							bg_fit[ibin][ibin2][ibin3] = new TF1("bg_fit%d%d%d","[0]+x-x",-3.,3.);
							//double level = (left_level+right_level)/2.;
double level;
							 level = (left_level+right_level)/2.;
							 bg_fit[ibin][ibin2][ibin3]->SetParameter(0,(left_level+right_level)/2.);
					yield_inc_proj[ibin][ibin2][ibin3]->SetAxisRange(-3., 3., "X");
					yield_inc_proj[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-3., 3.);
					//yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(level- (3* yield_inc_proj[ibin][ibin2][ibin3]->GetRMS());
					//yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(level- (3* yield_inc_proj[ibin][ibin2][ibin3]->GetRMS()));
					//yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(level- 1000);
					//yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(level*1.2);
					//yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(level + (3* yield_inc_proj[ibin][ibin2][ibin3]->GetRMS() ));
					//yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(level + 1000);
					
					yield_inc_proj[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");
					yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(0.85 * level );
					yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(1.15 * level );

					if(ibin == 0 && ibin3 == 0)
					{
					yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(195. );					
					yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(175. );
					}

					if(ibin == 0 && ibin3 == 1)
					{
					yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(28. );
					yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(23. );	
					}

					if(ibin == 0 && ibin3 == 2)
					{
					yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(5. );
					yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(3. );	
					}
				
					
					if(ibin == 0 && ibin3 == 3)
					{
					yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(2.0 );
					yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(0.5 );	
					}

					
					if(ibin == 1 && ibin3 == 0)
					{
					yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(110. );
					yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(100. );	
					}

					if(ibin == 1 && ibin3 == 1)
					{
					yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(17. );
					yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(14. );	
					}

					
					if(ibin == 1 && ibin3 == 2)
					{
					yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(3. );
					yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(2. );	
					}

					
					if(ibin == 1 && ibin3 == 3)
					{
					yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(1.0 );
					yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(0.5 );	
					}


					if(ibin == 2 && ibin3 == 0)
					{
					yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(47. );
					yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(42. );	
					}

					
					if(ibin == 2 && ibin3 == 1)
					{
					yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(7.5 );
					yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(5.7 );	
					}
					
					
					if(ibin == 2 && ibin3 == 2)
					{
					yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(1.8 );
					yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(0.8 );	
					}

					
					if(ibin == 2 && ibin3 == 3)
					{
					yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(1. );
					yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(0. );	
					}

					
					if(ibin == 3 && ibin3 == 0)
					{
					yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(15. );
					yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(11. );	
					}


					if(ibin == 3 && ibin3 == 1)
					{
					yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(3.5 );
					yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(0.5 );	
					}

					
					if(ibin == 3 && ibin3 == 2)
					{
					yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(1.5 );
					yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(0. );	
					}

					
					if(ibin == 3 && ibin3 == 3)
					{
					yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(1. );
					yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(0. );	
					}
					//if(ibin==0&&ibin3==0)	yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(80.);
					//if(ibin==1&&ibin3==0)	yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(50.);
					//if(ibin==2&&ibin3==0)	yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(20.);
			yield_inc_proj[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(0.8);
			yield_inc_proj[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
			yield_inc_proj[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
			yield_inc_proj[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
			yield_inc_proj[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.05);


					yield_inc_proj[ibin][ibin2][ibin3]->Draw();
					bg_fit[ibin][ibin2][ibin3]->Draw("same");

						
				{
					TLatex *centtex = new TLatex(0.15,0.85,CBin_labels[ibin]);
					centtex->SetNDC();
					centtex->Draw();
					TLatex *pttex = new TLatex(0.15,0.78,TrkPtBin_labels[ibin3]);
					pttex->SetNDC();
					pttex->Draw();
					TLatex *Ajtex = new TLatex(0.75,0.85,AjBin_labels[0]);
					Ajtex->SetNDC();
					Ajtex->Draw();
					TLatex *Jettex = new TLatex(0.70, 0.80, JetLabel_labels[0]);
					Jettex->SetNDC();
					Jettex->Draw();
					TLatex *COMtex = new TLatex(0.15, 0.93, COM_label[0]);
					COMtex->SetNDC();
				//	if(ibin == 3)
					COMtex->Draw();
					TLatex *CMStex = new TLatex(0.15, 0.70, CMS_label[0]);
					CMStex->SetNDC();
				//	if(ibin == 2)
					CMStex->Draw();
					TLatex *Pbtex = new TLatex(0.65, 0.93, Pb_label[0]);
					Pbtex->SetNDC();
					Pbtex->Draw();

				}

//			cc12->SaveAs("cc12_RightBandLeftBandComparison_mmagnified_Leading_Under22.png");
			
				cc13->cd(4*(ibin3+1)-ibin);
	
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
			//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Draw("LEGO2");   //testing that this is the S/ME from cc6

		//	yield_inc[ibin][ibin2][ibin3] =  (TH2D*) hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Clone((TString)("Yield_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));	

			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

			
			hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

			hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(100./me00[ibin][ibin2][ibin3]);	
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() * PhiBinWidth*EtaBinWidth));	
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Divide(hJetTrackMELeading[ibin][ibin2][ibin3]);

				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> CenterTitle();
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> CenterTitle();
				//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(0.8);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(0.8);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.05);


//			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Divide(hJetTrackMELeading[ibin][ibin2][ibin3]);
			yield_inc_proj[ibin][ibin2][ibin3] = (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionX((Form("ResultEta%d%d%d",ibin,ibin2,ibin3)),1, 50);	
			//yield_inc_proj[ibin][ibin2][ibin3] = (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionX((Form("ResultEta%d%d%d",ibin,ibin2,ibin3)),1, 25);	
					yield_inc_proj[ibin][ibin2][ibin3]->Rebin(4);
					yield_inc_proj[ibin][ibin2][ibin3]->Scale(EtaBinWidth/4.);   //
						
					//yield_inc_proj[ibin][ibin2][ibin3]->Draw();
					cout << " ############################# yield RMS = " <<	yield_inc_proj[ibin][ibin2][ibin3]->GetRMS() << endl;
						yield_inc_proj[ibin][ibin2][ibin3]->Fit("flat_fit","0","",-3.,-1.5);
					left_level = flat_fit->GetParameter(0);
					left_level_error = flat_fit->GetParError(0);
					cout << "###############L############# Left Error = " << left_level_error << endl;
							yield_inc_proj[ibin][ibin2][ibin3]->Fit("flat_fit","0","",1.5,3.);
						right_level =  flat_fit->GetParameter(0);
						right_level_error =  flat_fit->GetParError(0);

							bg_fit[ibin][ibin2][ibin3] = new TF1("bg_fit%d%d%d","[0]+x-x",-3.,3.);
							//double level = (left_level+right_level)/2.;
							 level = (left_level+right_level)/2.;
							 bg_fit[ibin][ibin2][ibin3]->SetParameter(0,(left_level+right_level)/2.);
					yield_inc_proj[ibin][ibin2][ibin3]->SetAxisRange(-3., 3., "X");
					yield_inc_proj[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-3., 3.);
					//yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(level- (3* yield_inc_proj[ibin][ibin2][ibin3]->GetRMS());
					//yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(level- (3* yield_inc_proj[ibin][ibin2][ibin3]->GetRMS()));
					//yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(level- 1000);
					//yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(level*1.2);
					//yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(level + (3* yield_inc_proj[ibin][ibin2][ibin3]->GetRMS() ));
					//yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(level + 1000);
					
					yield_inc_proj[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");
					yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(3*right_level_error + level );
					yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(-3*right_level_error + level );
					//if(ibin==0&&ibin3==0)	yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(80.);
					//if(ibin==1&&ibin3==0)	yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(50.);
					//if(ibin==2&&ibin3==0)	yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(20.);
			yield_inc_proj[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(0.8);
			yield_inc_proj[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
			yield_inc_proj[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
			yield_inc_proj[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
			yield_inc_proj[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.05);


					yield_inc_proj[ibin][ibin2][ibin3]->Draw();
					bg_fit[ibin][ibin2][ibin3]->Draw("same");
				//	cout <<  " RMS = " << bg_fit[ibin][ibin2][ibin3]->GetRMS() << endl;
						
				{
					TLatex *centtex = new TLatex(0.15,0.85,CBin_labels[ibin]);
					centtex->SetNDC();
					centtex->Draw();
					TLatex *pttex = new TLatex(0.15,0.78,TrkPtBin_labels[ibin3]);
					pttex->SetNDC();
					pttex->Draw();
					TLatex *Ajtex = new TLatex(0.75,0.85,AjBin_labels[0]);
					Ajtex->SetNDC();
					Ajtex->Draw();
					TLatex *Jettex = new TLatex(0.70, 0.80, JetLabel_labels[0]);
					Jettex->SetNDC();
					Jettex->Draw();

				}

//				cc13->SaveAs("cc13_RightBandLeftBandComparison_mmagnified_Leading_Under22.png");

///////////////////////
				cc14->cd(4*(ibin3 + 1)- ibin);
				hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

				hJetTrackMELeading[ibin][ibin2][ibin3]->Draw("LEGO2");
	{
					TLatex *centtex = new TLatex(0.15,0.85,CBin_labels[ibin]);
					centtex->SetNDC();
					centtex->Draw();
					TLatex *pttex = new TLatex(0.15,0.78,TrkPtBin_labels[ibin3]);
					pttex->SetNDC();
					pttex->Draw();
					TLatex *Ajtex = new TLatex(0.75,0.85,AjBin_labels[0]);
					Ajtex->SetNDC();
					Ajtex->Draw();
					TLatex *Jettex = new TLatex(0.70, 0.80, JetLabel_labels[0]);
					Jettex->SetNDC();
					Jettex->Draw();
		}
				cc15->cd(4*(ibin3+1)-ibin);
				
				me_projeta[ibin][ibin2][ibin3] =  hJetTrackMELeading[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_eta_cc15_%d%d%d",ibin,ibin2,ibin3),1,100); //what is the 1,100 for?
				me_projeta[ibin][ibin2][ibin3]->Scale(1./100.);
				me_projeta[ibin][ibin2][ibin3]->Draw() ;
				me_projeta[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("Delta eta");
				NEtaBins =  me_projeta[ibin][ibin2][ibin3]->GetNbinsX() ;
				me_projeta[ibin][ibin2][ibin3] -> Fit("pol0","","",-.2,.2);
				me00[ibin][ibin2][ibin3] = me_projeta[ibin][ibin2][ibin3]->GetFunction("pol0")->GetParameter(0);		
				//me_projeta[ibin][ibin2][ibin3]->SetAxisRange(0, 500,"Y");
				me_projeta[ibin][ibin2][ibin3]->SetAxisRange(-2, 2,"X");
				me_projeta[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-2.,2.);	
				//me_projeta[ibin][ibin2][ibin3]->SetMinimum(100);
				me_projeta[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");
			//	me_projeta[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(0.8);
				me_projeta[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(0.8);
				me_projeta[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				me_projeta[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				me_projeta[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.05);
				me_projeta[ibin][ibin2][ibin3]->GetXaxis() -> CenterTitle();

				//me_projeta[ibin][ibin2][ibin3]->Scale( 1.1*me_projeta[ibin][ibin2][ibin3]->GetMaximum());
				me_projeta[ibin][ibin2][ibin3]->DrawCopy() ;
		
				//gPad->Update();
					
			//	if (((4*ibin3 + 1) - ibin) < 4)
			//		return 1;	
				// For now we have issue setting different values for the axes

				//yield_inc[ibin][ibin2][ibin3]->Draw("Lego2");	
				//if(ibin<3)
				{
					TLatex *centtex = new TLatex(0.15,0.85,CBin_labels[ibin]);
					centtex->SetNDC();
					centtex->Draw();
					TLatex *pttex = new TLatex(0.15,0.78,TrkPtBin_labels[ibin3]);
					pttex->SetNDC();
					pttex->Draw();
					TLatex *Ajtex = new TLatex(0.76,0.85,AjBin_labels[0]);
					Ajtex->SetNDC();
					Ajtex->Draw();
					TLatex *Jettex = new TLatex(0.715, 0.80, JetLabel_labels[0]);
					Jettex->SetNDC();
					Jettex->Draw();

				}
			//	if(ibin==3){
			//		TLatex *centtex = new TLatex(0.15,0.85,CBin_labels[ibin]);
			//		centtex->SetNDC();
			//		centtex->Draw();
			//		TLatex *pttex = new TLatex(0.15,0.78,TrkPtBin_labels[ibin3]);
			//		pttex->SetNDC();
			//		pttex->Draw();
			//		TLatex *Ajtex = new TLatex(0.76,0.85,AjBin_labels[0]);
			//		Ajtex->SetNDC();
			//		Ajtex->Draw();
			//		TLatex *Jettex = new TLatex(0.70, 0.80, JetLabel_labels[0]);
			//		Jettex->SetNDC();
			//		Jettex->Draw();
			//	

			//	}
//	cc14->SaveAs("cc14_ME_EtaProj_Leading_Under22.png");

			cc16->cd(4*(ibin3+1) - ibin);

				hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


				me_projeta[ibin][ibin2][ibin3] =  hJetTrackMELeading[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_eta_cc16_%d%d%d",ibin,ibin2,ibin3),1,100); //what is the 1,100 for?
				me_projeta[ibin][ibin2][ibin3]->Scale(1./100.);
				me_projeta[ibin][ibin2][ibin3]->Draw() ;
				me_projeta[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("Delta eta");
				NEtaBins =  me_projeta[ibin][ibin2][ibin3]->GetNbinsX() ;
			//	me_projeta[ibin][ibin2][ibin3] -> Fit("pol0","","",-.2,.2);
			//	me00[ibin][ibin2][ibin3] = me_projeta[ibin][ibin2][ibin3]->GetFunction("pol0")->GetParameter(0);	
	//##///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////				
	//##//
	//##//				cout << "Eta Flat Fitting Parameter value = " << me00[ibin][ibin2][ibin3] << endl;
				//me_projeta[ibin][ibin2][ibin3] -> Fit("pol0","0","",-.2,.2);
				//me00[ibin][ibin2][ibin3] = me_projeta[ibin][ibin2][ibin3]->GetFunction("pol0")->GetParameter(0);	
				me_projeta[ibin][ibin2][ibin3]->Scale (1./me00[ibin][ibin2][ibin3]);
				//me_projeta[ibin][ibin2][ibin3]->Scale (10./me00[ibin][ibin2][ibin3]);
	//##//				me_projeta[ibin][ibin2][ibin3]->Scale (1./me00[ibin][ibin2][ibin3]);
			//		hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1./me00[ibin][ibin2][ibin3]);
	//##//				hJetTrackMELeading[ibin][ibin2][ibin3]->Scale (NEtaBins /me00[ibin][ibin2][ibin3]);
	//##//				//cout<< "me00 = " << me00 <<endl;
				
				me_projeta[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");
				me_projeta[ibin][ibin2][ibin3]->SetAxisRange(-2,2);
				me_projeta[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-2.,2.);	
				me_projeta[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(0.8);
				me_projeta[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				me_projeta[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				me_projeta[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.06);
				me_projeta[ibin][ibin2][ibin3]->GetXaxis() -> CenterTitle();

				me_projeta[ibin][ibin2][ibin3]->Draw();
	
				{
					TLatex *centtex = new TLatex(0.15,0.85,CBin_labels[ibin]);
					centtex->SetNDC();
					centtex->Draw();
					TLatex *pttex = new TLatex(0.15,0.78,TrkPtBin_labels[ibin3]);
					pttex->SetNDC();
					pttex->Draw();
					TLatex *Ajtex = new TLatex(0.76,0.85,AjBin_labels[0]);
					Ajtex->SetNDC();
					Ajtex->Draw();
					TLatex *Jettex = new TLatex(0.715, 0.80, JetLabel_labels[0]);
					Jettex->SetNDC();
					Jettex->Draw();

					}

//					cc16->SaveAs("cc16_MEO_EtaProj_ScaledTo1_Leading_Under22.png");

////////////////////////////////////////
				
				cc17->cd(4*(ibin3+1)-ibin);
				hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

			hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1./me00[ibin][ibin2][ibin3]);	
			//hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(0.01/me00[ibin][ibin2][ibin3]);	
			
		hJetTrackMELeading[ibin][ibin2][ibin3]->Draw("LEGO2");
					{
					TLatex *centtex = new TLatex(0.15,0.85,CBin_labels[ibin]);
					centtex->SetNDC();
					centtex->Draw();
					TLatex *pttex = new TLatex(0.15,0.78,TrkPtBin_labels[ibin3]);
					pttex->SetNDC();
					pttex->Draw();
					TLatex *Ajtex = new TLatex(0.75,0.85,AjBin_labels[0]);
					Ajtex->SetNDC();
					Ajtex->Draw();
					TLatex *Jettex = new TLatex(0.70, 0.80, JetLabel_labels[0]);
					Jettex->SetNDC();
					Jettex->Draw();

				}
	
//				cc17->SaveAs("cc17_MEOverMEO_Leading_Under22.png");


				cc18->cd(4*(ibin3+1)-ibin);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

				hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
				
	
				me_projeta[ibin][ibin2][ibin3] =  hJetTrackMELeading[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_eta_cc18_%d%d%d",ibin,ibin2,ibin3),1,100); //what is the 1,100 for?
				me_projeta[ibin][ibin2][ibin3]->Scale(1./100.);
				//me_projeta[ibin][ibin2][ibin3]->Draw() ;
				me_projeta[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("Delta eta");
				NEtaBins =  me_projeta[ibin][ibin2][ibin3]->GetNbinsX() ;
				me_projeta[ibin][ibin2][ibin3] -> Fit("pol0","","",-.2,.2);
				me00[ibin][ibin2][ibin3] = me_projeta[ibin][ibin2][ibin3]->GetFunction("pol0")->GetParameter(0);	
				hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1./me00[ibin][ibin2][ibin3]);	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Divide(hJetTrackMELeading[ibin][ibin2][ibin3]);	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() * EtaBinWidth * PhiBinWidth ));
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-3.,3.);	
	
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Draw("LEGO2");


	{
					TLatex *centtex = new TLatex(0.15,0.85,CBin_labels[ibin]);
					centtex->SetNDC();
					centtex->Draw();
					TLatex *pttex = new TLatex(0.15,0.78,TrkPtBin_labels[ibin3]);
					pttex->SetNDC();
					pttex->Draw();
					TLatex *Ajtex = new TLatex(0.75,0.85,AjBin_labels[0]);
					Ajtex->SetNDC();
					Ajtex->Draw();
					TLatex *Jettex = new TLatex(0.70, 0.80, JetLabel_labels[0]);
					Jettex->SetNDC();
					Jettex->Draw();

				}
	
      			cc18->SaveAs("cc18_MEOverMEO_JetScaled_perEtaperPhi_Leading_AjUnder22.png");


				cc19->cd(4*(ibin3+1)-ibin);
		hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

		hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
				

		hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1./me00[ibin][ibin2][ibin3]);	
		hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Divide(hJetTrackMELeading[ibin][ibin2][ibin3]);	
		hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() * EtaBinWidth * PhiBinWidth ));	
		
		me_projeta[ibin][ibin2][ibin3] =  hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_eta_cc19_%d%d%d",ibin,ibin2,ibin3),18,33); //what is the 1,100 for?
		me_projeta[ibin][ibin2][ibin3]->Scale(PhiBinWidth);
		me_projeta[ibin][ibin2][ibin3]->Rebin(2);
		me_projeta[ibin][ibin2][ibin3]->Scale(1./2.);	
		me_projeta[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-3.,3.);
		me_projeta[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");		
		me_projeta[ibin][ibin2][ibin3]->Draw();
//				cc19->SaveAs("cc19_MEOverMEO_JetScaled_EtaProj_PhiNarrow_Leading_Under22.png");

				cc20->cd(4*(ibin3+1)-ibin);
	
				
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

				hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
			

				me_projeta[ibin][ibin2][ibin3] =  hJetTrackMELeading[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_eta_cc20_%d%d%d",ibin,ibin2,ibin3),1,100); //what is the 1,100 for?
				me_projeta[ibin][ibin2][ibin3]->Scale(1./100.);
				me_projeta[ibin][ibin2][ibin3]->Draw() ;
				me_projeta[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("Delta eta");
				NEtaBins =  me_projeta[ibin][ibin2][ibin3]->GetNbinsX() ;
				me_projeta[ibin][ibin2][ibin3] -> Fit("pol0","","",-.2,.2);
				me00[ibin][ibin2][ibin3] = me_projeta[ibin][ibin2][ibin3]->GetFunction("pol0")->GetParameter(0);	

				hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1./me00[ibin][ibin2][ibin3]);
				//hJetTrackMELeading[ibin][ibin2][ibin3]->Draw("Lego2");
				//return 20;	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Divide(hJetTrackMELeading[ibin][ibin2][ibin3]);	
				
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / (  Aj[ibin]->GetEntries() * EtaBinWidth * PhiBinWidth ));	
	
				PhiProj_L[ibin][ibin2][ibin3] = (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionY(Form("ResultPhi_L%d%d%d",ibin,ibin2,ibin3),21, 35);
				PhiProj_R[ibin][ibin2][ibin3] =  (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionY(Form("ResultPhi_R%d%d%d",ibin,ibin2,ibin3),66, 80);	
				PhiProj_R[ibin][ibin2][ibin3]->Add(PhiProj_L[ibin][ibin2][ibin3], 1);
				PhiProj_R[ibin][ibin2][ibin3]->Rebin(4);
				PhiProj_R[ibin][ibin2][ibin3]->Scale(1./4.);		
				PhiProj_R[ibin][ibin2][ibin3]->Scale(1./30.);
				PhiProj_R[ibin][ibin2][ibin3]->Draw();

//				cc20->SaveAs("cc20_MEOverMEO_JetScaled_PhiProj_EtaSide_perEtaperPhi_Leading_Under22.png");

				cc21->cd(4*(ibin3+1)-ibin);


			SignalMinusRL[ibin][ibin2][ibin3] = (TH1D*)Signal[ibin][ibin2][ibin3]->Clone((TString)("signal_MinusRL"+CBin_strs[ibin]+"_"+TrkPtBin_strs[ibin3]));
			SignalMinusRL[ibin][ibin2][ibin3]->Add(	PhiProj_R[ibin][ibin2][ibin3], -1);
			SignalMinusRL[ibin][ibin2][ibin3]->SetAxisRange(-1.5, 1.5,"X");
			SignalMinusRL[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("1/N_{jet} dN/(d#Delta#phi dp_{T})");
		//	SignalMinusRL[ibin][ibin2][ibin3]->GetYaxis()->SetTitleSize(0.06);
			SignalMinusRL[ibin][ibin2][ibin3]->GetYaxis()->CenterTitle();
			SignalMinusRL[ibin][ibin2][ibin3]->GetYaxis()->SetTitleOffset(0.6);
			
		//	SignalMinusRL[ibin][ibin2][ibin3]->Scale(1/(3.0 *  Aj[ibin]->GetEntries() * PhiBinWidth));

			if(ibin3==0)	SignalMinusRL[ibin][ibin2][ibin3]->SetMinimum(-2.);
			if(ibin3==0)	SignalMinusRL[ibin][ibin2][ibin3]->SetMaximum(8.);


			SignalMinusRL[ibin][ibin2][ibin3]->SetMarkerStyle(20);
			SignalMinusRL[ibin][ibin2][ibin3]->SetMarkerColor(kBlue);
			SignalMinusRL[ibin][ibin2][ibin3]->Draw();

//				cc21->SaveAs("cc21_MEOverMEO_JetScaled_PhiProj_EtaSide_perEtaperPhi_Leading_Under22.png");

				cc22->cd(4*(ibin3+1)-ibin);

	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

				PhiBinWidth =  me_projphi[ibin][ibin2][ibin3]->GetBinWidth(1) ;
				EtaBinWidth =  hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->GetBinWidth(1) ;
				PhiBinWidth =  hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->GetBinWidth(1) ;
				cout << "Eta Bint Width = " << 	 hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->GetBinWidth(1) << endl;


				//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / (3.0 *  Aj[ibin]->GetEntries() * PhiBinWidth));
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / (  Aj[ibin]->GetEntries() ));
			//	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / (EtaBinWidth *  Aj[ibin]->GetEntries() * PhiBinWidth));
			//	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / (3.0 *  Aj[ibin]->GetEntries()));
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> CenterTitle();
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> CenterTitle();
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(1.6);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(1.6);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->SetAxisRange(-3., 3.,"X");
		
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-3.,3.);	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.06);
	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Draw("LEGO2");	



				{
					TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
					centtex->SetNDC();
					centtex->Draw();
					TLatex *pttex = new TLatex(0.05,0.83,TrkPtBin_labels[ibin3]);
					pttex->SetNDC();
					pttex->Draw();
					TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
					Ajtex->SetNDC();
				//	Ajtex->Draw();
					TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
					Jettex->SetNDC();
				//	Jettex->Draw();


				}
			


			
//				cc22->SaveAs("cc22_S_Scaled_NotPerEtaNotPerPhi.png");


				cc23->cd(4*(ibin3+1)-ibin);
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


			hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


			//hJetTrackMELeading[ibin][ibin2][ibin3]->Draw("LEGO2");	
			//return 23;


			hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1./me00[ibin][ibin2][ibin3]);	
			//hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(100 /me00[ibin][ibin2][ibin3]);
		
	//	hJetTrackMELeading[ibin][ibin2][ibin3]->Draw("LEGO2");	
	//		return 23;
		//	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() ));
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() ));
//			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() * PhiBinWidth*EtaBinWidth));	
				
			//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] -> Draw("LEGO2");
			//return 23;
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Divide(hJetTrackMELeading[ibin][ibin2][ibin3]);
			//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] -> Draw("LEGO2");
			//return 23;	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> CenterTitle();
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> CenterTitle();
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(1.6);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(1.6);
				//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->SetAxisRange(-3., 3.,"X");
				//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->SetAxisRange(-3., 3.,"X");

				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-3.,3.);	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.06);



	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Draw("LEGO2");	

	//				if(ibin<3){
					TLatex *centtex23 = new TLatex(0.05,0.9,CBin_labels[ibin]);
					centtex23->SetNDC();
					centtex23->Draw();
					TLatex *pttex23 = new TLatex(0.05,0.83,TrkPtBin_labels[ibin3]);
					pttex23->SetNDC();
					pttex23->Draw();
					TLatex *Ajtex23 = new TLatex(0.80,0.90,AjBin_labels[0]);
					Ajtex23->SetNDC();
					Ajtex23->Draw();
					TLatex *Jettex23 = new TLatex(0.75, 0.85, JetLabel_labels[0]);
					Jettex23->SetNDC();
					Jettex23->Draw();


	//			}
	//			if(ibin==3){
	//				TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	//				centtex->SetNDC();
	//				centtex->Draw();
	//				TLatex *pttex = new TLatex(0.15,0.83,TrkPtBin_labels[ibin3]);
	//				pttex->SetNDC();
	//				pttex->Draw();
	//				TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
	//				Ajtex->SetNDC();
	//				Ajtex->Draw();
	//				TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
	//				Jettex->SetNDC();
	//				Jettex->Draw();


	//			}
	
//				cc23->SaveAs("cc23_S_Scaled_overMEoverM0_NotScaled_NotPerEtaNotPerPhi.png");
				
			
				cc24->cd(4*(ibin3+1)-ibin);
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


			hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


			//hJetTrackMELeading[ibin][ibin2][ibin3]->Draw("LEGO2");	
		

			me_projeta[ibin][ibin2][ibin3] =  hJetTrackMELeading[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_eta_cc24_%d%d%d",ibin,ibin2,ibin3),1,100);	
			me_projeta[ibin][ibin2][ibin3]->Scale(1. / 100.);	
			me_projeta[ibin][ibin2][ibin3] -> Fit("pol0","","",-.2,.2);
			hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1./me00[ibin][ibin2][ibin3]);	
			//hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(100 /me00[ibin][ibin2][ibin3]);
		
	//	hJetTrackMELeading[ibin][ibin2][ibin3]->Draw("LEGO2");	
	//		return 23;
		//	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() ));
		//	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() ));
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() * PhiBinWidth * EtaBinWidth));	
				
			//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] -> Draw("LEGO2");
			//return 23;
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Divide(hJetTrackMELeading[ibin][ibin2][ibin3]);
			//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] -> Draw("LEGO2");
			//return 23;	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> CenterTitle();
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> CenterTitle();
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(1.6);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(1.6);
		      	//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->SetAxisRange(-3., 3.,"X");
		      	//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->SetAxisRange(-3., 3.,"X");

				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-3.,3.);	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.06);

	
//				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Draw("LEGO2");	
				//return 24;
				// make a clone of the output
				result[ibin][ibin2][ibin3] = (TH2D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Clone((TString)("result_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])); 
				
			//	result[ibin][ibin2][ibin3]->Draw("LEGO2");
	
				 llimiteta = result[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-3.0 + 0.0001); 
				 rlimiteta = result[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-1.5 - 0.0001); 

//				cout << "llimiteta = " << llimiteta << " rlimiteta = " << rlimiteta << endl;
//return 24;
				background_left[ibin][ibin2][ibin3] = (TH1D*)result[ibin][ibin2][ibin3]->ProjectionY((TString)("LeftSideband"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]), llimiteta, rlimiteta);

				//background_left[ibin][ibin2][ibin3] -> Draw();
		
				llimiteta =  result[ibin][ibin2][ibin3]->GetXaxis()->FindBin(1.5 + 0.0001); 
				rlimiteta =  result[ibin][ibin2][ibin3]->GetXaxis()->FindBin(3.0 - 0.0001); 
							
				cout << "For the right Eta band, llimiteta = " << llimiteta << " rlimiteta = " << rlimiteta << endl;
				

				background_proj[ibin][ibin2][ibin3] = (TH1D*)result[ibin][ibin2][ibin3]->ProjectionY((TString)("Sideband"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]), llimiteta, rlimiteta);

				background_proj[ibin][ibin2][ibin3]->Add(	background_left[ibin][ibin2][ibin3] )	;
				background_proj[ibin][ibin2][ibin3]->Scale(1. / (2*(rlimiteta - llimiteta+1)))	;

//				background_proj[ibin][ibin2][ibin3]->Draw();

         			background[ibin][ibin2][ibin3] = (TH2D*)result[ibin][ibin2][ibin3]->Clone((TString)("Background"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	
//				background[ibin][ibin2][ibin3]->Draw("LEGO2");		
		cout << "How many Y bins in result = " << result[ibin][ibin2][ibin3]->GetNbinsY() << endl;		

		for (int k = 1; k <= result[ibin][ibin2][ibin3]->GetNbinsY(); k++){
		temp1 = background_proj[ibin][ibin2][ibin3]->GetBinContent(k);	
		err1  = background_proj[ibin][ibin2][ibin3]->GetBinError(k);

		for (int m = 1; m <= result[ibin][ibin2][ibin3]->GetNbinsX(); m++){
		//cout << "temp1 = " << temp1 << endl; 
		background[ibin][ibin2][ibin3]->SetBinContent(m, k, temp1);
		background[ibin][ibin2][ibin3]->SetBinError(m, k, err1);
		//cout << "I am working " << endl;
		}		
		}

		result[ibin][ibin2][ibin3]->Add(background[ibin][ibin2][ibin3], -1.);
		background[ibin][ibin2][ibin3]->Draw("LEGO2");
		//result[ibin][ibin2][ibin3]->Draw("LEGO2");
		
	//				if(ibin<3){
					TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
					centtex->SetNDC();
					centtex->Draw();
					TLatex *pttex = new TLatex(0.05,0.83,TrkPtBin_labels[ibin3]);
					pttex->SetNDC();
					pttex->Draw();
					TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
					Ajtex->SetNDC();
					Ajtex->Draw();
					TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
					Jettex->SetNDC();
					Jettex->Draw();


	//			}
	//			if(ibin==3){
	//				TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	//				centtex->SetNDC();
	//				centtex->Draw();
	//				TLatex *pttex = new TLatex(0.15,0.83,TrkPtBin_labels[ibin3]);
	//				pttex->SetNDC();
	//				pttex->Draw();
	//				TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
	//				Ajtex->SetNDC();
	//				Ajtex->Draw();
	//				TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
	//				Jettex->SetNDC();
	//				Jettex->Draw();


	//			}
	
		
				

				cc24->SaveAs("CC24_S_Scaled_overMEoverM0_NotScaled_PerEtaPerPhi_ProjectionY_etaSideBand_AjUnder22.png");	
	//hJetBackgroundLeading[ibin][ibin2][ibin3]->Fill(PhiProj_R[ibin][ibin2][ibin3], 1);

	cc25->cd(4*(ibin3+1)-ibin);
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


			hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


			//hJetTrackMELeading[ibin][ibin2][ibin3]->Draw("LEGO2");	
			//return 23;


			hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1./me00[ibin][ibin2][ibin3]);	
			//hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(100 /me00[ibin][ibin2][ibin3]);
		
	//	hJetTrackMELeading[ibin][ibin2][ibin3]->Draw("LEGO2");	
	//		return 23;
		//	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() ));
		//	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() ));
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() * PhiBinWidth * EtaBinWidth));	
				
			//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] -> Draw("LEGO2");
			//return 23;
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Divide(hJetTrackMELeading[ibin][ibin2][ibin3]);
			//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] -> Draw("LEGO2");
			//return 23;	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> CenterTitle();
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> CenterTitle();
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(1.6);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(1.6);
		      	//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->SetAxisRange(-3., 3.,"X");
		      	//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->SetAxisRange(-3., 3.,"X");

				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-3.,3.);	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.06);

					
				//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Draw("LEGO2");	
				//return 24;
				// make a clone of the output
				result[ibin][ibin2][ibin3] = (TH2D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Clone((TString)("result_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])); 

				 llimiteta = result[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-3.+0.0001); 
				 rlimiteta = result[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-1.5 - 0.0001); 

				cout << "llimiteta = " << llimiteta << " rlimiteta = " << rlimiteta << endl;
//return 24;
				background_left[ibin][ibin2][ibin3] = (TH1D*)result[ibin][ibin2][ibin3]->ProjectionY((TString)("LeftSideband"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]), llimiteta, rlimiteta);

				llimiteta =  result[ibin][ibin2][ibin3]->GetXaxis()->FindBin(1.5 + 0.0001); 
				rlimiteta =  result[ibin][ibin2][ibin3]->GetXaxis()->FindBin(3.0 - 0.0001); 

				background_proj[ibin][ibin2][ibin3] = (TH1D*)result[ibin][ibin2][ibin3]->ProjectionY((TString)("Sideband"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]), llimiteta, rlimiteta);

				background_proj[ibin][ibin2][ibin3]->Add(	background_left[ibin][ibin2][ibin3] )	;
				background_proj[ibin][ibin2][ibin3]->Scale(1. / (2*(rlimiteta - llimiteta+1)))	;

//				background_proj[ibin][ibin2][ibin3]->Draw();

			background[ibin][ibin2][ibin3] = (TH2D*)result[ibin][ibin2][ibin3]->Clone((TString)("Background"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
			

		for (int k = 1; k <= result[ibin][ibin2][ibin3]->GetNbinsY(); k++){
		temp1 = background_proj[ibin][ibin2][ibin3]->GetBinContent(k);	
		err1  = background_proj[ibin][ibin2][ibin3]->GetBinError(k);

		for (int m = 1; m <= result[ibin][ibin2][ibin3]->GetNbinsX(); m++){
		//cout << "temp1 = " << temp1 << endl; 
		background[ibin][ibin2][ibin3]->SetBinContent(m, k, temp1);
		background[ibin][ibin2][ibin3]->SetBinError(m, k, err1);
		//cout << "I am working " << endl;
		}		
		}

		result[ibin][ibin2][ibin3]->Add(background[ibin][ibin2][ibin3], -1.);
		//background[ibin][ibin2][ibin3]->Draw("LEGO2");
		result[ibin][ibin2][ibin3]->Draw("LEGO2");
		
	//				if(ibin<3){
				TLatex *centtex1 = new TLatex(0.05,0.9,CBin_labels[ibin]);
					centtex1->SetNDC();
					centtex1->Draw();
					//TLatex *pttex1 = new TLatex(0.05,0.83,TrkPtBin_labels[ibin3]);
					//pttex1->SetNDC();
					pttex->Draw();
//					TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
//					Ajtex->SetNDC();
					Ajtex->Draw();
//					TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
//					Jettex->SetNDC();
					Jettex->Draw();


	//			}
	//			if(ibin==3){
	//				TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	//				centtex->SetNDC();
	//				centtex->Draw();
	//				TLatex *pttex = new TLatex(0.15,0.83,TrkPtBin_labels[ibin3]);
	//				pttex->SetNDC();
	//				pttex->Draw();
	//				TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
	//				Ajtex->SetNDC();
	//				Ajtex->Draw();
	//				TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
	//				Jettex->SetNDC();
	//				Jettex->Draw();


	//			}
	
		
				

				cc25->SaveAs("CC25_S_Scaled_overMEoverM0_NotScaled_PerEtaPerPhi_ProjectionY_etaSideBand_AjUnder22.png");	


					cc26->cd(4*(ibin3+1)-ibin);
				cc26->SetRightMargin(0.01);

			hJetTrackSignalBackgroundLeading_cc26[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_cc26_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
			
			//hJetTrackSignalBackgroundLeading_cc26[ibin][ibin2][ibin3]->Draw("LEGO2");
			hJetTrackMELeading_cc26[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_cc26_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));




			hJetTrackMELeading_cc26[ibin][ibin2][ibin3]->Scale(1./me00[ibin][ibin2][ibin3]);
//			hJetTrackMELeading_cc26[ibin][ibin2][ibin3]->Draw("LEGO2");	
				
			//hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(100 /me00[ibin][ibin2][ibin3]);
		
	//	hJetTrackMELeading[ibin][ibin2][ibin3]->Draw("LEGO2");	
	//		return 23;
		//	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() ));
		//	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() ));
			hJetTrackSignalBackgroundLeading_cc26[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() * PhiBinWidth * EtaBinWidth));	
				
		//	hJetTrackSignalBackgroundLeading_cc26[ibin][ibin2][ibin3] -> Draw("LEGO2");
			//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] -> Draw("LEGO2");
			//return 23;
			hJetTrackSignalBackgroundLeading_cc26[ibin][ibin2][ibin3]->Divide(hJetTrackMELeading_cc26[ibin][ibin2][ibin3]);

					
			//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] -> Draw("LEGO2");
			//return 23;	
			hJetTrackSignalBackgroundLeading_cc26[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
			hJetTrackSignalBackgroundLeading_cc26[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
			hJetTrackSignalBackgroundLeading_cc26[ibin][ibin2][ibin3]->GetXaxis() -> CenterTitle();
			hJetTrackSignalBackgroundLeading_cc26[ibin][ibin2][ibin3]->GetYaxis() -> CenterTitle();
			hJetTrackSignalBackgroundLeading_cc26[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(1.6);
			hJetTrackSignalBackgroundLeading_cc26[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(1.6);
			hJetTrackSignalBackgroundLeading_cc26[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-3.,3.);				
			hJetTrackSignalBackgroundLeading_cc26[ibin][ibin2][ibin3] -> Draw("LEGO2");

		      	//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->SetAxisRange(-3., 3.,"X");
		      	//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->SetAxisRange(-3., 3.,"X");

				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-3.,3.);	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.06);

	
				//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Draw("LEGO2");	
				//return 24;
				// make a clone of the output
				result_cc26[ibin][ibin2][ibin3] = (TH2D*)hJetTrackSignalBackgroundLeading_cc26[ibin][ibin2][ibin3]->Clone((TString)("result_cc26"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])); 

				 llimiteta = result[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-3.0 + 0.0001); 
				 rlimiteta = result[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-1.5 - 0.0001); 

				cout << "llimiteta = " << llimiteta << " rlimiteta = " << rlimiteta << endl;
//return 24;
				background_left_cc26[ibin][ibin2][ibin3] = (TH1D*)result_cc26[ibin][ibin2][ibin3]->ProjectionY((TString)("LeftSideband_cc26"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]), llimiteta, rlimiteta);
				
				llimiteta =  result[ibin][ibin2][ibin3]->GetXaxis()->FindBin(1.5 + 0.0001); 
				rlimiteta =  result[ibin][ibin2][ibin3]->GetXaxis()->FindBin(3.0 - 0.0001); 

				cout << "llimiteta = " << llimiteta << " rlimiteta = " << rlimiteta << endl;				

				background_proj_cc26[ibin][ibin2][ibin3] = (TH1D*)result_cc26[ibin][ibin2][ibin3]->ProjectionY((TString)("Sideband_cc26"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]), llimiteta, rlimiteta);

				
				background_proj_cc26[ibin][ibin2][ibin3]->Draw();

				std::cout << "background_proj_cc26: Number of bins (size) = " << background_proj_cc26[ibin][ibin2][ibin3]->GetSize() << std::endl; 
				std::cout << "background_left_cc26: Number of bins (size) = " << background_left_cc26[ibin][ibin2][ibin3]->GetSize() << std::endl; 
			
//return 26;	
				background_proj_cc26[ibin][ibin2][ibin3]->Add(	background_left_cc26[ibin][ibin2][ibin3] )	;
				background_proj_cc26[ibin][ibin2][ibin3]->Scale(1. / (2*(rlimiteta - llimiteta+1)))	;

//				background_proj_cc26[ibin][ibin2][ibin3]->Draw();

        	background_cc26[ibin][ibin2][ibin3] = (TH2D*)result_cc26[ibin][ibin2][ibin3]->Clone((TString)("Background_cc26"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
			

		for (int k = 1; k <= result_cc26[ibin][ibin2][ibin3]->GetNbinsY(); k++){
		temp1 = background_proj_cc26[ibin][ibin2][ibin3]->GetBinContent(k);	
		err1  = background_proj_cc26[ibin][ibin2][ibin3]->GetBinError(k);

		for (int m = 1; m <= result_cc26[ibin][ibin2][ibin3]->GetNbinsX(); m++){
		//cout << "temp1 = " << temp1 << endl; 
		background_cc26[ibin][ibin2][ibin3]->SetBinContent(m, k, temp1);
		background_cc26[ibin][ibin2][ibin3]->SetBinError(m, k, err1);
		//cout << "I am working " << endl;
		}		
		}
//return 26;


		std::cout << "result_cc26: Number of bins (size) = "     << result_cc26[ibin][ibin2][ibin3]->GetSize() << std::endl; 
		std::cout << "background_cc26: Number of bins (size) = " << background_cc26[ibin][ibin2][ibin3]->GetSize() << std::endl; 
			

		result_cc26[ibin][ibin2][ibin3]->Add(background_cc26[ibin][ibin2][ibin3], -1.);
		//background[ibin][ibin2][ibin3]->Draw("LEGO2");
//return 26;
		llimiteta = result_cc26[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-1.5 + 0.0001); 
		rlimiteta = result_cc26[ibin][ibin2][ibin3]->GetXaxis()->FindBin(1.5 - 0.0001); 
	

        	//Signal[ibin][ibin2][ibin3] =    (TH1D*)result[ibin][ibin2][ibin3]->ProjectionY(Form("Phi_Signal%d%d%d",ibin,ibin2,ibin3),36, 65);	
        	Signal_cc26[ibin][ibin2][ibin3] =    (TH1D*)result_cc26[ibin][ibin2][ibin3]->ProjectionY(Form("Phi_Signal_cc26_%d%d%d",ibin,ibin2,ibin3),llimiteta, rlimiteta);	
			//	double MaxSignal[6][6];
			//	MaxSignal[ibin][ibin3] = Signal[ibin][ibin2][ibin3]->GetMaximum();
			//	Signal[ibin][ibin2][ibin3]->SetMaximum(MaxSignal[ibin][ibin3]+30);
			
				Signal_cc26[ibin][ibin2][ibin3]->SetMarkerStyle(20);
				Signal_cc26[ibin][ibin2][ibin3]->SetMarkerColor(kRed);
				Signal_cc26[ibin][ibin2][ibin3]->Rebin(2);
				Signal_cc26[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#phi");
				Signal_cc26[ibin][ibin2][ibin3]->Scale(EtaBinWidth/2.);
				Signal_cc26[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(0.8);
				Signal_cc26[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(0.8);
				Signal_cc26[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				Signal_cc26[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				Signal_cc26[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(0.8);
				Signal_cc26[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				Signal_cc26[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
				Signal_cc26[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				Signal_cc26[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.05);
				Signal_cc26[ibin][ibin2][ibin3]->GetXaxis()-> SetRangeUser(-1.5,1.5);
				Signal_cc26[ibin][ibin2][ibin3]->GetYaxis()-> SetRangeUser(-0.5, 7);
				Signal_cc26[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#frac{1}{N_{jet}}#frac{d^{2}N}{ d#Delta#phi dp_{T}}");
				Signal_cc26[ibin][ibin2][ibin3]->SetMaximum(TMath::Max(Signal_cc26[ibin][ibin2][ibin3]->GetMaximum()+5.0, Signal_cc26[ibin][ibin2][ibin3]->GetMaximum())*1.05);
				
				Signal_cc26[ibin][ibin2][ibin3]->Draw();


	//	result[ibin][ibin2][ibin3]->Draw("LEGO2");
		
					if(ibin<3){
				//	TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
				//	centtex->SetNDC();
					centtex->Draw();
				//	TLatex *pttex = new TLatex(0.05,0.83,TrkPtBin_labels[ibin3]);
				//	pttex->SetNDC();
					pttex->Draw();
				//	TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
				//	Ajtex->SetNDC();
					Ajtex->Draw();
				//	TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
				//	Jettex->SetNDC();
					Jettex->Draw();


				}
				if(ibin==3){
					TLatex *centtex2 = new TLatex(0.15,0.9,CBin_labels[ibin]);
					centtex2->SetNDC();
					centtex2->Draw();
					TLatex *pttex2 = new TLatex(0.15,0.83,TrkPtBin_labels[ibin3]);
					pttex2->SetNDC();
					pttex2->Draw();
					TLatex *Ajtex2 = new TLatex(0.80,0.90,AjBin_labels[0]);
					Ajtex2->SetNDC();
					Ajtex2->Draw();
					TLatex *Jettex2 = new TLatex(0.75, 0.85, JetLabel_labels[0]);
					Jettex2->SetNDC();
					Jettex2->Draw();


				}	

				cc26->SaveAs("CC26_S_Scaled_overMEoverM0_NotScaled_PerEtaPerPhi_MinusEtaSideBand_ProjectionY_AjUnder22.png");	
				cc26->SaveAs("PbPb_Leading_PhiProjections.pdf");	
			
			cc27->cd(4*(ibin3+1)-ibin);
			//gStyle->SetPadRightMargin (0.001);

			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


			hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


			//hJetTrackMELeading[ibin][ibin2][ibin3]->Draw("LEGO2");	
			//return 23;


			hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1./me00[ibin][ibin2][ibin3]);	
			//hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(100 /me00[ibin][ibin2][ibin3]);
		
	//	hJetTrackMELeading[ibin][ibin2][ibin3]->Draw("LEGO2");	
	//		return 23;
		//	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() ));
		//	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() ));
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() * PhiBinWidth * EtaBinWidth));	
				
			//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] -> Draw("LEGO2");
			//return 23;
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Divide(hJetTrackMELeading[ibin][ibin2][ibin3]);
			//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] -> Draw("LEGO2");
			//return 23;	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> CenterTitle();
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> CenterTitle();
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(1.6);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(1.6);
		      	//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->SetAxisRange(-3., 3.,"X");
		      	//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->SetAxisRange(-3., 3.,"X");

				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-3.,3.);	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.06);

	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Draw("LEGO2");	
				//return 24;
				// make a clone of the output
				result[ibin][ibin2][ibin3] = (TH2D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Clone((TString)("result_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])); 

				 llimiteta = result[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-3.+0.0001); 
				 rlimiteta = result[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-1.5 - 0.0001); 

				cout << "llimiteta = " << llimiteta << " rlimiteta = " << rlimiteta << endl;
//return 24;
				background_left[ibin][ibin2][ibin3] = (TH1D*)result[ibin][ibin2][ibin3]->ProjectionY((TString)("LeftSideband"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]), llimiteta, rlimiteta);

				llimiteta =  result[ibin][ibin2][ibin3]->GetXaxis()->FindBin(1.5 + 0.0001); 
				rlimiteta =  result[ibin][ibin2][ibin3]->GetXaxis()->FindBin(3.0 - 0.0001); 

				background_proj[ibin][ibin2][ibin3] = (TH1D*)result[ibin][ibin2][ibin3]->ProjectionY((TString)("Sideband"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]), llimiteta, rlimiteta);

				background_proj[ibin][ibin2][ibin3]->Add(	background_left[ibin][ibin2][ibin3] )	;
				background_proj[ibin][ibin2][ibin3]->Scale(1. / (2*(rlimiteta - llimiteta+1)))	;

//				background_proj[ibin][ibin2][ibin3]->Draw();

			background[ibin][ibin2][ibin3] = (TH2D*)result[ibin][ibin2][ibin3]->Clone((TString)("Background"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
			

		for (int k = 1; k < result[ibin][ibin2][ibin3]->GetNbinsY(); k++){
		temp1 = background_proj[ibin][ibin2][ibin3]->GetBinContent(k);	
		err1  = background_proj[ibin][ibin2][ibin3]->GetBinError(k);

		for (int m = 1; m < result[ibin][ibin2][ibin3]->GetNbinsX(); m++){
		//cout << "temp1 = " << temp1 << endl; 
		background[ibin][ibin2][ibin3]->SetBinContent(m, k, temp1);
		background[ibin][ibin2][ibin3]->SetBinError(m, k, err1);
		//cout << "I am working " << endl;
		}		
		}

		result[ibin][ibin2][ibin3]->Add(background[ibin][ibin2][ibin3], -1.);
		//background[ibin][ibin2][ibin3]->Draw("LEGO2");



        	//Signal[ibin][ibin2][ibin3] =    (TH1D*)result[ibin][ibin2][ibin3]->ProjectionY(Form("Phi_Signal%d%d%d",ibin,ibin2,ibin3),36, 65);	
			//	double MaxSignal[6][6];
			//	MaxSignal[ibin][ibin3] = Signal[ibin][ibin2][ibin3]->GetMaximum();
			//	Signal[ibin][ibin2][ibin3]->SetMaximum(MaxSignal[ibin][ibin3]+30);
			
				Signal[ibin][ibin2][ibin3]->SetMarkerStyle(20);
				Signal[ibin][ibin2][ibin3]->SetMarkerColor(kRed);
				Signal[ibin][ibin2][ibin3]->Rebin(4);
				Signal[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#phi");
				Signal[ibin][ibin2][ibin3]->Scale(EtaBinWidth/4.);
				Signal[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(0.8);
				Signal[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(0.8);
				Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				Signal[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(0.8);
				Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				Signal[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
				Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				Signal[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.05);
				Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetRangeUser(-1.5,1.5);
				
				//	Signal[ibin][ibin2][ibin3]->Draw();
		//From H	check_old_phi[g][i][j][l] = result2[g][i][j][l]->ProjectionY(summedchecknamePhi,llimiteta,rlimiteta);
				/*double*/ etalim = 1.5;
				/*double*/ dx_phi = 0.;
			  llimiteta = result[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-etalim+.001);
			  rlimiteta = result[ibin][ibin2][ibin3]->GetXaxis()->FindBin(etalim-.001);

			cout << "llimiteta = " << llimiteta << " rlimiteta = " << rlimiteta << endl; 

			          check_old_phi[ibin][ibin2][ibin3] =  (TH1D*)result[ibin][ibin2][ibin3]->ProjectionY(Form("Phi_Signal%d%d%d",ibin,ibin2,ibin3),llimiteta,rlimiteta);	
				 dx_phi = check_old_phi[ibin][ibin2][ibin3]->GetBinWidth(1);
				 // check_old_phi[ibin][ibin2][ibin3]->Scale(1./30);
				  check_old_phi[ibin][ibin2][ibin3]->Scale(dx_phi);
 
			//	 check_old_phi[ibin][ibin2][ibin3]->Draw();
			//From H	check_old_phi_rebin[g][i][j][l] =(TH1D*)Rebin_dPhi(check_old_phi[g][i][j][l]);
				check_old_phi_rebin[ibin][ibin2][ibin3] =(TH1D*)Rebin_dPhi(check_old_phi[ibin][ibin2][ibin3]);
				//	 dx_phi = check_old_phi_rebin[ibin][ibin2][ibin3]->GetBinWidth(1);
				//	  check_old_phi_rebin[ibin][ibin2][ibin3]->Scale(1./dx_phi);
		 		
				//TString newchecknamePhi = in_name;
	
				  check_new_phi[ibin][ibin2][ibin3] = result[ibin][ibin2][ibin3]->ProjectionY(Form("Phi_Signal2%d%d%d",ibin,ibin2,ibin3),llimiteta,rlimiteta);
				  check_new_phi[ibin][ibin2][ibin3]->Scale(1/dx_phi);

				
				check_old_phi_rebin[ibin][ibin2][ibin3]->SetMarkerStyle(20);
				check_old_phi_rebin[ibin][ibin2][ibin3]->SetMarkerColor(kRed);
				check_old_phi_rebin[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#phi");
				check_old_phi_rebin[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#frac{1}{N_{jet}}#frac{d^{2}N}{ d#Delta#phi dp_{T}}");
				check_old_phi_rebin[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(0.8);
				check_old_phi_rebin[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(0.8);
				check_old_phi_rebin[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				check_old_phi_rebin[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				check_old_phi_rebin[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(0.8);
				check_old_phi_rebin[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				check_old_phi_rebin[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);   			
				check_old_phi_rebin[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				check_old_phi_rebin[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.05);   	//	result[ibin][ibin2][ibin3]->Draw("LEGO2");
				check_old_phi_rebin[ibin][ibin2][ibin3]->GetXaxis()-> SetRangeUser(-1.5,1.5);	
				check_old_phi_rebin[ibin][ibin2][ibin3]->GetYaxis()-> SetRangeUser(-0.5, 7);
				//check_old_phi_rebin[ibin][ibin2][ibin3]->Rebin(3);
				//check_old_phi_rebin[ibin][ibin2][ibin3]->Scale(1/3.);
				//check_old_phi_rebin[ibin][ibin2][ibin3]->Draw();
			//	check_old_phi_rebin[ibin][ibin2][ibin3]->Draw("Y+");
			//	check_old_phi_rebin[ibin][ibin2][ibin3]->Draw("Y-", "SAME");
		
				check_old_phi_rebin[ibin][ibin2][ibin3]->Draw();
			
			
					if(ibin<3){
				//	TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
				//	centtex->SetNDC();
					centtex->Draw();
				//	TLatex *pttex = new TLatex(0.05,0.83,TrkPtBin_labels[ibin3]);
				//	pttex->SetNDC();
					pttex->Draw();
				//	TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
				//	Ajtex->SetNDC();
					Ajtex->Draw();
				//	TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
				//	Jettex->SetNDC();
					Jettex->Draw();


				}
				if(ibin==3){
					//TLatex *centtex2 = new TLatex(0.15,0.9,CBin_labels[ibin]);
					TLatex *centtex2 = new TLatex(0.20,0.90,CBin_labels[ibin]);
					centtex2->SetNDC();
					centtex2->Draw();
					//TLatex *pttex2 = new TLatex(0.15,0.83,TrkPtBin_labels[ibin3]);
					TLatex *pttex2 = new TLatex(0.20,0.83,TrkPtBin_labels[ibin3]);
					pttex2->SetNDC();
					pttex2->Draw();
					TLatex *Ajtex2 = new TLatex(0.85,0.90,AjBin_labels[0]);
					Ajtex2->SetNDC();
					Ajtex2->Draw();
					TLatex *Jettex2 = new TLatex(0.80, 0.85, JetLabel_labels[0]);
					Jettex2->SetNDC();
					Jettex2->Draw();


				}
	
		
				
				cc27->SaveAs("CC27_S_ScaledMEoverM0_NotScaled_PerEtaPerPhi_MinusEtaSideBand_ProjectionY_AjUnder22.png");
			

				cc28->cd(4*(ibin3+1)-ibin);
					
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


			hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


			//hJetTrackMELeading[ibin][ibin2][ibin3]->Draw("LEGO2");	
			//return 23;


			hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1./me00[ibin][ibin2][ibin3]);	
			//hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(100 /me00[ibin][ibin2][ibin3]);
		
	//	hJetTrackMELeading[ibin][ibin2][ibin3]->Draw("LEGO2");	
	//		return 23;
		//	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() ));
		//	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() ));
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() * PhiBinWidth * EtaBinWidth));	
				
			//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] -> Draw("LEGO2");
			//return 23;
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Divide(hJetTrackMELeading[ibin][ibin2][ibin3]);
			//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] -> Draw("LEGO2");
			//return 23;	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> CenterTitle();
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> CenterTitle();
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(1.6);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(1.6);
		      	//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->SetAxisRange(-3., 3.,"X");
		      	//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->SetAxisRange(-3., 3.,"X");

				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-3.,3.);	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.06);

	
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Draw("LEGO2");	
	//			return 28;
				// make a clone of the output
				result[ibin][ibin2][ibin3] = (TH2D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Clone((TString)("result_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])); 

				 llimiteta = result[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-3.+0.0001); 
				 rlimiteta = result[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-1.5 - 0.0001); 

				cout << "llimiteta = " << llimiteta << " rlimiteta = " << rlimiteta << endl;
//return 24;
				background_left[ibin][ibin2][ibin3] = (TH1D*)result[ibin][ibin2][ibin3]->ProjectionY((TString)("LeftSideband"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]), llimiteta, rlimiteta);

				llimiteta =  result[ibin][ibin2][ibin3]->GetXaxis()->FindBin(1.5 + 0.0001); 
				rlimiteta =  result[ibin][ibin2][ibin3]->GetXaxis()->FindBin(3.0 - 0.0001); 

				background_proj[ibin][ibin2][ibin3] = (TH1D*)result[ibin][ibin2][ibin3]->ProjectionY((TString)("Sideband"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]), llimiteta, rlimiteta);

				background_proj[ibin][ibin2][ibin3]->Add(	background_left[ibin][ibin2][ibin3] )	;
		//		background_proj[ibin][ibin2][ibin3]->Scale(1. / (2*(rlimiteta - llimiteta+1)))	;
				background_proj[ibin][ibin2][ibin3]->Scale(1. / (2*(rlimiteta - llimiteta)))	;

			//	background_proj[ibin][ibin2][ibin3]->Draw();  //good!
	
				background[ibin][ibin2][ibin3] = (TH2D*)result[ibin][ibin2][ibin3]->Clone((TString)("Background"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));	
//return 27;


				llimiteta = result[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-1.5 + 0.0001); 
         			rlimiteta = result[ibin][ibin2][ibin3]->GetXaxis()->FindBin( 1.5 - 0.0001); 


				Signal[ibin][ibin2][ibin3] =    (TH1D*)result[ibin][ibin2][ibin3]->ProjectionY(Form("Phi_Signal%d%d%d",ibin,ibin2,ibin3),llimiteta, rlimiteta);
	//			Signal[ibin][ibin2][ibin3] -> Scale(1. / (rlimiteta - llimiteta+1))	;
				Signal[ibin][ibin2][ibin3] -> Scale(1. / (rlimiteta - llimiteta))	;
		//		Signal[ibin][ibin2][ibin3] -> Draw();
				
				//Signal[ibin][ibin2][ibin3] -> Scale(1. / (30))	;

				Signal[ibin][ibin2][ibin3] -> Add(	background_proj[ibin][ibin2][ibin3], -1);

	
			//	double MaxSignal[6][6];
			//	MaxSignal[ibin][ibin3] = Signal[ibin][ibin2][ibin3]->GetMaximum();
			//	Signal[ibin][ibin2][ibin3]->SetMaximum(MaxSignal[ibin][ibin3]+30);
			
				Signal[ibin][ibin2][ibin3]->SetMarkerStyle(20);
				Signal[ibin][ibin2][ibin3]->SetMarkerColor(kRed);
				Signal[ibin][ibin2][ibin3]->Rebin(3);
				Signal[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#phi");
			
				//Signal[ibin][ibin2][ibin3]->Scale(EtaBinWidth/4.);
				Signal[ibin][ibin2][ibin3]->Scale(1./3.);
				Signal[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(0.8);
				Signal[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(0.8);
				Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				Signal[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(0.8);
				Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				Signal[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
				Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				Signal[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.05);
				Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetRangeUser(-1.5,1.5);
				
				Signal[ibin][ibin2][ibin3]->Draw();

					
				cc28->SaveAs("CC28_SubtractBeforeFilling_AjUnder22.png");

	
			}
		}
	}

	//me_proj_canvas->SaveAs((TString)("All_ME_Projections_Inclusive_"+datalabel+"_Leading_Under22.png"));
	//result_proj_canvas->SaveAs((TString)("Corrected_dEta_Incluisve_Projections_"+datalabel+"_Leading_Under22.png"));

	//cout << "Fitting the central part of eta" << endl;
	//ProEta->Fit("pol0","","",-0.2,0.2);
	//float s = pol0->GetParameter(0); 
	//cout << "The value we need to normalize for is = " << s << endl;

	//cout << "Normalizing the Mixed Events" << endl;
	//hists_hJetTrackMECent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Scale(EtaBins/s);

	//hists_hJetTrackSignalBackgroundCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Draw("LEGO2");
	//cout << "Cloning Mixed Event to devide by it" << endl;
	//TH2D* ME = hists_hJetTrackMECent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Clone("ME");

	//cout << "Divide Signal by ME ...." <<  endl;
	//hists_hJetTrackSignalBackgroundCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Divide(ME);
	//hists_hJetTrackSignalBackgroundCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Draw("Lego2");

	//hists_hJetTrackSignalBackgroundCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("ResultPhi",20, 80)->Draw();

	//So lets remember, we have 100 bins for eta. and the range of eta is |eta| < [-5, 5] (obtained from Proj Y -> it was created as such)
	//so 10 / 100 = .1 (eta/bin) or 10 bins per eta 
	//1.5 < dEta < 3 = 1.5 eta, and 1.5*10 = 15
	//but we cut on 3, do we have
	//so side band 1.5 < |Deta| < 3
	//and signal |Deta| < 1.5}
	//
	//It is important to know which bin "exactly" is within the interval and which bin is not. You can do the math:
	//[0,20] removed    ----->  dEta < -3 
	//[21,35] side-band       -3 < dEta < -1.5 <--------------------------------- T
	//[36-65] Signal           -1.5 < dEta < 1.5                                  |
	//[66-80] side-band        1.5 < dEta < 3                                     |
	//[81-100] removed         3 < dEta                                           |
	//                                                                            |
	//for a check, run this command:                                              | 
	// ProEta->GetBinLowEdge(21)                                                  |
	//which gives                                                                 |
	//(const Double_t)(-3.00000000000000000e+00) ---------------------------------|       

	fout_inc->cd();
	fout_inc->Write();
	cout << "i = " << i << " j = " << j << endl;
	return i;
	}                                           
