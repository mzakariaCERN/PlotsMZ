
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

	int macro_Plots_AjOver22_Leading(int i=0)
	{
	TFile  *finbg, *fout_inc ;

	//TFile *fin = TFile::Open("Data2011_trkPtCut1_All_AJ11.root");  // The input file
	//TFile *fin = TFile::Open("Data2011_All_AJ22.root");  // The input file
	//TFile *fin = TFile::Open("Data2011_All_AJOver22.root");  // The input file
	//TFile *fin = TFile::Open("Data2011_AjUnder22.root");  // The input file
	TFile *fin = TFile::Open("Data2011_AjUnder22_Only13.root");  // The input file
	//TFile *fin = TFile::Open("FromH/PbPb_Leading_Correlations.root");  // The input file
	//TFile *fin = TFile::Open("Data2011_All_AJ33.root");  // The input file
	//TFile *fin = TFile::Open("/home/mzakaria3/Documents/Research/MyGitProject/PlotsMZ/Data2011_trkPtCut1_All_AJ11.root");  // The input file
	//	TFile *fin = TFile::Open("/home/mzakaria3/Documents/Research/MyGitProject/PlotsMZ/Data2011_trkPtCut1_All_Aj11_20150318.root");  // The input file
	//TFile *fin = TFile::Open("/home/mzakaria3/Documents/Research/MyGitProject/PlotsMZ/Data2011_trkPtCut1_All_Aj11_20150319.root");  // The input file
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
	float left_level = 0, right_level = 0, left_level_error = 0, right_level_error = 0;
	//	int me00bin;

	float PtBins[nPtBins+1] = {100, 300};
	int NDijets[nCBins] = 0;
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


	//	TCanvas *c1 = new TCanvas("c1", "C1_SignalOverME", 0, 0, 1500, 1500);	
	//	TCanvas *c2 = new TCanvas("c2", "C2_SignalOverME_projPhi_LB", 0, 0, 1500, 1500);	
	//	TCanvas *c3 = new TCanvas("c3", "C3_SignalOVerME_projPhi_RB", 0, 0, 1500, 1500);	
	//	TCanvas *c4 = new TCanvas("c4", "C4_SignalOverME_projPhi_LB", 0, 0, 1500, 1500);	
	//	TCanvas *c5 = new TCanvas("c5", "C5_SignalOverME_projPhi_LB+RB", 0, 0, 1500, 1500);	
	//	TCanvas *c6 = new TCanvas("c6", "C6_SignalOverME_projPhi_Signal", 0, 0, 1500, 1500);	
	//	TCanvas *c7 = new TCanvas("c7", "C7_SignalOverME_projPhi_SignalNormalized", 0, 0, 1500, 1500);	
	//	TCanvas *c8 = new TCanvas("c8", "C8_TestSideBinninL", 0, 0, 1500, 1500);	
	//	TCanvas *c9 = new TCanvas("c9", "C9_TestSideBinninR", 0, 0, 1500, 1500);	
	//	TCanvas *c10 = new TCanvas("c10", "C10_SignalOverME_projPhi_BinFit", 0, 0, 1500, 1500);	


	TCanvas *cc1 = new TCanvas("cc1", "CC1_me_projeta", 0, 0, 1500, 1500);	
	TCanvas *cc2 = new TCanvas("cc2", "CC2_me_projphi", 0, 0, 1500, 1500);	
	TCanvas *cc3 = new TCanvas("cc3", "CC3_SignalNoScale", 0, 0, 1500, 1500);	
	TCanvas *cc4 = new TCanvas("cc4", "CC4_S", 0, 0, 1500, 1500);	
	TCanvas *cc5 = new TCanvas("cc5", "CC5_ME", 0, 0, 1500, 1500);	
	TCanvas *cc6 = new TCanvas("cc6", "CC6_MEnormalized", 0, 0, 1500, 1500);	
	TCanvas *cc7 = new TCanvas("cc7", "CC6j_LBoverRB", 0, 0, 1500, 1500);	
	TCanvas *cc8 = new TCanvas("cc8", "CC7_EtaBandSymmetry", 0, 0, 1500, 1500);	
	TCanvas *cc9 = new TCanvas("cc9", "CC8_EtaBandSymmetry_HT", 0, 0, 1500, 1500);	
	TCanvas *cc10 = new TCanvas("cc10", "CC9_EtaBandSymmetry_HT", 0, 0, 1500, 1500);	
	TCanvas *cc11 = new TCanvas("cc11", "CC10_EtaBandSymmetry_HT", 0, 0, 1500, 1500);	
	TCanvas *cc12 = new TCanvas("cc12", "CC12_CC6j_LBoverRB_magnified", 0, 0, 1500, 1500);	
	TCanvas *cc13 = new TCanvas("cc13", "CC13_SignalLBoverRB_magnified", 0, 0, 1500, 1500);	


	//	c1->Divide(4,4,0.0,0.0);
	//	c6->Divide(4,4,0.0,0.0);
	//	c8->Divide(4,4,0.0,0.0);
	//	c9->Divide(4,4,0.0,0.0);
	//	c10->Divide(4,4,0.0,0.0);

	cc1->Divide(4,4,0.0001,0.0001);   // 4 Raw by 4 Columns and the separation is 0.1 and 0.1
	cc2->Divide(4,4,0.0001,0.0001);
	cc3->Divide(4,4,0.0001,0.0001);
	cc4->Divide(4,4,0.0001,0.0001);
	cc5->Divide(4,4,0.0001,0.0001);
	cc6->Divide(4,4,0.0001,0.0001);
	cc7->Divide(4,4,0.0001,0.0001);
	cc8->Divide(4,4,0.0001,0.0001);
	cc9->Divide(4,4,0.0001,0.0001);
	cc10->Divide(4,4,0.0001,0.0001);
	cc11->Divide(4,4,0.0001,0.0001);

	cc11->Divide(4,4,0.000,0.000);
	cc12->Divide(4,4,0.0001,0.0001);
	cc13->Divide(4,4,0.0001,0.0001);


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

	//if(!Logical)
	//  printf("Not logical value at line number %d in file %s\n", __LINE__, __FILE__);  // used to print the line (debugging)

	TString  datalabel;

	gStyle->SetOptStat(0);   
	gStyle->SetPadBottomMargin(0.006);
	gStyle->SetPadTopMargin   (0.006);
	gStyle->SetPadLeftMargin  (0.006);
	gStyle->SetPadRightMargin (0.006);
	gStyle->SetTitleSize(0.06,"xyz");
	gStyle->SetLabelSize(0.06);
				//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(0.8);
				//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(0.8);
				//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);


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
	double A[nCBins][nPtBins][nTrkPtBins] = 0;
	double B[nCBins][nPtBins][nTrkPtBins] = 0;
	double C[nCBins][nPtBins][nTrkPtBins] = 0;

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
				
			//	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

				hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
				
				//hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

				Correlation[ibin][ibin2][ibin3] = (TH2D*)  fin->Get((TString) (desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

						//	hJetTrackSignalBackground[ibin][ibin2][ibin3]->Draw("LEGO@");
			//	cout << "Plot of the Signal Events" << endl;
			//	hJetTrackMELeading[ibin][ibin2][ibin3]->Draw("LEGO@");
			//	cout << "Plot of the Mixed Events" << endl;


			//	hJetTrackME[ibin][ibin2][ibin3]->ProjectionX("ProEta")->Draw();
			//	cout << "Eta Projection" <<  endl;

			//	hJetTrackME[ibin][ibin2][ibin3]->ProjectionY("me_projphi")->GetBinWidth(1);
			//	hJetTrackME[ibin][ibin2][ibin3]->ProjectionY("me_projphi");
			//	hJetTrackME[ibin][ibin2][ibin3]->ProjectionY("me_projphi");
			//	cout << "Phi Projection" << endl;
			//	>ProjectionY("me_projphi")->Draw();
			//	cout << "dPhi Bin Width = " << hJetTrackME[ibin][ibin2][ibin3]   << endl; 
				//cout << "Number of Eta bins:" << ProEta->GetNbinsX() << endl;
				//float EtaBins = ProEta->GetNbinsX() ;
		//		cout << "Number of Phi bins:" << me_projphi->GetNbinsX() << endl;


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

				fout_inc->cd();
				//cout << "Canvas number = " << 4*(ibin3+1)-ibin << endl;
				//me_proj_canvas->cd(4*(ibin3+1)-ibin);
				
				yield_inc[ibin][ibin2][ibin3] =  (TH2D*) hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Clone((TString)("Yield_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

			   yield_inc_copy[ibin][ibin2][ibin3] =  (TH2D*) hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Clone((TString)("Yield_copy"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

				
				//ProEta->Fit("pol0","","",-0.2,0.2);
				//float s = pol0->GetParameter(0);  // need to make this an array
				me_projeta[ibin][ibin2][ibin3] =  hJetTrackMELeading[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_eta%d%d%d",ibin,ibin2,ibin3),1,100); //what is the 1,100 for?
				me_projphi[ibin][ibin2][ibin3] =  hJetTrackMELeading[ibin][ibin2][ibin3]->ProjectionY(Form("me_proj_phi%d%d%d",ibin,ibin2,ibin3),1,100); //what is the 1,100 for?
				//me_projeta[ibin][ibin2][ibin3]->Scale (1./100);
				EtaProj[ibin][ibin2][ibin3] = (TH1D*)Correlation[ibin][ibin2][ibin3]->ProjectionX((TString)("Eta_Proj_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),lbin,rbin);
	  
				EtaProjDummy[ibin][ibin2][ibin3] = (TH1D*) EtaProj[ibin][ibin2][ibin3]->Clone((TString)("Eta_Proj_Dummy_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

			

				cout << "Trying to get Number of eta bins = " << endl;
				cout << me_projeta[ibin][ibin2][ibin3]->GetNbinsX() << endl;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				cc1->cd(4*(ibin3+1)-ibin);
				// cc1->SetTickx(1);
				//cc1->SetGridy();
				//cc1->SetTopMargin(0.001);
				//cc1->SetBottomMargin(0.001);
				//cc1->Modified();
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
					TLatex *pttex = new TLatex(0.15,0.80,TrkPtBin_labels[ibin3]);
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
			//		TLatex *pttex = new TLatex(0.15,0.80,TrkPtBin_labels[ibin3]);
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
	cc1->SaveAs("cc1_ME_EtaProj_Leading_Over22.png");
	//##//
	//##//			//	gPad->Update();
	//##//				
	//##//			//	return 0;
	//##/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//##//				
	cc2->cd(4*(ibin3+1) - ibin);
	me_projeta[ibin][ibin2][ibin3] =  hJetTrackMELeading[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_eta%d%d",ibin,ibin2,ibin3),1,100); //what is the 1,100 for?
					
	//##////			me_projphi[ibin][ibin2][ibin3]->DrawCopy();
	//##//				//me_projphi[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("Delta Phi");
	//##//				
			//	PhiBinWidth =  me_projphi[ibin][ibin2][ibin3]->GetBinWidth(1) ;
	//##//				cout << "Phi Projection = " << me_projphi[ibin][ibin2][ibin3]->GetBinWidth(1) << endl ;
				//me00[ibin][ibin2][ibin3] = me_projeta[ibin][ibin2][ibin3]->GetFunction("pol0")->GetParameter(0);
			//	eta_level = pol0->GetParameter(0);
	//##//				//me_proj_canvas->cd(4*(ibin3+1)-ibin);				
	//##//				//me_projphi[ibin][ibin2][ibin3]->Draw();
	//##//				//me_proj_canvas->Update();
	//##//
	//##//
	//##//
	//##///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////				
	//##//
	//##//				cout << "Eta Flat Fitting Parameter value = " << me00[ibin][ibin2][ibin3] << endl;
				//me_projeta[ibin][ibin2][ibin3] -> Fit("pol0","0","",-.2,.2);
				//me00[ibin][ibin2][ibin3] = me_projeta[ibin][ibin2][ibin3]->GetFunction("pol0")->GetParameter(0);	
	me_projeta[ibin][ibin2][ibin3]->Scale (100./me00[ibin][ibin2][ibin3]);
				//me_projeta[ibin][ibin2][ibin3]->Scale (10./me00[ibin][ibin2][ibin3]);
	//##//				me_projeta[ibin][ibin2][ibin3]->Scale (1./me00[ibin][ibin2][ibin3]);
					hJetTrackMELeading[ibin][ibin2][ibin3]->Scale (100./me00[ibin][ibin2][ibin3]);
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
	//##///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//##//				
	//##//				yield_inc[ibin][ibin2][ibin3]->Divide(hJetTrackMELeading[ibin][ibin2][ibin3]);
	//##//				//yield_inc[ibin][ibin2][ibin3]->Draw();
	//##//				//yield_inc_proj[ibin][ibin2][ibin3] = (TH1D*) yield_inc[ibin][ibin2][ibin3]->ProjectionX((Form("ResultEta%d%d%d",ibin,ibin2,ibin3)),20, 80);
	//##//				yield_inc_proj[ibin][ibin2][ibin3] = (TH1D*) yield_inc[ibin][ibin2][ibin3]->ProjectionX((Form("ResultEta%d%d%d",ibin,ibin2,ibin3)), 1, 50);
	//##//				//yield_inc_proj[ibin][ibin2][ibin3] = (TH1D*) yield_inc[ibin][ibin2][ibin3]->ProjectionX((Form("ResultEta%d%d%d",ibin,ibin2,ibin3)));
	//##//				yield_inc_proj[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");
	//##//
	//##//				yield_inc_proj[ibin][ibin2][ibin3]->Draw();
	//##//				
	//##//				//yield_inc_proj_copy[ibin][ibin2][ibin3] = (TH1D*) yield_inc[ibin][ibin2][ibin3]->ProjectionX((Form("ResultEta_copy%d%d%d",ibin,ibin2,ibin3)),20, 80);
	//##//				//yield_inc_proj_copy[ibin][ibin2][ibin3] = (TH1D*) yield_inc[ibin][ibin2][ibin3]->ProjectionX((Form("ResultEta_copy%d%d%d",ibin,ibin2,ibin3)), 1, 50); // this range is to get
	//##//				yield_inc_proj_copy[ibin][ibin2][ibin3] = (TH1D*) yield_inc[ibin][ibin2][ibin3]->ProjectionX((Form("ResultEta_copy%d%d%d",ibin,ibin2,ibin3))); // this range is to get
	//##//				// the jet and not the ridge
	//##//					yield_inc_proj_copy[ibin][ibin2][ibin3]->SetMarkerStyle(20);
	//##//						yield_inc_proj_copy[ibin][ibin2][ibin3]->SetMarkerColor(kRed);
	//##//				//yield_inc_proj_copy[ibin][ibin2][ibin3] -> Draw("same");
	//##//				//return 0;	
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
					TLatex *Jettex = new TLatex(0.715, 0.80, JetLabel_labels[0]);
					Jettex->SetNDC();
					Jettex->Draw();

					}

					cc2->SaveAs("cc2_MEOverME00_Leading_Over22.png");

	//##//
	//##//							
				cc3->cd(4*(ibin3+1)-ibin);

	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> CenterTitle();
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> CenterTitle();
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(1.6);
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(1.6);
	//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(0.8);
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.06);

	//			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]-> SetTitleSize(.1, "X");
	//			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis() -> SetLabelSize(0.2); // for the number not title
		hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->SetAxisRange(-3., 3.,"X");
		hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-3.,3.);				
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Draw("LEGO2");
	//##//				yield_inc[ibin][ibin2][ibin3]->Draw("Lego2");
	//return 3;
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
					TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
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
					TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
					Jettex->SetNDC();
					Jettex->Draw();


				}
	cc3->SaveAs("cc3_SignalNoScale_Leading_Over22.png");
//return 3;
				cc4->cd(4*(ibin3+1)-ibin);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

				PhiBinWidth =  me_projphi[ibin][ibin2][ibin3]->GetBinWidth(1) ;
				EtaBinWidth =  hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->GetBinWidth(1) ;
				PhiBinWidth =  hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->GetBinWidth(1) ;
				cout << "Eta Bint Width = " << 	 hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->GetBinWidth(1) << endl;


				//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / (3.0 *  Aj[ibin]->GetEntries() * PhiBinWidth));
				//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / (0.1 *  Aj[ibin]->GetEntries() * PhiBinWidth));
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / (EtaBinWidth *  Aj[ibin]->GetEntries() * PhiBinWidth));
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
					TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
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
					TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
					Jettex->SetNDC();
					Jettex->Draw();


				}
			
			cc4->SaveAs("cc4_S_Leading_Over22.png");
//return 4;
			cc5->cd(4*(ibin3+1)-ibin);
		
			//	hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
			

				//hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1.0 / (3.0 *  Aj[ibin]->GetEntries() * PhiBinWidth));
			//	hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1.0 / (EtaBinWidth *  Aj[ibin]->GetEntries() * PhiBinWidth));
			//	hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1.0/me00[ibin][ibin2][ibin3]);
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis() -> CenterTitle();
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetYaxis() -> CenterTitle();
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(1.6);
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(1.6);
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.06);

			//	hJetTrackMELeading[ibin][ibin2][ibin3]->SetAxisRange(-3., 3.,"X");
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-3.,3.);	
			//	hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(100./me00[ibin][ibin2][ibin3]);
				
				hJetTrackMELeading[ibin][ibin2][ibin3]->Draw("LEGO2");

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
					TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
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
					TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
					Jettex->SetNDC();
					Jettex->Draw();


				}


					cc5->SaveAs("cc5_MEOverME00_Notscaled_Leading_Over22.png");
//return 5;

	cc6->cd(4*(ibin3+1)-ibin);
		
			//	hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
			

				hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1.0 / (0.1 *  Aj[ibin]->GetEntries() * PhiBinWidth));
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis() -> CenterTitle();
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetYaxis() -> CenterTitle();
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(1.6);
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(1.6);
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.06);

				hJetTrackMELeading[ibin][ibin2][ibin3]->SetAxisRange(-3., 3.,"X");
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-3.,3.);				
				
				//hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1./me_projeta[ibin][ibin2][ibin3]);				
				//hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(100./me_projeta[ibin][ibin2][ibin3]);				
			//	hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(100./me00[ibin][ibin2][ibin3]);

				hJetTrackMELeading[ibin][ibin2][ibin3]->Draw("LEGO2");

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
					TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
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
					TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
					Jettex->SetNDC();
					Jettex->Draw();


				}



			cc6->SaveAs("cc6_MEOverME00_Scaled_Leading_Over22.png");

//return 6;

		cc7->cd(4*(ibin3+1)-ibin);
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


			hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

			hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(100./me00[ibin][ibin2][ibin3]);	
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
					TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
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
					TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
					Jettex->SetNDC();
					Jettex->Draw();


				}


				cc7->SaveAs("cc7_S_Scaled_OverMEOverME00_NotScaled_Leading_Over22.png");
//	return 7;
				cc8->cd(4*(ibin3+1)-ibin);

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
							double level = (left_level+right_level)/2.;
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
				cc8->SaveAs("cc8_RightBandLeftBandComparison_Leading_Over22.png");
								
//return 8;
				cc9->cd((4*(ibin3+1)-ibin));
			//	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

				
				
			//	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Divide(hJetTrackMELeading[ibin][ibin2][ibin3]);

				 //hists_hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Rebin2D(5,2);
				

				//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Draw("LEGO2");

				//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionX("Check%d%d%d",1, 50)->Draw();
				//yield_inc_proj[ibin][ibin2][ibin3] = (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionX((Form("ResultEta%d%d%d",ibin,ibin2,ibin3)),1, 50);

				//TH1D *Signal = (TH1D*)hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("Signal",36, 65);
				Signal[ibin][ibin2][ibin3] =    (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionY(Form("Phi_Signal%d%d%d",ibin,ibin2,ibin3),36, 65);	
				Signal[ibin][ibin2][ibin3]->SetMarkerStyle(20);
				Signal[ibin][ibin2][ibin3]->SetMarkerColor(kRed);
				Signal[ibin][ibin2][ibin3]->Rebin(2);
				Signal[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#phi");
				Signal[ibin][ibin2][ibin3]->Scale(EtaBinWidth/2.);
				Signal[ibin][ibin2][ibin3]->GetYaxis() -> SetTitleOffset(0.8);
				Signal[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(0.8);
				Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				Signal[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(0.8);
				Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				Signal[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
				Signal[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				Signal[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.05);

	
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
				//yield_inc_proj[ibin][ibin2][ibin3]->Draw();	


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
			cc9->SaveAs("cc9_PhiProjectionofSignalandBands_Leading_Over22.png")			;			
//return 9;
			cc10->cd(4*(ibin3+1)-ibin);
			SignalMinusRL[ibin][ibin2][ibin3] = (TH1D*)Signal[ibin][ibin2][ibin3]->Clone((TString)("signal_MinusRL"+CBin_strs[ibin]+"_"+TrkPtBin_strs[ibin3]));
			SignalMinusRL[ibin][ibin2][ibin3]->Add(	PhiProj_R[ibin][ibin2][ibin3], -1);
			SignalMinusRL[ibin][ibin2][ibin3]->SetAxisRange(-1.5, 1.5,"X");
		//	SignalMinusRL[ibin][ibin2][ibin3]->Scale(1/(3.0 *  Aj[ibin]->GetEntries() * PhiBinWidth));

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
			cc10->SaveAs("cc10_resultOnPhi_Leading_Over22.png");


				cc11->cd(4*(ibin3+1)-ibin);
			//SignalMinusRL[ibin][ibin2][ibin3] = (TH1D*)Signal[ibin][ibin2][ibin3]->Clone((TString)("signal_MinusRL"+CBin_strs[ibin]+"_"+TrkPtBin_strs[ibin3]));
			//SignalMinusRL[ibin][ibin2][ibin3]->Add(	PhiProj_R[ibin][ibin2][ibin3], -1);
			//SignalMinusRL[ibin][ibin2][ibin3]->SetAxisRange(-1.5, 1.5,"X");
			SignalMinusRL[ibin][ibin2][ibin3]->Scale(1/(Aj[ibin]->GetEntries()));

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
			cc11->SaveAs("cc11_resultOnPhi_Leading_Over22.png");	
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
							double level = (left_level+right_level)/2.;
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

			cc12->SaveAs("cc12_RightBandLeftBandComparison_mmagnified_Leading_Over22.png");
			
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
							double level = (left_level+right_level)/2.;
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

				cc13->SaveAs("cc13_RightBandLeftBandComparison_mmagnified_Leading_Over22.png");



	//##//				//yield_inc[ibin][ibin2][ibin3]->Scale(me00);
	//##//				//me_proj_canvas->Update();
	//##//			
	//##////				yield_inc[ibin][ibin2][ibin3]->Draw("LEGO2");    //This is the yield (corrected for ME) that we need.
	//##//
	//##//
	//##////return 0;
	//##//				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Write();
	//##//				yield_inc[ibin][ibin2][ibin3]->Write();
	//##//				yield_inc_copy[ibin][ibin2][ibin3]->Write();
	//##//				hJetTrackMELeading[ibin][ibin2][ibin3]->Write();
	//##//
	//##//
	//##//				//cout << "The value we need to normalize for is = " << s << endl;
	//##//
	//##//				//cout << "Normalizing the Mixed Events" << endl;
	//##//				//hJetTrackME[ibin][ibin2][ibin3]->Scale(EtaBins/s);
	//##//				//	hJetTrackME[ibin][ibin2][ibin3]->Scale(1./me00);
	//##//
	//##//				//TH2D* ME = hJetTrackME[ibin][ibin2][ibin3]->Clone("ME"); //make is an array later
	//##//
	//##//				//hJetTrackSignalBackground[ibin][ibin2][ibin3]->Divide(ME); //done below
	//##//
	//##//				//hJetTrackSignalBackground[ibin][ibin2][ibin3]->Draw("Lego2");
	//##//				///hJetTrackSignalBackground[ibin][ibin2][ibin3]->ProjectionY("ResultPhi",20, 80)->Draw();
	//##//				//me_projeta[ibin][ibin2][ibin3] = hJetTrackME[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_temp_%d%d%d",ibin,ibin2,ibin3),1,100);
	//##//				
	//##//				cc4->cd(4*(ibin3+1)-ibin);
	//##//				//yield_inc[ibin][ibin2][ibin3]->ProjectionY("ResultPhi",20, 80)->Draw();
	//##//				//yield_inc_proj[ibin][ibin2][ibin3]->Draw();
	//##//				EtaProjDummy[ibin][ibin2][ibin3] = (TH1D*) yield_inc_proj[ibin][ibin2][ibin3]->Clone((TString)("Eta_Proj_Dummy_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));
	//##////					EtaProjDummy[ibin][ibin2][ibin3]->Rebin(4);				
	//##////				EtaProjDummy[ibin][ibin2][ibin3]->Scale(1.0 / (3.0 *  Aj[ibin]->GetEntries() * PhiBinWidth));
	//##//	
	//##//			//		EtaProjDummy[ibin][ibin2][ibin3]->Draw();
	//##////					EtaProjDummy[ibin][ibin2][ibin3]->Fit("flat_fit","0","",-3.,-1.5);
	//##////			left_level = flat_fit->GetParameter(0);
	//##////					EtaProjDummy[ibin][ibin2][ibin3]->Fit("flat_fit","0","",1.5,3.);
	//##//	  
	//##////	  			right_level =  flat_fit->GetParameter(0);
	//##////				bg_fit[ibin][ibin2][ibin3] = new TF1("bg_fit%d%d%d","[0]+x-x",-3.,3.);
	//##////	  		double level = (left_level+right_level)/2.;
	//##////	  bg_fit[ibin][ibin2][ibin3]->SetParameter(0,(left_level+right_level)/2.);
	//##////					//EtaProjDummy[ibin][ibin2][ibin3]->SetMinimum((left_level+right_level)/2.)*0.95);
	//##////					EtaProjDummy[ibin][ibin2][ibin3]->SetAxisRange(-3, 3, "X");
	//##////					EtaProjDummy[ibin][ibin2][ibin3]->SetMinimum(level*0.9);
	//##////					EtaProjDummy[ibin][ibin2][ibin3]->SetMaximum(level*1.1);
	//##//					EtaProjDummy[ibin][ibin2][ibin3]->Draw();
	//##//
	//##//			// yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(10.);
	//##//			// yield_inc_proj[ibin][ibin2][ibin3]->SetAxisRange(0, 100, "Y");
	//##//		//	yield_inc_proj[ibin][ibin2][ibin3]->Draw();
	//##//					
	//##////					bg_fit[ibin][ibin2][ibin3]->Draw("same");
	//##//	//			yield_inc[ibin][ibin2][ibin3]->ProjectionX("ResultEta%d%d",20, 80)->DrawCopy();
	//##//	//			yield_inc[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("Delta eta");
	//##//
	//##//	//			//yield_inc[ibin][ibin2][ibin3]->GetXaxis()->SetRange(20, 80);					
	//##//	//			yield_inc[ibin][ibin2][ibin3]->GetXaxis()->SetRange(20, 80);					
	//##//			
	//##//									if(ibin<3){
	//##//					TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
	//##//					centtex->SetNDC();
	//##//					centtex->Draw();
	//##//					TLatex *pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
	//##//					pttex->SetNDC();
	//##//					pttex->Draw();
	//##//					TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
	//##//					Ajtex->SetNDC();
	//##//					Ajtex->Draw();
	//##//					TLatex *Jettex = new TLatex(0.70, 0.85, JetLabel_labels[0]);
	//##//					Jettex->SetNDC();
	//##//					Jettex->Draw();
	//##//
	//##//
	//##//				}
	//##//				if(ibin==3){
	//##//					TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	//##//					centtex->SetNDC();
	//##//					centtex->Draw();
	//##//					TLatex *pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	//##//					pttex->SetNDC();
	//##//					pttex->Draw();
	//##//					TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
	//##//					Ajtex->SetNDC();
	//##//					Ajtex->Draw();
	//##//					TLatex *Jettex = new TLatex(0.70, 0.85, JetLabel_labels[0]);
	//##//					Jettex->SetNDC();
	//##//					Jettex->Draw();
	//##//
	//##//
	//##//				}
	//##//
	//##//				
	//##//			//	return 1;
	//##//				//yield_inc_SBL[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionX(Form("ResultSBL%d",ibin,ibin2,ibin3),21,35 );
	//##//				
	//##//				cc5->cd(4*(ibin3+1)-ibin);
	//##//				//PhiProj_L[ibin][ibin2][ibin3]->Add(PhiProj_R[ibin][ibin2][ibin3]);
	//##//				//PhiProj_L[ibin][ibin2][ibin3]->Add(PhiProj_R[ibin][ibin2][ibin3]);
	//##//			        //PhiProj_L[ibin][ibin2][ibin3]->Draw();	
	//##//				
	//##//				//Signal[ibin][ibin2][ibin3]->Draw();
	//##//				//PhiProj_L[ibin][ibin2][ibin3]->Draw("SAME");
	//##//				PhiProj_L[ibin][ibin2][ibin3] =  (TH1D*) yield_inc[ibin][ibin2][ibin3]->ProjectionY(Form("ResultPhi_L%d%d%d",ibin,ibin2,ibin3),21, 35);
	//##//		//		PhiProj_L[ibin][ibin2][ibin3] =  (TH1D*) yield_inc_copy[ibin][ibin2][ibin3]->ProjectionY(Form("ResultPhi_L%d%d%d",ibin,ibin2,ibin3),36, 65);
	//##//		//		PhiProj_L[ibin][ibin2][ibin3]->Draw();
	//##//				//PhiProj_L[ibin][ibin2][ibin3]->GetXaxis()->SetRange(0,15 );		
	//##//				
	//##//				PhiProj_R[ibin][ibin2][ibin3] =  (TH1D*) yield_inc[ibin][ibin2][ibin3]->ProjectionY(Form("ResultPhi_R%d%d%d",ibin,ibin2,ibin3),66, 80);
	//##//		//		PhiProj_R[ibin][ibin2][ibin3]->Draw();
	//##//				//PhiProj_R[ibin][ibin2][ibin3]->GetXaxis()->SetRange(46,60);
	//##//				//Signal[ibin][ibin2][ibin3] = yield_inc_copy[ibin][ibin2][ibin3]->ProjectionY(Form("Signal%d%d%d",ibin,ibin2,ibin3),36, 65);
	//##//				Signal[ibin][ibin2][ibin3] =  (TH1D*) yield_inc[ibin][ibin2][ibin3]->ProjectionY(Form("Signal%d%d%d",ibin,ibin2,ibin3),36, 65);
	//##//				//Signal[ibin][ibin2][ibin3]->GetXaxis()->SetRange(16,45);
	//##//
	//##//				PhiProj_L[ibin][ibin2][ibin3]->Add(PhiProj_R[ibin][ibin2][ibin3]);
	//##//				PhiProj_L[ibin][ibin2][ibin3]->Scale(1.0 / (3.0 *  Aj[ibin]->GetEntries() * PhiBinWidth));
	//##//				PhiProj_L[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#phi");
	//##//				PhiProj_L[ibin][ibin2][ibin3]->SetMarkerStyle(20);
	//##//				PhiProj_L[ibin][ibin2][ibin3]->SetMarkerColor(kRed);
	//##//				PhiProj_L[ibin][ibin2][ibin3]->Rebin(4);
	//##//				PhiProj_L[ibin][ibin2][ibin3]->Draw();
	//##//
	//##//
	//##//				Signal[ibin][ibin2][ibin3]->Rebin(4);
	//##//				Signal[ibin][ibin2][ibin3]->Scale(1.0 / (3.0 *  Aj[ibin]->GetEntries() * PhiBinWidth));
	//##//				Signal[ibin][ibin2][ibin3]->SetMarkerColor(kBlack);
	//##//				Signal[ibin][ibin2][ibin3]->SetMarkerStyle(20);
	//##//			//	Signal[ibin][ibin2][ibin3]->SetMinimum(0);
	//##//			//	Signal[ibin][ibin2][ibin3]->SetMaximum(20000);
	//##//				Signal[ibin][ibin2][ibin3]->Draw("same");
	//##//			//	return 0;
	//##//							//	PhiProj_R[ibin][ibin2][ibin3]->Draw("SAME");
	//##//				//PhiProj_L[ibin][ibin2][ibin3]->Draw();
	//##//
	//##//	
	//##//				//yield_inc_proj[ibin][ibin2][ibin3]->Draw();
	//##//				//yield_inc_proj[ibin][ibin2][ibin3]->GetYaxis()->SetRangeUser(bg_fit[ibin][ibin2][ibin3]*0.9, bg_fit[ibin][ibin2][ibin3]*1.1);
	//##//				//yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(bg_fit[ibin][ibin2][ibin3]*1.1);
	//##//				//yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(bg_fit[ibin][ibin2][ibin3]*0.9);
	//##//				//cc5->Update();
	//##//				//yield_inc_proj[ibin][ibin2][ibin3]->Draw();
	//##//
	//##//						
	//##//				// TH1F *hint1 = new TH1F("hint1","h1 bins integral",100,-3,3);
	//##//				//yield_inc_proj[ibin][ibin2][ibin3]->SetAxisRange(-3, 3,"X");	
	//##//				//yield_inc_proj[ibin][ibin2][ibin3]->Draw();
	//##//				//return 5;	
	//##//			        //yield_inc_SBL[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionX("ResultSBL%d",21,80,"ed" );  // setting the range doesn't make much difference
	//##//			        //yield_inc_SBL[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionX( );  // setting the range doesn't make much difference
	//##//	//			yield_inc[ibin][ibin2][ibin3]->Draw("lego2");
	//##//				
	//##//	//			yield_inc[ibin][ibin2][ibin3]->ProjectionX("ResultEta%d%d",0, 100)->DrawCopy();
	//##//	//			yield_inc[ibin][ibin2][ibin3]->GetXaxis()->SetRange(20, 80);					
	//##//	//if(ibin<3){
	//##//	//				TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
	//##//	//				centtex->SetNDC();
	//##//	//				centtex->Draw();
	//##//	//				TLatex *pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
	//##//	//				pttex->SetNDC();
	//##//	//				pttex->Draw();
	//##//	//			}
	//##//	//			if(ibin==3){
	//##//	//				TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	//##//	//				centtex->SetNDC();
	//##//	//				centtex->Draw();
	//##//	//				TLatex *pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	//##//	//				pttex->SetNDC();
	//##//	//				pttex->Draw();
	//##//	//			}
	//##//
	//##//				//yield_inc[ibin][ibin2][ibin3]->ProjectionX("ResultEta%d",0, 100)->DrawCopy();	
	//##//				//yield_inc[ibin][ibin2][ibin3]->ProjectionX("ResultEta%d",0, 100)->DrawCopy();	
	//##//				//yield_inc[ibin][ibin2][ibin3]->ProjectionX("ResultEta%d",30, 70)->DrawCopy();
	//##//			//	ResultEta%d->Draw();	
	//##//			//	yield_inc[ibin][ibin2][ibin3]->GetXaxis()->SetRange(20, 80);					
	//##//				
	//##//	//			yield_inc_SBL[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionX("ResultSBL%d",21,80,"ed" );  // setting the range doesn't make much difference
	//##//	//			yield_inc_SBL[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_temp_%d%d%d",ibin,ibin2,ibin3),36,65);
	//##//				//yield_inc_SBL[ibin][ibin2][ibin3] = yield_inc_proj[ibin][ibin2][ibin3];			
	//##//			//	//yield_inc_SBL[ibin][ibin2][ibin3]->GetXaxis()->SetRange(21,35 );			
	//##//	//			//yield_inc_SBL[ibin][ibin2][ibin3]->SetAxisRange(-5.0, 1.5,"X");	
	//##//				//yield_inc_SBL[ibin][ibin2][ibin3]->SetMarkerStyle(20);
	//##//				//yield_inc_SBL[ibin][ibin2][ibin3]->SetMarkerColor(kRed);
	//##//				//yield_inc_SBL[ibin][ibin2][ibin3]->SetLineColor(kRed);
	//##//			//	//yield_inc_SBL[ibin][ibin2][ibin3]	=new TH1D(histname,"",100,5,5);				   
	//##//			//	//yield_inc_SBL[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_temp_%d%d%d",ibin,ibin2,ibin3),36,65);	
	//##//	//			//yield_inc_SBL[ibin][ibin2][ibin3]->GetXaxis()->SetRange(0,15 );
	//##//	//			//yield_inc_SBL[ibin][ibin2][ibin3]->Reset();
	//##//	//			//yield_inc_proj[ibin][ibin2][ibin3]->SetBins(15, -3.0, -1.5);
	//##//				//yield_inc_SBL[ibin][ibin2][ibin3]->Draw("same");
	//##//
	//##//
	//##//				//yield_inc_SBR[ibin][ibin2][ibin3] = yield_inc_proj[ibin][ibin2][ibin3];	
	//##//				////yield_inc_SBR[ibin][ibin2][ibin3]-> GetXaxis()->SetRange(65,80) ; // setting the range doesn't make much difference
	//##//	//			////yield_inc_SBR[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionX("ResultSBR%d",21,80,"ed" );  // setting the range doesn't make much difference
	//##//	//			//yield_inc_SBR[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionX("ResultSBR%d",21,80,"ed" );  // setting the range doesn't make much difference
	//##//	//			//yield_inc_SBR[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_temp_%d%d%d",ibin,ibin2,ibin3),36,65);
	//##//				//yield_inc_SBL[ibin][ibin2][ibin3]->SetMarkerStyle(20);
	//##//				//yield_inc_SBL[ibin][ibin2][ibin3]->SetMarkerColor(kBlue);
	//##//				//yield_inc_SBL[ibin][ibin2][ibin3]->SetLineColor(kBlue);
	//##//	//			//yield_inc_SBL[ibin][ibin2][ibin3]->GetXaxis()->SetRange(65,80);
	//##//	//			//yield_inc_SBL[ibin][ibin2][ibin3]->Reset();
	//##//	//			//yield_inc_SBL[ibin][ibin2][ibin3]->SetBins(15, 1.5, 3.0);
	//##//	//			//yield_inc_SBR[ibin][ibin2][ibin3]->SetAxisRange(-1.5, 5,"X");	
	//##//	//			//yield_inc_SBR[ibin][ibin2][ibin3]->GetXaxis()->SetRange(46,60);
	//##//				//yield_inc_SBL[ibin][ibin2][ibin3]->Draw("same");
	//##//				//yield_inc_SBR[ibin][ibin2][ibin3]->Draw();
	//##//
	//##//		if(ibin<3){
	//##//					TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
	//##//					centtex->SetNDC();
	//##//					centtex->Draw();
	//##//					TLatex *pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
	//##//					pttex->SetNDC();
	//##//					pttex->Draw();
	//##//					TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
	//##//					Ajtex->SetNDC();
	//##//					Ajtex->Draw();
	//##//					TLatex *Jettex = new TLatex(0.70, 0.85, JetLabel_labels[0]);
	//##//					Jettex->SetNDC();
	//##//					Jettex->Draw();
	//##//
	//##//				}
	//##//				if(ibin==3){
	//##//					TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	//##//					centtex->SetNDC();
	//##//					centtex->Draw();
	//##//					TLatex *pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	//##//					pttex->SetNDC();
	//##//					pttex->Draw();
	//##//					TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
	//##//					Ajtex->SetNDC();
	//##//					Ajtex->Draw();
	//##//					TLatex *Jettex = new TLatex(0.70, 0.85, JetLabel_labels[0]);
	//##//					Jettex->SetNDC();
	//##//					Jettex->Draw();
	//##//
	//##//				}
	//##////return 0;
	//##/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//##//				cc6->cd(4*(ibin3+1) - ibin);
	//##//			//	TH1F *ForPlotting = new TH1F ("ForPlotting", "ForPlotting%d", 100, -5, 5);
	//##//			//	ForPlotting->Draw();
	//##//	//			return 0;
	//##//				yield_inc[ibin][ibin2][ibin3]->ProjectionX("ResultEta%d",20, 80)->DrawCopy();
	//##//				yield_inc[ibin][ibin2][ibin3]->GetXaxis()->SetRange(20, 80);
	//##//				yield_inc_SBL[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionX("ResultSBL%d",21,80 );  // setting the range doesn't make much difference
	//##//				//yield_inc_SBL[ibin][ibin2][ibin3]->SetAxisRange(-5.0, 1.5,"X");	
	//##//				yield_inc_SBL[ibin][ibin2][ibin3]->GetXaxis()->SetRange(0,15 );							
	//##//				//yield_inc_SBL_corrected[ibin][ibin2][ibin3]->GetXaxis()->SetRange(21,35);	
	//##//				yield_inc_SBL[ibin][ibin2][ibin3]->SetMarkerStyle(20);
	//##//				yield_inc_SBL[ibin][ibin2][ibin3]->SetMarkerColor(kRed);
	//##//				yield_inc_SBL[ibin][ibin2][ibin3]->SetLineColor(kRed);	
	//##//				yield_inc_SBL[ibin][ibin2][ibin3]->DrawCopy();
	//##//
	//##//				//yield_inc_SBL[ibin][ibin2][ibin3]->Fit("pol0","0","",-3.0,-1.5);
	//##//				yield_inc_SBL[ibin][ibin2][ibin3]->Fit("pol0","0","",-3.0,-1.5); // MZ_Problem cannot update 
	//##//			//	yield_inc_SBL[ibin][ibin2][ibin3]->DrawCopy("same");
	//##//			        A[ibin][ibin2][ibin3] = pol0->GetParameter(0);
	//##//				cout << "A = " << A[ibin][ibin2][ibin3] << endl;
	//##//
	//##//				yield_inc_SBR[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionX("ResultSBR%d",21,80 );  // setting the range doesn't make much difference				
	//##//				yield_inc_SBR[ibin][ibin2][ibin3]->SetAxisRange(-1.5, 5,"X");						 
	//##//			        yield_inc_SBR[ibin][ibin2][ibin3]->Draw("same");				
	//##//			        yield_inc_SBR[ibin][ibin2][ibin3]->SetMarkerStyle(20);
	//##//			        yield_inc_SBR[ibin][ibin2][ibin3]->SetMarkerColor(kBlue);
	//##//				yield_inc_SBR[ibin][ibin2][ibin3]->GetXaxis()->SetRange(46,60);				
	//##//				yield_inc_SBR[ibin][ibin2][ibin3]->Fit("pol0","0","",1.5,3.0);
	//##//				yield_inc_SBR[ibin][ibin2][ibin3]->DrawCopy("same");
	//##//				
	//##//			         B[ibin][ibin2][ibin3] = pol0->GetParameter(0);				
	//##//				cout << "B = " << B[ibin][ibin2][ibin3] << endl;
	//##//
	//##//				C[ibin][ibin2][ibin3] = (B[ibin][ibin2][ibin3] + A[ibin][ibin2][ibin3]) / 2. ;
	//##//				cout << " C = " << C[ibin][ibin2][ibin3] << endl;
	//##//	if(ibin<3){
	//##//					TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
	//##//					centtex->SetNDC();
	//##//					centtex->Draw();
	//##//					TLatex *pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
	//##//					pttex->SetNDC();
	//##//					pttex->Draw();
	//##//				}
	//##//				if(ibin==3){
	//##//					TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	//##//					centtex->SetNDC();
	//##//					centtex->Draw();
	//##//					TLatex *pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	//##//					pttex->SetNDC();
	//##//					pttex->Draw();
	//##//				}
	//##//
	//##//
	//##//				cc8->cd((4*(ibin3+1) - ibin));
	//##//				//Signal[ibin][ibin2][ibin3] =  (TH1D*) yield_inc[ibin][ibin2][ibin3]->ProjectionY(Form("Signal_copy%d%d%d",ibin,ibin2,ibin3),36, 65);
	//##//				//Signal[ibin][ibin2][ibin3]->Rebin(4);
	//##//				//Signal[ibin][ibin2][ibin3]->SetMarkerColor(kBlack);
	//##//				//Signal[ibin][ibin2][ibin3]->SetMarkerStyle(20);
	//##//	
	//##//				SignalMinusRL[ibin][ibin2][ibin3]= (TH1D*)Signal[ibin][ibin2][ibin3]->Clone((TString)("signal_MinusRL"+CBin_strs[ibin]+"_"+TrkPtBin_strs[ibin3]));
	//##//	
	//##//				 SignalMinusRL[ibin][ibin2][ibin3]->Add( PhiProj_L[ibin][ibin2][ibin3], -1);
	//##//   				 SignalMinusRL[ibin][ibin2][ibin3]->SetAxisRange(-1.5, 1.5,"X");
	//##//				 SignalMinusRL[ibin][ibin2][ibin3]->SetMinimum(0);
	//##//				// SignalMinusRL[ibin][ibin2][ibin3]->SetMaximum(30000);
	//##//				//SignalMinusRL[ibin][ibin2][ibin3]->Scale(1.0 / (Aj[ibin]->GetEntries()));   //MZ it was scaled in cc5
	//##//				
	//##//				 SignalMinusRL[ibin][ibin2][ibin3]->SetMaximum(10);
	//##//				 SignalMinusRL[ibin][ibin2][ibin3]->SetMinimum(-.50);
	//##//				 SignalMinusRL[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#phi");
	//##//				 SignalMinusRL[ibin][ibin2][ibin3]->Draw();
	//##//				
	//##//				if(ibin<3){
	//##//					TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
	//##//					centtex->SetNDC();
	//##//					centtex->Draw();
	//##//					TLatex *pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
	//##//					pttex->SetNDC();
	//##//					pttex->Draw();
	//##//					TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
	//##//					Ajtex->SetNDC();
	//##//					Ajtex->Draw();
	//##//					TLatex *Jettex = new TLatex(0.70, 0.85, JetLabel_labels[0]);
	//##//					Jettex->SetNDC();
	//##//					Jettex->Draw();
	//##//
	//##//				}
	//##//				if(ibin==3){
	//##//					TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	//##//					centtex->SetNDC();
	//##//					centtex->Draw();
	//##//					TLatex *pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	//##//					pttex->SetNDC();
	//##//					pttex->Draw();
	//##//					TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
	//##//					//TLatex *Ajtex = new TLatex(0.80,0.90,AjMax);
	//##//					Ajtex->SetNDC();
	//##//					Ajtex->Draw();
	//##//					TLatex *Jettex = new TLatex(0.70, 0.85, JetLabel_labels[0]);
	//##//					Jettex->SetNDC();
	//##//					Jettex->Draw();
	//##//
	//##//				}
	//##//			  //EtaProjDummy[ibin][ibin2][ibin3]->Fit("flat_fit","0","",-3.,-1.5);
	//##//			 //left_level = flat_fit->GetParameter(0);
	//##//
	//##//			  //EtaProjDummy[ibin][ibin2][ibin3]->Fit("flat_fit","0","",1.5,3.);
	//##//	  
	//##//	  		//	right_level =  flat_fit->GetParameter(0);
	//##//
	//##//	  		//	bg_fit[ibin][ibin2][ibin3] = new TF1("bg_fit","[0]+x-x",-3.,3.);
	//##//	
	//##//	  //bg_fit[ibin][ibin2][ibin3]->SetParameter(0,(left_level+right_level)/2.);
	//##//	//		bg_fit[ibin][ibin2][ibin3]->Draw();
	//##//			//	if((4*(ibin3+1) - ibin) == 16)
	//##//			//	return 6;
	//##//
	//##//				//if((4*(ibin3+1) - ibin) > 1)
	//##//			//	cout << 4*(ibin3+1) - ibin << endl;
	//##//				//if((4*(ibin3+1) - ibin) < 8)
	//##//			//	if((4*(ibin3+1) - ibin) == 16)
	//##//				//if((4*(ibin3+1) - ibin) % 4 == 0)
	//##//				//break;
	//##//			//	return 0;
	//##//	//			c2->cd();					
	//##//				// yield_inc_SBL[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionY("ResultSBL%d",21,35 )->Draw();
	//##//				//yield_inc_SBL[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionX("ResultSBL%d",1,50 )->Draw();
	//##//		//		return 2;
	//##//				//yield_inc_SBL[ibin][ibin2][ibin3]->GetXaxis()->SetRange(21, 35);
	//##//				//yield_inc_SBL[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-3.0, -1.5);
	//##//	//			yield_inc_SBL[ibin][ibin2][ibin3]->Draw();
	//##//	//			//return 3;				
	//##//	//			//yield_inc_SBR[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionY(Form("ResultSBR%d",ibin,ibin2,ibin3),66,80 );
	//##//	//			yield_inc_SBR[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionX(Form("ResultSBR%d",ibin,ibin2,ibin3),66,80 );
	//##//	//			yield_inc_SBR[ibin][ibin2][ibin3]->GetXaxis()->SetRange(66,80 );
	//##//
	//##//				//yield_inc_SBR[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(1.5,3.0 );
	//##//				//yield_inc_SBR[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionY("ResultSBR%d",66,80 );
	//##//				//return 0;				
	//##//
	//##//				c3->cd();
	//##//				yield_inc_SBR[ibin][ibin2][ibin3]->Draw();
	//##//			//	c1->Update();
	//##//	
	//##//				//return 0;				
	//##//				yield_inc_SBL_ForRatio[ibin][ibin2][ibin3] = yield_inc[ibin][ibin2][ibin3]->ProjectionX(Form("ResultSBL_Ratio%d",ibin,ibin2,ibin3),21,35 );
	//##//				yield_inc_SBL_ForRatio[ibin][ibin2][ibin3]->GetXaxis()->SetRange(21,35 );
	//##//				c4->cd();	
	//##//				yield_inc_SBL_ForRatio[ibin][ibin2][ibin3]->Draw();		
	//##//				
	//##//				//return 0;
	//##// 	
	//##//				yield_inc_SBL[ibin][ibin2][ibin3]->Add(yield_inc_SBR[ibin][ibin2][ibin3]);
	//##//				c5->cd();
	//##//				yield_inc_SBL[ibin][ibin2][ibin3]->Draw();
	//##//				c5->Update();
	//##//				//return 5; 
	//##//					
	//##//				yield_inc_SBL_ForRatio[ibin][ibin2][ibin3]->Divide(yield_inc_SBR[ibin][ibin2][ibin3]);
	//##//				yield_inc_SBL_ForRatio[ibin][ibin2][ibin3]->Draw();
	//##//				//return 6;
	//##//				//yield_inc_SBL_ForRatio[ibin][ibin2][ibin3]->Fit("pol0","","",-0.5,0.5); 	
	//##//				   //	yield_inc_SBL_ForRatio[ibin][ibin2][ibin3]->GetFunction("pol0")->GetParameter(0) ;
	//##//				//SideBandRatio[ibin][ibin2][ibin3] =  	yield_inc_SBL_ForRatio[ibin][ibin2][ibin3]->GetFunction("pol0")->GetParameter(0); 	
	//##//				//cout << "SideBandRatio = " << SideBandRatio[ibin][ibin2][ibin3] << endl;	
	//##//
	//##//                                  //Tais m7allel!!// yield_inc_Signal[ibin][ibin2][ibin3]       =new TH1D(histname,"",100,5,5);
	//##//
	//##//				yield_inc_Signal[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionY(Form("me_proj_temp_%d%d%d",ibin,ibin2,ibin3),36,65);
	//##//				yield_inc_Signal_ForBinCheck[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_temp_%d%d%d",ibin,ibin2,ibin3),36,65);
	//##//				yield_corrected_inc_SBR_corrected[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionY(Form("me_proj_temp_%d%d%d",ibin,ibin2,ibin3),36,65);
	//##//				yield_corrected_inc_SBL_corrected[ibin][ibin2][ibin3]	=new TH1D(histname,"",100,5,5);
	//##//				yield_corrected_inc_SBL_corrected[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionY(Form("me_proj_temp_%d%d%d",ibin,ibin2,ibin3),36,65);
	//##//				
	//##//				yield_inc_Signal[ibin][ibin2][ibin3]->GetXaxis()->SetRange(36,65 );
	//##//				yield_corrected_inc_SBR_corrected[ibin][ibin2][ibin3]->GetXaxis()->SetRange(21,35 );
	//##//				yield_corrected_inc_SBR_corrected[ibin][ibin2][ibin3]->Draw();
	//##//				yield_corrected_inc_SBL_corrected[ibin][ibin2][ibin3]->GetXaxis()->SetRange(66,80 );
	//##//				yield_corrected_inc_SBL_corrected[ibin][ibin2][ibin3]->Draw("same");
	//##//				
	//##//				//yield_inc_Signal[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionY("Result_Signal%d",36,65);
	//##//				//result_proj_OnlyLeading->cd(4*(ibin3+1) - ibin);
	//##//				//yield_inc_Signal[ibin][ibin2][ibin3]->Draw();
	//##//				//return 0;
	//##//				yield_inc[ibin][ibin2][ibin3]->ProjectionY("yield_inc_Signal%d",36,65 );
	//##//
	//##//				//c6->Update();	
	//##//				
	//##//				yield_inc_Signal[ibin][ibin2][ibin3]->Add(yield_inc_SBL[ibin][ibin2][ibin3],-1);
	//##//				yield_inc_Signal[ibin][ibin2][ibin3]->Scale(1.0 / (3.0 *  Aj[ibin]->GetEntries() * PhiBinWidth));
	//##//				yield_inc_Signal[ibin][ibin2][ibin3]->SetAxisRange(-1.5, 1.5,"X");
	//##//		//		c6->cd(4*(ibin3+1) - ibin);
	//##//				//result_proj_OnlyLeading->cd(4*(ibin3+1) - ibin)	;											
	//##//				 yield_inc_Signal[ibin][ibin2][ibin3]->Draw();
	//##//	if(ibin<3){
	//##//					TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
	//##//					centtex->SetNDC();
	//##//					centtex->Draw();
	//##//					TLatex *pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
	//##//					pttex->SetNDC();
	//##//					pttex->Draw();
	//##//				}
	//##//				if(ibin==3){
	//##//					//cout << ibin << endl;
	//##//					TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	//##//					centtex->SetNDC();
	//##//					centtex->Draw();
	//##//					TLatex *pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	//##//					pttex->SetNDC();
	//##//					pttex->Draw();
	//##//				}
	//##//
	//##//			
	//##//				//cc7->cd();
	//##//				cc7->cd(4*(ibin3+1)-ibin);
	//##//				 fa->SetParameter(0, C[ibin][ibin2][ibin3] );	
	//##//			
	//##//
	//##//				//yield_inc_proj[ibin][ibin2][ibin3]->DrawCopy();	
	//##//				//yield_inc[ibin][ibin2][ibin3]->GetXaxis()->SetRange(20, 80);										
	//##//				yield_inc_proj[ibin][ibin2][ibin3]->SetAxisRange(C[ibin][ibin2][ibin3] * 1.1,  C[ibin][ibin2][ibin3] * 0.9,"Y");
	//##//				//yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(C[ibin][ibin2][ibin3] * 1.1);
	//##//				//yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(C[ibin][ibin2][ibin3] * 0.9);
	//##//				yield_inc_proj[ibin][ibin2][ibin3]->Draw();
	//##//			fa->Draw("SAME");
	//##//	
	//##//			
	//##//	if(ibin<3){
	//##//					TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
	//##//					centtex->SetNDC();
	//##//					centtex->Draw();
	//##//					TLatex *pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
	//##//					pttex->SetNDC();
	//##//					pttex->Draw();
	//##//				}
	//##//				if(ibin==3){
	//##//					//cout << ibin << endl;
	//##//					TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	//##//					centtex->SetNDC();
	//##//					centtex->Draw();
	//##//					TLatex *pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	//##//					pttex->SetNDC();
	//##//					pttex->Draw();
	//##//				}
	//##//				
	//##//		//	if((4*(ibin3+1) - ibin) == 16)		
	//##//		//return 0;
	//##//	
	//##//	//			c8->cd(4*(ibin3+1)-ibin);
	//##//	//			//yield_corrected_inc_SBL_corrected[ibin][ibin2][ibin3]->Scale(1.0 / (1.5 *  Aj[ibin]->GetEntries() * PhiBinWidth));	
	//##//	//			yield_corrected_inc_SBL_corrected[ibin][ibin2][ibin3]->SetAxisRange(-5.0, 1.5,"X");				
	//##//	//			//yield_corrected_inc_SBL_corrected[ibin][ibin2][ibin3]->GetXaxis()->SetRange(21,35);				
	//##//	//			yield_corrected_inc_SBL_corrected[ibin][ibin2][ibin3]->Draw();
	//##//	//			yield_corrected_inc_SBL_corrected[ibin][ibin2][ibin3]->Fit("pol0","0","",-3.0,-1.5);
	//##//	//		         A[ibin][ibin2][ibin3] = pol0->GetParameter(0);
	//##//	//			cout << "A = " << A[ibin][ibin2][ibin3] << endl;
	//##//
	//##//	//			c9->cd(4*(ibin3+1)-ibin);
	//##//	//			yield_corrected_inc_SBR_corrected[ibin][ibin2][ibin3]->SetAxisRange(-1.5, 5,"X");				
	//##//	//			yield_corrected_inc_SBL_corrected[ibin][ibin2][ibin3]->Draw();
	//##//	//			yield_corrected_inc_SBL_corrected[ibin][ibin2][ibin3]->Fit("pol0","0","",1.5,3.0);
	//##//	//			yield_corrected_inc_SBL_corrected[ibin][ibin2][ibin3]->Draw();
	//##//	//			
	//##//	//		         B[ibin][ibin2][ibin3] = pol0->GetParameter(0);				
	//##//	//			cout << "B = " << B[ibin][ibin2][ibin3] << endl;
	//##//
	//##//	//			C[ibin][ibin2][ibin3] = (B[ibin][ibin2][ibin3] + A[ibin][ibin2][ibin3]) / 2. ;
	//##//			//	   fa->SetParameter(0, C[ibin][ibin2][ibin3] );				
	//##//
	//##//	//			c10->cd(4*(ibin3+1)-ibin);
	//##//	//			yield_inc_Signal[ibin][ibin2][ibin3]->SetAxisRange(-5, 5,"X");
	//##//	//			yield_inc_Signal[ibin][ibin2][ibin3]->SetAxisRange(C[ibin][ibin2][ibin3] * 1.1,  C[ibin][ibin2][ibin3] * 0.9,"Y");
	//##//	//			//yield_inc_Signal[ibin][ibin2][ibin3]->SetAxisRange(C[ibin][ibin2][ibin3] + 0.1,  C[ibin][ibin2][ibin3] + 0.1,"Y");
	//##//	//			//yield_inc_Signal[ibin][ibin2][ibin3]->SetAxisRange( -0.1,  0.1,"Y");
	//##//	//			//yield_inc_Signal[ibin][ibin2][ibin3]->SetAxisRange( -0.2,  0.2,"Y");
	//##//	//			yield_inc_Signal[ibin][ibin2][ibin3]->Draw();	
	//##//	//		//	fa->Draw("SAME");
	//##//
	//##//
	//##//	//if(ibin<3){
	//##//	//				TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
	//##//	//				centtex->SetNDC();
	//##//	//				centtex->Draw();
	//##//	//				TLatex *pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
	//##//	//				pttex->SetNDC();
	//##//	//				pttex->Draw();
	//##//	//			}
	//##//	//			if(ibin==3){
	//##//	//				//cout << ibin << endl;
	//##//	//				TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	//##//	//				centtex->SetNDC();
	//##//	//				centtex->Draw();
	//##//	//				TLatex *pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	//##//	//				pttex->SetNDC();
	//##//	//				pttex->Draw();
	//##//	//			}
	//##//
	//##//	
	//##//				//yield_inc_Signal[ibin][ibin2][ibin3]->Draw();	
	//##//				//me_proj_canvas->Update();
	//##//				//yield_inc_Signal[ibin][ibin2][ibin3]->Rebin(5);				
	//##//	///			c7->cd();
	//##//				//yield_inc_Signal[ibin][ibin2][ibin3]->Draw();
	//##//				//result_proj_OnlyLeading->Update();
	//##//				yield_inc_Signal[ibin][ibin2][ibin3]->Write();
	//##//				
	//##//
	//##//				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Write();
	//##//				yield_inc[ibin][ibin2][ibin3]->Write();
	//##//				hJetTrackMELeading[ibin][ibin2][ibin3]->Write();
	//##//			//	cout << "Canvas panel" << 4*(ibin3+1)-ibin << endl;
	//##//					//me_proj_canvas->Modified();	
	//##//		//	result_proj_canvas->cd(4*(ibin3+1)-ibin);
	//##//			
	//##//			//cout << ibin << endl:
	//##//			//result_proj_canvas->cd(ibin3,ibin);
	//##//
	//##//				int phil = yield_inc[ibin][ibin2][ibin3]->GetYaxis()->FindBin(-.99999);
	//##//				int phir = yield_inc[ibin][ibin2][ibin3]->GetYaxis()->FindBin(.99999);
	//##//
	//##//				yield_inc_proj[ibin][ibin2][ibin3]= (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionX((TString)("yield_inc_proj_"+CBin_strs[ibin]+"_"+TrkPtBin_strs[ibin3]),phil,phir);
	//##//				yield_inc_proj_signal[ibin][ibin2][ibin3]= (TH1D*)yield_inc_Signal[ibin][ibin2][ibin3]->Clone((TString)("yield_inc_proj_signal"+CBin_strs[ibin]+"_"+TrkPtBin_strs[ibin3]));
	//##//				test_me_inc_proj[ibin][ibin2][ibin3]= (TH1D*)hJetTrackMELeading[ibin][ibin2][ibin3]->ProjectionX((TString)("test_me_inc_proj_"+CBin_strs[ibin]+"_"+TrkPtBin_strs[ibin3]),phil,phir);
	//##//			
	//##//				
	//##//
	//##//				//yield_inc_proj[ibin][ibin2][ibin3]->SetMarkerStyle(10);
	//##//	//			me_projeta[ibin][ibin2][ibin3]->SetMarkerStyle(10);
	//##//				//yield_inc_proj[ibin][ibin2][ibin3]->SetMarkerSize(1);
	//##//	//			me_projeta[ibin][ibin2][ibin3]->SetMarkerColor(kBlack);
	//##//
	//##//	//			result_proj_OnlyLeading->cd(4*(ibin3+1)-ibin);			
	//##//				
	//##//	//			yield_inc_proj[ibin][ibin2][ibin3]->Draw();
	//##//	//			result_proj_OnlyLeading->Update();
	//##//	//			me_projeta[ibin][ibin2][ibin3]->SetMarkerSize(1);
	//##//				//yield_inc_proj[ibin][ibin2][ibin3]->SetMarkerColor(kBlack);
	//##//	//			me00[ibin][ibin2][ibin3];
	//##//				//me_projeta[ibin][ibin2][ibin3]->GetYaxis()->SetRange(0.9 * me00[ibin][ibin2][ibin3], 1.1 * me00[ibin][ibin2][ibin3]); // didn't work
	//##//	//			me_projeta[ibin][ibin2][ibin3]->SetAxisRange(0.5 * me00[ibin][ibin2][ibin3], 1.5 * me00[ibin][ibin2][ibin3],"Y");
	//##//	//			me_projeta[ibin][ibin2][ibin3]->SetAxisRange(-5.0, 5.0,"X");
	//##//				//me_proj_canvas->Update();
	//##//				
	//##//	//			me_projeta[ibin][ibin2][ibin3]->Draw();
	//##//				//me_proj_canvas->Update();
	//##//
	//##//				//test_me_inc_proj[ibin][ibin2][ibin3]->Scale(1/50.);
	//##//
	//##//	//			test_me_inc_proj[ibin][ibin2][ibin3]->SetMarkerStyle(10);
	//##//	//			test_me_inc_proj[ibin][ibin2][ibin3]->SetMarkerSize(1);
	//##//	//			test_me_inc_proj[ibin][ibin2][ibin3]->SetMarkerColor(kRed);
	//##//	//			test_me_inc_proj[ibin][ibin2][ibin3]->SetLineColor(kRed);
	//##//	//			test_me_inc_proj[ibin][ibin2][ibin3]->Draw("same");
	//##//				//me_proj_canvas->Update();
	//##//
	//##//	//			if(ibin<3){
	//##//	//				centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
	//##//	//				centtex->SetNDC();
	//##//	//				centtex->Draw();
	//##//	//				pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
	//##//	//				pttex->SetNDC();
	//##//	//				pttex->Draw();
	//##//	//			}
	//##//	//			if(ibin==3){
	//##//	//				centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	//##//	//				centtex->SetNDC();
	//##//	//				centtex->Draw();
	//##//	//				pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	//##//	//				pttex->SetNDC();
	//##//	//				pttex->Draw();
	//##//	//			}
	//##//
	//##//				//me_proj_sideBands->cd(4*(ibin3+1)-ibin);
	//##//			
	//##//	//			yield_inc_SBL_ForRatio[ibin][ibin2][ibin3]->SetMarkerStyle(1);
	//##//	//			yield_inc_SBL_ForRatio[ibin][ibin2][ibin3]->SetMarkerSize(10);
	//##//	//			yield_inc_SBL_ForRatio[ibin][ibin2][ibin3]->SetMarkerColor(kBlue);
	//##//	//			yield_inc_SBL_ForRatio[ibin][ibin2][ibin3]->SetLineColor(kBlack);
	//##//				//yield_inc_SBL_ForRatio[ibin][ibin2][ibin3]->SetAxisRange(0.9* SideBandRatio[ibin][ibin2][ibin3], 1.1 * SideBandRatio[ibin][ibin2][ibin3],"Y");
	//##//	//			yield_inc_SBL_ForRatio[ibin][ibin2][ibin3]->SetAxisRange(0.8 , 1.2 ,"Y");
	//##//	//			yield_inc_SBL_ForRatio[ibin][ibin2][ibin3]->SetAxisRange(-5, 5,"X");
	//##//
	//##//				//yield_inc_SBL_ForRatio[ibin][ibin2][ibin3]->Draw("hist p");
	//##//	//			yield_inc_SBL_ForRatio[ibin][ibin2][ibin3]->Draw();
	//##//	//if(ibin<3){
	//##//	//				TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
	//##//	//				centtex->SetNDC();
	//##//	//				centtex->Draw();
	//##//	//				TLatex *pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
	//##//	//				pttex->SetNDC();
	//##//	//				pttex->Draw();
	//##//	//			}
	//##//	//			if(ibin==3){
	//##//	//				TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	//##//	//				centtex->SetNDC();
	//##//	//				centtex->Draw();
	//##//	//				TLatex *pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	//##//	//				pttex->SetNDC();
	//##//	//				pttex->Draw();
	//##//	//			}
	//##//			//	result_proj_OnlyLeading->cd(4*(ibin3+1)-ibin);
	//##//		//	yield_inc_Signal[ibin][ibin2][ibin3]->SetAxisRange->(-1.5, 0, "X");
	//##//		//	yield_inc_Signal[ibin][ibin2][ibin3]->Draw();
	//##//				

			}
		}
	}

	//me_proj_canvas->SaveAs((TString)("All_ME_Projections_Inclusive_"+datalabel+"_Leading_Over22.png"));
	//result_proj_canvas->SaveAs((TString)("Corrected_dEta_Incluisve_Projections_"+datalabel+"_Leading_Over22.png"));

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

	//hists_hJetTrackSignalBackgroundCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("SBL",21,35 )->Draw();
	//hists_hJetTrackSignalBackgroundCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("SBR",66,80 )->Draw();
	//SBL->Add(SBR);
	//SBL->Draw();

	//hists_hJetTrackSignalBackgroundCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("Signal",36,65 )->Draw();
	//Signal->Add(SBL,-1);
	// Signal->Draw();

	//hists_hJetTrackSignalBackgroundCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("Signal",36,65 )->Draw();

	//TCanvas *me_proj_canvas = new TCanvas("me_proj_canvas","",0,0,1500,1600);
	//  me_proj_canvas->Divide(4,4,0.0000,0.0000);

	return i;
	}                                           
