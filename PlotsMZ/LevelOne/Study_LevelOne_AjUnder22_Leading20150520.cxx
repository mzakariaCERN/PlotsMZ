
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

	int Study_LevelOne_AjUnder22_Leading20150520(int i=0)
	{
	TFile  *finbg, *fout_inc ;

	//TFile *fin = TFile::Open("Data2011_trkPtCut1_All_AJ11.root");  // The input file
	//TFile *fin = TFile::Open("Data2011_All_AJ22.root");  // The input file
	//TFile *fin = TFile::Open("Data2011_All_AJOver22.root");  // The input file
	//TFile *fin = TFile::Open("Data2011_AjUnder22.root");  // The input file
	//TFile *fin = TFile::Open("Data2011_AjUnder22_Only13.root");  // The input file
	TFile *fin = TFile::Open("../Data2011_AjUnder22_Only13_20150416.root");  // The input file
	//TFile *fin = TFile::Open("FromH/PbPb_Leading_Correlations.root");  // The input file
	//TFile *fin = TFile::Open("Data2011_All_AJ33.root");  // The input file
	//TFile *fin = TFile::Open("/home/mzakaria3/Documents/Research/MyGitProject/PlotsMZ/Data2011_trkPtCut1_All_AJ11.root");  // The input file
	//	TFile *fin = TFile::Open("/home/mzakaria3/Documents/Research/MyGitProject/PlotsMZ/Data2011_trkPtCut1_All_Aj11_20150318.root");  // The input file
	//TFile *fin = TFile::Open("/home/mzakaria3/Documents/Research/MyGitProject/PlotsMZ/Data2011_trkPtCut1_All_Aj11_20150319.root");  // The input file
	fout_inc = new TFile("LevelOne_PbPb_AjUnder22_Correlations.root","RECREATE"); // The output file

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
	TString AjBin_labels[1] = {"Aj < 0.22"};
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


	TCanvas *cc1 = new TCanvas("cc1", "CC1_me_projeta", 0, 0, 1500, 1500);	
	TCanvas *cc2 = new TCanvas("cc2", "CC2_me_projphi", 0, 0, 1500, 1500);	
	TCanvas *cc3 = new TCanvas("cc3", "CC3_SignalNoScale", 0, 0, 1500, 1500);	
	TCanvas *cc4 = new TCanvas("cc4", "CC4_S", 0, 0, 1500, 1500);	
	TCanvas *cc5 = new TCanvas("cc5", "CC5_ME", 0, 0, 1500, 1500);	
	TCanvas *cc6 = new TCanvas("cc6", "CC6_MEnormalized", 0, 0, 1500, 1500);	
	TCanvas *cc7 = new TCanvas("cc7", "CC6j_LBoverRB", 0, 0, 1500, 1500);	
	TCanvas *cc8 = new TCanvas("cc8", "CC7_EtaBandSymmetry", 0, 0, 1500, 1500);	
	TCanvas *cc9 = new TCanvas("cc9", "CC8_EtaBandSymmetry_HT", 0, 0, 1500, 1500);	


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
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Signal_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
				
			//	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

				hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
				
				//hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

				Correlation[ibin][ibin2][ibin3] = (TH2D*)  fin->Get((TString) (desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Scaled_Inv_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				cc1->cd(4*(ibin3+1)-ibin);
				// cc1->SetTickx(1);
				//cc1->SetGridy();
				//cc1->SetTopMargin(0.001);
				//cc1->SetBottomMargin(0.001);
				//cc1->Modified();
				me_projeta[ibin][ibin2][ibin3] =  hJetTrackMELeading[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_eta_cc1_%d%d%d",ibin,ibin2,ibin3),1,100); 				
				//me_projeta[ibin][ibin2][ibin3]->Draw() ;
	
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
				//me_projeta[ibin][ibin2][ibin3]->DrawCopy() ;
				me_projeta[ibin][ibin2][ibin3]->Draw() ;
		
		
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

	cc1->SaveAs("cc1_ME_EtaProj_Leading_Under22.png");
// return 1;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	
	//##/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//##//				
				cc2->cd(4*(ibin3+1) - ibin);
				me_projeta[ibin][ibin2][ibin3] =  hJetTrackMELeading[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_eta_cc2_%d%d%d",ibin,ibin2,ibin3),1,100); // the %d means a 4bit integer
			//	me_projeta[ibin][ibin2][ibin3]->Draw();
				NEtaBins =  me_projeta[ibin][ibin2][ibin3]->GetNbinsX() ;
				//me_projeta[ibin][ibin2][ibin3] -> Fit("pol0","","",-.2,.2);
				//me00[ibin][ibin2][ibin3] = me_projeta[ibin][ibin2][ibin3]->GetFunction("pol0")->GetParameter(0);				
				me_projeta[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");
				me_projeta[ibin][ibin2][ibin3]->SetAxisRange(-2,2, "X");
				me_projeta[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-2.,2.);	
				me_projeta[ibin][ibin2][ibin3]->GetXaxis() -> SetTitleOffset(0.8);
				me_projeta[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				me_projeta[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				me_projeta[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.06);
				me_projeta[ibin][ibin2][ibin3]->GetXaxis() -> CenterTitle();
				

				me_projeta[ibin][ibin2][ibin3]->Scale (1./me00[ibin][ibin2][ibin3]);

				me_projeta[ibin][ibin2][ibin3]->SetMaximum(me_projeta[ibin][ibin2][ibin3]->GetMaximum()+0.15);
				me_projeta[ibin][ibin2][ibin3]->Clone((TString)("Eta_Proj_Normalized"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

				me_projeta[ibin][ibin2][ibin3]->Draw();
				
							//me_projeta[ibin][ibin2][ibin3]->Draw();
	
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

					cc2->SaveAs("cc2_MEOverME00_Leading_Under22.png");
//if((4*(ibin3+1) - ibin) == 5)
//{
//return 2;
//}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//##//
	//##//		
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////					
				cc3->cd(4*(ibin3+1)-ibin);

	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->CenterTitle();
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->CenterTitle();
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitleOffset(1.6);
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitleOffset(1.6);
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
					TLatex *pttex = new TLatex(0.05,0.83,TrkPtBin_labels[ibin3]);
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
					TLatex *pttex = new TLatex(0.15,0.83,TrkPtBin_labels[ibin3]);
					pttex->SetNDC();
					pttex->Draw();
					TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
					Ajtex->SetNDC();
					Ajtex->Draw();
					TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
					Jettex->SetNDC();
					Jettex->Draw();


				}
	cc3->SaveAs("cc3_SignalNoScale_Leading_Under22.png");
//return 3;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

				cc4->cd(4*(ibin3+1)-ibin);
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Scaled_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

				PhiBinWidth =  me_projphi[ibin][ibin2][ibin3]->GetBinWidth(1) ;
				EtaBinWidth =  hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->GetBinWidth(1) ;
				PhiBinWidth =  hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->GetBinWidth(1) ;
				cout << "Eta Bin Width = " << 	 hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->GetBinWidth(1) << endl;
				cout << "####################ibin = " << ibin << endl;

				//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / (3.0 *  Aj[ibin]->GetEntries() * PhiBinWidth));
				//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / (0.1 *  Aj[ibin]->GetEntries() * PhiBinWidth));
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1.0 / (EtaBinWidth *  Aj[ibin]->GetEntries() * PhiBinWidth));
				hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Clone((TString) ("Scaled_Signal_Inv_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
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
					TLatex *pttex = new TLatex(0.05,0.83,TrkPtBin_labels[ibin3]);
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
					TLatex *pttex = new TLatex(0.15,0.83,TrkPtBin_labels[ibin3]);
					pttex->SetNDC();
					pttex->Draw();
					TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
					Ajtex->SetNDC();
					Ajtex->Draw();
					TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
					Jettex->SetNDC();
					Jettex->Draw();


				}
			//This is an unvariant histogram			
			cc4->SaveAs("cc4_S_Leading_Under22.png");
					//if(ibin != 0)
			//return 4;
			
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			cc5->cd(4*(ibin3+1)-ibin);
		
				hJetTrackMELeading_cc5[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
		
				me_projeta[ibin][ibin2][ibin3] =  hJetTrackMELeading[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_eta_cc5_%d%d%d",ibin,ibin2,ibin3),1,100);	
				///me_projeta[ibin][ibin2][ibin3]->Scale(1. / 100.);	
				me_projeta[ibin][ibin2][ibin3]->Scale(1. / 99.);
				me_projeta[ibin][ibin2][ibin3]->Draw();	
				me_projeta[ibin][ibin2][ibin3]->Fit("pol0","","",-.2,.2);
				me_projeta[ibin][ibin2][ibin3]->Draw();	
				me00[ibin][ibin2][ibin3] = me_projeta[ibin][ibin2][ibin3]->GetFunction("pol0")->GetParameter(0);	
				//hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1.0 / (3.0 *  Aj[ibin]->GetEntries() * PhiBinWidth));
			//	hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1.0 / (EtaBinWidth *  Aj[ibin]->GetEntries() * PhiBinWidth));
			//	hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(100./me00[ibin][ibin2][ibin3]);
				hJetTrackMELeading_cc5[ibin][ibin2][ibin3]->Scale(1./me00[ibin][ibin2][ibin3]);
				hJetTrackMELeading_cc5[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
				hJetTrackMELeading_cc5[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
				hJetTrackMELeading_cc5[ibin][ibin2][ibin3]->GetXaxis()-> CenterTitle();
				hJetTrackMELeading_cc5[ibin][ibin2][ibin3]->GetYaxis()-> CenterTitle();
				hJetTrackMELeading_cc5[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleOffset(1.6);
				hJetTrackMELeading_cc5[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleOffset(1.6);
				hJetTrackMELeading_cc5[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				hJetTrackMELeading_cc5[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
				hJetTrackMELeading_cc5[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				hJetTrackMELeading_cc5[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.06);

			//	hJetTrackMELeading_cc5[ibin][ibin2][ibin3]->SetAxisRange(-3., 3.,"X");
				hJetTrackMELeading_cc5[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-3.,3.);	
			//	hJetTrackMELeading_cc5[ibin][ibin2][ibin3]->Scale(100./me00[ibin][ibin2][ibin3]);
				hJetTrackMELeading_cc5[ibin][ibin2][ibin3]->Clone((TString) ("Normalized_Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	
				hJetTrackMELeading_cc5[ibin][ibin2][ibin3]->Draw("LEGO2");

			//	if(ibin<3){
			//		TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
			//		centtex->SetNDC();
			//		centtex->Draw();
			//		TLatex *pttex = new TLatex(0.05,0.83,TrkPtBin_labels[ibin3]);
			//		pttex->SetNDC();
			//		pttex->Draw();
			//		TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
			//		Ajtex->SetNDC();
			//		Ajtex->Draw();
			//		TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
			//		Jettex->SetNDC();
			//		Jettex->Draw();


			//	}
			//	if(ibin==3){
			//		TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
			//		centtex->SetNDC();
			//		centtex->Draw();
			//		TLatex *pttex = new TLatex(0.15,0.83,TrkPtBin_labels[ibin3]);
			//		pttex->SetNDC();
			//		pttex->Draw();
			//		TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
			//		Ajtex->SetNDC();
			//		Ajtex->Draw();
			//		TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
			//		Jettex->SetNDC();
			//		Jettex->Draw();


			//	}
					cc5->Modified();
					cc5->Update();
				
					cc5->SaveAs("cc5_MEOverME00_Notscaled_Leading_Under22.png");
//return 5;

	cc6->cd(4*(ibin3+1)-ibin);
		
				hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
		
	

				hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1.0 / (EtaBinWidth *  Aj[ibin]->GetEntries() * PhiBinWidth));
				hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1./me00[ibin][ibin2][ibin3]);

				hJetTrackMELeading[ibin][ibin2][ibin3]->Clone((TString) ("MEoverME00_Scaled_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis()-> CenterTitle();
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetYaxis()-> CenterTitle();
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleOffset(1.6);
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleOffset(1.6);
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis()-> SetTitleSize(0.06);
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetYaxis()-> SetTitleSize(0.06);
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis()-> SetLabelSize(0.06);
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetYaxis()-> SetLabelSize(0.06);

				hJetTrackMELeading[ibin][ibin2][ibin3]->SetAxisRange(-3., 3.,"X");
				hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-3.,3.);				
				hJetTrackMELeading[ibin][ibin2][ibin3]->Clone((TString) ("Normalized_Scaled_Mixed_Event_Inv_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	
				//hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1./me_projeta[ibin][ibin2][ibin3]);				
				//hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(100./me_projeta[ibin][ibin2][ibin3]);				
			//	hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(100./me00[ibin][ibin2][ibin3]);

				hJetTrackMELeading[ibin][ibin2][ibin3]->Draw("LEGO2");

				if(ibin<3){
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


				}
				if(ibin==3){
					TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
					centtex->SetNDC();
					centtex->Draw();
					TLatex *pttex = new TLatex(0.15,0.83,TrkPtBin_labels[ibin3]);
					pttex->SetNDC();
					pttex->Draw();
					TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
					Ajtex->SetNDC();
					Ajtex->Draw();
					TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
					Jettex->SetNDC();
					Jettex->Draw();


				}

	//	cc6->Modified();
	//	cc6->Update();

		cc6->SaveAs("cc6_MEOverME00_Scaled_Leading_Under22.png");

//return 6;

		cc7->cd(4*(ibin3+1)-ibin);
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


			hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

			//hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(100./me00[ibin][ibin2][ibin3]);	
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



	
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Draw("LEGO2");	
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Clone((TString) ("Scaled_Normalized_Yield_Inv_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

					if(ibin<3){
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


				}
				if(ibin==3){
					TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
					centtex->SetNDC();
					centtex->Draw();
					TLatex *pttex = new TLatex(0.15,0.83,TrkPtBin_labels[ibin3]);
					pttex->SetNDC();
					pttex->Draw();
					TLatex *Ajtex = new TLatex(0.80,0.90,AjBin_labels[0]);
					Ajtex->SetNDC();
					Ajtex->Draw();
					TLatex *Jettex = new TLatex(0.75, 0.85, JetLabel_labels[0]);
					Jettex->SetNDC();
					Jettex->Draw();


				}


				cc7->SaveAs("cc7_S_Scaled_OverMEOverME00_NotScaled_Leading_Under22.png");
//	return 7;
				cc8->cd(4*(ibin3+1)-ibin);

			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");	
	      	//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Draw("LEGO2");   //testing that this is the S/ME from cc6

	      //	yield_inc[ibin][ibin2][ibin3] =  (TH2D*) hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Clone((TString)("Yield_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));	

			hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	      	
			hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

			hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1./me00[ibin][ibin2][ibin3]);	
			//hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1.0 / ( Aj[ibin]->GetEntries() * PhiBinWidth*EtaBinWidth));	
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
			yield_inc_proj[ibin][ibin2][ibin3] = (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionX((Form("ResultEta_cc8_%d%d%d",ibin,ibin2,ibin3)),1, 50);	
			//yield_inc_proj[ibin][ibin2][ibin3] = (TH1D*)hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->ProjectionX((Form("ResultEta%d%d%d",ibin,ibin2,ibin3)),1, 25);	
					yield_inc_proj[ibin][ibin2][ibin3]->Rebin(4);
				//	yield_inc_proj[ibin][ibin2][ibin3]->Scale(EtaBinWidth/4.);   //
					yield_inc_proj[ibin][ibin2][ibin3]->Scale(1./(4. *49.));   //
				//	yield_inc_proj[ibin][ibin2][ibin3]->Scale(EtaBinWidth / ( 4. * 49));   //
						
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

				}
				cc8->SaveAs("cc8_RightBandLeftBandComparison_Leading_Under22.png");
								
//return 8;
				cc9->cd((4*(ibin3+1)-ibin));
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


			hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

			//hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(100./me00[ibin][ibin2][ibin3]);	
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
			//	double MaxSignal[6][6];
			//	MaxSignal[ibin][ibin3] = Signal[ibin][ibin2][ibin3]->GetMaximum();
			//	Signal[ibin][ibin2][ibin3]->SetMaximum(MaxSignal[ibin][ibin3]+30);
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
			//	Signal[ibin][ibin2][ibin3]->SetMinimum(0);
			//	Signal[ibin][ibin2][ibin3]->SetMaximum(1500);
			//	if(ibin == 1 && ibin3 == 0)
			//	Signal[ibin][ibin2][ibin3]->SetMaximum(30);
			//	if(ibin == 2 && ibin3 == 0)
			//	Signal[ibin][ibin2][ibin3]->SetMaximum(15);
			//	if(ibin == 3 && ibin3 == 0)
			//	Signal[ibin][ibin2][ibin3]->SetMaximum(6);
			//	if(ibin == 0 && ibin3 == 1)
			//	Signal[ibin][ibin2][ibin3]->SetMaximum(20);
			//	if(ibin == 1 && ibin3 == 1)
			//	Signal[ibin][ibin2][ibin3]->SetMaximum(16);
			//	if(ibin == 2 && ibin3 == 1)
			//	Signal[ibin][ibin2][ibin3]->SetMaximum(8);
			//	if(ibin == 3 && ibin3 == 1)
			//	Signal[ibin][ibin2][ibin3]->SetMaximum(5);
			//	if(ibin == 0 && ibin3 == 2)
			//	Signal[ibin][ibin2][ibin3]->SetMaximum(6);
			//	if(ibin == 1 && ibin3 == 2)
			//	Signal[ibin][ibin2][ibin3]->SetMaximum(5);
			//	if(ibin ==2 && ibin3 == 2)
			//	Signal[ibin][ibin2][ibin3]->SetMaximum(4);
			//	if(ibin ==3 && ibin3 == 2)
			//	Signal[ibin][ibin2][ibin3]->SetMaximum(4);
			//	if(ibin ==0 && ibin3 ==3)
			//	Signal[ibin][ibin2][ibin3]->SetMaximum(8);
			//	if(ibin ==1 && ibin3 ==3)
			//	Signal[ibin][ibin2][ibin3]->SetMaximum(7.5);
			//	if(ibin ==2 && ibin3 ==3)
			//	Signal[ibin][ibin2][ibin3]->SetMaximum(7.5);
			//	if(ibin ==3 && ibin3 ==3)
			//	Signal[ibin][ibin2][ibin3]->SetMaximum(11);


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
		//hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Draw("LEGO2");
			cc9->SaveAs("cc9_PhiProjectionofSignalandBands_Leading_Under22.png");			
//return 9;
		
			}
		}
	}

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

	return i;
	}                                           
