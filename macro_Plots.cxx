
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

int macro_Plots(int i=0)
{
	TFile  *finbg, *fout_inc ;

	//TFile *fin = TFile::Open("Data2011_trkPtCut1_All_AJ11.root");  // The input file
	//TFile *fin = TFile::Open("Data2011_All_AJ22.root");  // The input file
	TFile *fin = TFile::Open("Data2011_All_AJ33.root");  // The input file
	fout_inc = new TFile("PbPb_Inclusive_Correlations.root","RECREATE"); // The output file

	const int nCBins = 4;   //Centrality
	const int nPtBins = 1;  // Hardest Jet range
	const int nTrkPtBins = 5; //Associated tracks range

	float me00; // used to save the mixing at (Deta,Dphi) = (0,0) (the output of fitting)
//	int me00bin;

	float PtBins[nPtBins+1] = {100, 300};
	TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};
	TLatex *centtex, *pttex; 

	float CBins[nCBins+1] = {0, 20, 60, 100, 200};
	TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
	TString CBin_labels[nCBins] = {"Cent. 0-10%","Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};

	float TrkPtBins[nTrkPtBins+1] = {1, 2, 3, 4, 8, 999};
	TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt999" };
	TString TrkPtBin_labels[nTrkPtBins] = {"1<pT<2","2<pT<3","3<pT<4","4<pT<8","pT>8"};


	TCanvas *me_proj_canvas = new TCanvas("me_proj_canvas","",0,0,1500,1600);
	me_proj_canvas->Divide(4,4,0.0000,0.0000);

	TCanvas *result_proj_canvas = new TCanvas("result_proj_canvas","",0,0,1500,1600);
	result_proj_canvas->Divide(4,4,0.0000,0.0000);



	cout << "the default parameter = " <<  i <<  endl;
	int Logical = 1;

	//if(!Logical)
	//  printf("Not logical value at line number %d in file %s\n", __LINE__, __FILE__);  // used to print the line (debugging)

	TString  datalabel;

	gStyle->SetOptStat(0);  
	gStyle->SetPadBottomMargin(0.15);
	gStyle->SetPadTopMargin   (0.05);
	gStyle->SetPadLeftMargin  (0.15);
	gStyle->SetPadRightMargin (0.05);

	TH2D* hJetTrackSignalBackground[nCBins][nPtBins][nTrkPtBins];  // The signal
	TH2D* hJetTrackME[nCBins][nPtBins][nTrkPtBins];                // The mixed events
	TH2D *yield_inc[nCBins][nPtBins][nTrkPtBins];      	       // used to save a clone of hJetTrackSignalBackground[][][] and then become the Histo of Signal over ME

	TH1D *yield_inc_SBL[nCBins][nPtBins][nTrkPtBins];              // The Y projection of the yield (phi) at the left range
	TH1D *yield_inc_SBR[nCBins][nPtBins][nTrkPtBins];		// The Y projection of the yield (phi) at the right range 
	TH1D *yield_inc_Signal[nCBins][nPtBins][nTrkPtBins];		// The Y projection of the yield (phi) at the center (the signal) 

	TH1D *me_proj[nCBins][nPtBins][nTrkPtBins];
	TH1D *yield_inc_proj[nCBins][nPtBins][nTrkPtBins];
	TH1D *yield_inc_proj_signal[nCBins][nPtBins][nTrkPtBins];
	TH1D *test_me_inc_proj[nCBins][nPtBins][nTrkPtBins];
	TH1D *ProEta;
	TH1D *ProPhi;

	TString desc = "hists";

	for (int ibin=0;ibin<nCBins;ibin++){
		for (int ibin2=0;ibin2<nPtBins;ibin2++){ 
			for (int ibin3=0;ibin3<nTrkPtBins-1;ibin3++){

				cout<<ibin<<" "<<ibin2<<" "<<ibin3<<endl;

				hJetTrackSignalBackground[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

				hJetTrackME[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackME"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

		//		hJetTrackSignalBackground[ibin][ibin2][ibin3]->Draw("LEGO@");
				cout << "Plot of the Signal Events" << endl;
		//		hJetTrackME[ibin][ibin2][ibin3]->Draw("LEGO@");
				cout << "Plot of the Mixed Events" << endl;


			//	hJetTrackME[ibin][ibin2][ibin3]->ProjectionX("ProEta")->Draw();
				cout << "Eta Projection" <<  endl;

			//	hJetTrackME[ibin][ibin2][ibin3]->ProjectionY("ProPhi")->Draw();
				cout << "Phi Projection" << endl;
				//cout << "Number of Eta bins:" << ProEta->GetNbinsX() << endl;
				//float EtaBins = ProEta->GetNbinsX() ;
				//cout << "Number of Phi bins:" << ProPhi->GetNbinsX() << endl;


			} // end of ibin3
		} //end of ibin2

		cout << "after 2 loops" << endl;
		//TFile *f = TFile::Open("/home/mzakaria3/Documents/Research/Data2011_trkPtCut1_part1.root");
		//TFile *f = TFile::Open("Data2011_trkPtCut1_All_AJ11.root");
		//cout << "took the file" << endl;

		//hists_hJetTrackMECent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Draw("LEGO2");
		//cout << "Plot of the Mixed Events" << endl;

		//hists_hJetTrackMECent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionX("ProEta")->Draw();
		//cout << "Eta Projection" <<  endl;

		//hists_hJetTrackMECent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("ProPhi")->Draw();
		//cout << "Phi Projection" << endl;

		//cout << "Number of Eta bins:" << ProEta->GetNbinsX() << endl;
		//float EtaBins = ProEta->GetNbinsX() ;
		//cout << "Number of Phi bins:" << ProPhi->GetNbinsX() << endl;

		//Normalizing

		for(int ibin2=0; ibin2<nPtBins;ibin2++){
			for (int ibin3=0; ibin3< 4; ibin3++){
		//		printf("the line number is: %d \n \n", __LINE__);

				cout << "hJetTrackSignalBackground " << ibin << " " << ibin2 << " " << ibin3 << endl;
			//	hJetTrackSignalBackground[ibin][ibin2][ibin3]->Draw("LEGO2");

				cout << "Fitting the central part of eta" << endl;

				fout_inc->cd();
				me_proj_canvas->cd(4*(ibin3+1)-ibin);

				yield_inc[ibin][ibin2][ibin3] =  (TH2D*) hJetTrackSignalBackground[ibin][ibin2][ibin3]->Clone((TString)("Yield_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));
				//ProEta->Fit("pol0","","",-0.2,0.2);
				//float s = pol0->GetParameter(0);  // need to make this an array
				me_proj[ibin][ibin2][ibin3] = hJetTrackME[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_temp_%d%d%d",ibin,ibin2,ibin3),1,100);
				//me_proj[ibin][ibin2][ibin3]->Scale (1./100);
				me_proj[ibin][ibin2][ibin3] -> Fit("pol0","","",-.2,.2);
				me00 = 	me_proj[ibin][ibin2][ibin3]->GetFunction("pol0")->GetParameter(0);
				//me_proj[ibin][ibin2][ibin3]->Scale (100./me00);
				hJetTrackME[ibin][ibin2][ibin3]->Scale (100./me00);
				cout<< "me00 = " << me00 <<endl;

				//me_proj[ibin][ibin2][ibin3]->SetAxisRange(-.299,.299);
				//me_proj[ibin][ibin2][ibin3]->Draw();

				//if(ibin<3){
				//	TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
				//	centtex->SetNDC();
				//	centtex->Draw();
				//	TLatex *pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
				//	pttex->SetNDC();
				//	pttex->Draw();
				//}
				//if(ibin==3){
				//	TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
				//	centtex->SetNDC();
				//	centtex->Draw();
				//	TLatex *pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
				//	pttex->SetNDC();
				//	pttex->Draw();
				//}
				yield_inc[ibin][ibin2][ibin3]->Divide(hJetTrackME[ibin][ibin2][ibin3]);
				//yield_inc[ibin][ibin2][ibin3]->Scale(me00);
			//	yield_inc[ibin][ibin2][ibin3]->Draw("LEGO2");

				hJetTrackSignalBackground[ibin][ibin2][ibin3]->Write();
				yield_inc[ibin][ibin2][ibin3]->Write();
				hJetTrackME[ibin][ibin2][ibin3]->Write();


				//cout << "The value we need to normalize for is = " << s << endl;

				//cout << "Normalizing the Mixed Events" << endl;
				//hJetTrackME[ibin][ibin2][ibin3]->Scale(EtaBins/s);
				//	hJetTrackME[ibin][ibin2][ibin3]->Scale(1./me00);

				//TH2D* ME = hJetTrackME[ibin][ibin2][ibin3]->Clone("ME"); //make is an array later

				//hJetTrackSignalBackground[ibin][ibin2][ibin3]->Divide(ME); //done below

				//hJetTrackSignalBackground[ibin][ibin2][ibin3]->Draw("Lego2");
				///hJetTrackSignalBackground[ibin][ibin2][ibin3]->ProjectionY("ResultPhi",20, 80)->Draw();
				//me_proj[ibin][ibin2][ibin3] = hJetTrackME[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_temp_%d%d%d",ibin,ibin2,ibin3),1,100);

				yield_inc[ibin][ibin2][ibin3]->ProjectionY("ResultPhi",20, 80);
				yield_inc_SBL[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionY(Form("ResultSBL",ibin,ibin2,ibin3),21,35 );
				yield_inc_SBR[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionY(Form("ResultSBR",ibin,ibin2,ibin3),66,80 );

				yield_inc_SBL[ibin][ibin2][ibin3]->Add(yield_inc_SBR[ibin][ibin2][ibin3]);
				//yield_inc_SBL[ibin][ibin2][ibin3]->Draw("LEGO2");

				yield_inc_Signal[ibin][ibin2][ibin3]= yield_inc[ibin][ibin2][ibin3]->ProjectionY(Form("me_proj_temp_%d%d%d",ibin,ibin2,ibin3),36,65);
				//yield_inc[ibin][ibin2][ibin3]->ProjectionY("yield_inc_Signal",36,65 )->Draw();
				yield_inc_Signal[ibin][ibin2][ibin3]->Add(yield_inc_SBL[ibin][ibin2][ibin3],-1);
				yield_inc_Signal[ibin][ibin2][ibin3]->Draw();
				yield_inc_Signal[ibin][ibin2][ibin3]->Write();

				hJetTrackSignalBackground[ibin][ibin2][ibin3]->Write();
				yield_inc[ibin][ibin2][ibin3]->Write();
				hJetTrackME[ibin][ibin2][ibin3]->Write();
			
				if(ibin<3){
					TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
					centtex->SetNDC();
					centtex->Draw();
					TLatex *pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
					pttex->SetNDC();
					pttex->Draw();
				}
				if(ibin==3){
					TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
					centtex->SetNDC();
					centtex->Draw();
					TLatex *pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
					pttex->SetNDC();
					pttex->Draw();
				}

				result_proj_canvas->cd(4*(ibin3+1)-ibin);

				int phil = yield_inc[ibin][ibin2][ibin3]->GetYaxis()->FindBin(-.99999);
				int phir = yield_inc[ibin][ibin2][ibin3]->GetYaxis()->FindBin(.99999);

				yield_inc_proj[ibin][ibin2][ibin3]= (TH1D*)hJetTrackSignalBackground[ibin][ibin2][ibin3]->ProjectionX((TString)("yield_inc_proj_"+CBin_strs[ibin]+"_"+TrkPtBin_strs[ibin3]),phil,phir);
				yield_inc_proj_signal[ibin][ibin2][ibin3]= (TH1D*)yield_inc_Signal[ibin][ibin2][ibin3]->Clone((TString)("yield_inc_proj_signal"+CBin_strs[ibin]+"_"+TrkPtBin_strs[ibin3]));
				test_me_inc_proj[ibin][ibin2][ibin3]= (TH1D*)hJetTrackME[ibin][ibin2][ibin3]->ProjectionX((TString)("test_me_inc_proj_"+CBin_strs[ibin]+"_"+TrkPtBin_strs[ibin3]),phil,phir);

				yield_inc_proj[ibin][ibin2][ibin3]->SetMarkerStyle(10);
				yield_inc_proj[ibin][ibin2][ibin3]->SetMarkerSize(1);
				yield_inc_proj[ibin][ibin2][ibin3]->SetMarkerColor(kBlack);
				yield_inc_proj[ibin][ibin2][ibin3]->Draw();

				//test_me_inc_proj[ibin][ibin2][ibin3]->Scale(1/50.);

				test_me_inc_proj[ibin][ibin2][ibin3]->SetMarkerStyle(10);
				test_me_inc_proj[ibin][ibin2][ibin3]->SetMarkerSize(1);
				test_me_inc_proj[ibin][ibin2][ibin3]->SetMarkerColor(kRed);
				test_me_inc_proj[ibin][ibin2][ibin3]->SetLineColor(kRed);
				test_me_inc_proj[ibin][ibin2][ibin3]->Draw("same");

				if(ibin<3){
					centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
					centtex->SetNDC();
					centtex->Draw();
					pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
					pttex->SetNDC();
					pttex->Draw();
				}
				if(ibin==3){
					centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
					centtex->SetNDC();
					centtex->Draw();
					pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
					pttex->SetNDC();
					pttex->Draw();
				}



			}
		}
	}

	me_proj_canvas->SaveAs((TString)("All_ME_Projections_Inclusive_"+datalabel+".png"));
	result_proj_canvas->SaveAs((TString)("Corrected_dEta_Incluisve_Projections_"+datalabel+".png"));

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
