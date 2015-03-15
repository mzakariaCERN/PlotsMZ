int macro_Plots(int i=0)
{

      const int nCBins = 4;
        const int nPtBins = 1;
          const int nTrkPtBins = 5;

            float me00;
              int me00bin;

                float PtBins[nPtBins+1] = {100, 300};
                  TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};
                    

                    float CBins[nCBins+1] = {0, 20, 60, 100, 200};
                      TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
                        TString CBin_labels[nCBins] = {"Cent. 0-10%","Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};

                          float TrkPtBins[nTrkPtBins+1] = {1, 2, 3, 4, 8, 999};
                            TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt999" };
                              TString TrkPtBin_labels[nTrkPtBins] = {"1<pT<2","2<pT<3","3<pT<4","4<pT<8","pT>8"}
    cout << "the default parameter = " <<  i <<  endl;



TString ref_name;
 
  gStyle->SetOptStat(0);  
    gStyle->SetPadBottomMargin(0.15);
      gStyle->SetPadTopMargin   (0.05);
        gStyle->SetPadLeftMargin  (0.15);
          gStyle->SetPadRightMargin (0.05);

//TFile *f = TFile::Open("/home/mzakaria3/Documents/Research/Data2011_trkPtCut1_part1.root");
TFile *f = TFile::Open("Data2011_trkPtCut1_All_AJ11.root");
cout << "took the file" << endl;

hists_hJetTrackMECent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Draw("LEGO2");
cout << "Plot of the Mixed Events" << endl;

hists_hJetTrackMECent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionX("ProEta")->Draw();
cout << "Eta Projection" <<  endl;

hists_hJetTrackMECent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("ProPhi")->Draw();
cout << "Phi Projection" << endl;

cout << "Number of Eta bins:" << ProEta->GetNbinsX() << endl;
float EtaBins = ProEta->GetNbinsX() ;
cout << "Number of Phi bins:" << ProPhi->GetNbinsX() << endl;

cout << "Fitting the central part of eta" << endl;
ProEta->Fit("pol0","","",-0.2,0.2);
float s = pol0->GetParameter(0); 
cout << "The value we need to normalize for is = " << s << endl;

cout << "Normalizing the Mixed Events" << endl;
hists_hJetTrackMECent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Scale(EtaBins/s);

hists_hJetTrackSignalBackgroundCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Draw("LEGO2");
cout << "Cloning Mixed Event to devide by it" << endl;
TH2D* ME = hists_hJetTrackMECent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Clone("ME");

cout << "Divide Signal by ME ...." <<  endl;
hists_hJetTrackSignalBackgroundCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Divide(ME);
hists_hJetTrackSignalBackgroundCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Draw("Lego2");

hists_hJetTrackSignalBackgroundCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("ResultPhi",20, 80)->Draw();

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

hists_hJetTrackSignalBackgroundCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("SBL",21,35 )->Draw();
hists_hJetTrackSignalBackgroundCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("SBR",66,80 )->Draw();
SBL->Add(SBR);
SBL->Draw();

hists_hJetTrackSignalBackgroundCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("Signal",36,65 )->Draw();
Signal->Add(SBL,-1);
 Signal->Draw();

hists_hJetTrackSignalBackgroundCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("Signal",36,65 )->Draw();

TCanvas *me_proj_canvas = new TCanvas("me_proj_canvas","",0,0,1500,1600);
  me_proj_canvas->Divide(4,4,0.0000,0.0000);


return i;
                   }                                           
