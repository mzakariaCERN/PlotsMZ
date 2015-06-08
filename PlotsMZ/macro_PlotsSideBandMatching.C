void macro_PlotsSideBandMatching()
{
int NumberDijets = 0;
double dPhiBinWidth= 0; 

//TFile *f = TFile::Open("/home/mzakaria3/Documents/Research/MyGitProject/PlotsMZ/Data2011_trkPtCut1_All_AJ11.root");
//TFile *f = TFile::Open("/home/mzakaria3/Documents/Research/MyGitProject/PlotsMZ/Data2011_trkPtCut1_All_AJ11.root");
//TFile *f = TFile::Open("/home/mzakaria3/Documents/Research/MyGitProject/PlotsMZ/Data2011_trkPtCut1_All_Aj11_20150319.root");
TFile *f = TFile::Open("/home/mzakaria3/Documents/Research/MyGitProject/PlotsMZ/Data2011_trkPtCut1_All_Aj11_20150318.root");
//TFile *f = TFile::Open("PbPb_AjBins_Cent0_Cent10_Pt1_Pt2.root");
//TFile *f = TFile::Open("/home/mzakaria3/Documents/Research/MyGitProject/PlotsMZ/PbPb_AjBins_Cent0_Cent10_Pt1_Pt2.root");
//TFile *f = TFile::Open("/home/mzakaria3/Documents/Research/MyGitProject/PlotsMZ/Data2011_All_AJ33.root");
cout << "took the file" << endl;

//cout << "hists_AjCent50_Cent100->GetEntries() = " << hists_AjCent50_Cent100->GetEntries()  << endl;

//NumberDijets = hists_AjCent0_Cent10->GetEntries() + hists_AjCent10_Cent30->GetEntries() + hists_AjCent30_Cent50->GetEntries() + hists_AjCent50_Cent100->GetEntries();
NumberDijets = hists_only_leadingjets_corrpTCent0_Cent10_Pt100_Pt300->GetEntries();
cout << "NumberDijets using corrpT= " <<  NumberDijets << endl;

cout << "Number of Dijets = " << NumberDijets << endl;
//hists_hJetTrackMELeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Draw("LEGO2");
hists_hJetTrackMELeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Draw("LEGO2");
cout << "Plot of the Mixed Events" << endl;

//hists_hJetTrackMELeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionX("ProEta")->Draw();
hists_hJetTrackMELeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionX("ProEta")->Draw();
cout << "Eta Projection" <<  endl;

//hists_hJetTrackMELeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("ProPhi")->Draw();
hists_hJetTrackMELeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("ProPhi")->Draw();
cout << "Phi Projection" << endl;
 
cout << "Number of Eta bins:" << ProEta->GetNbinsX() << endl;
float EtaBins = ProEta->GetNbinsX();
float dEtaBinWidth = ProEta->GetBinWidth(1) ;
cout << "Number of Phi bins:" << ProPhi->GetNbinsX() << endl;
dPhiBinWidth =  ProPhi->GetBinWidth(1);
cout << "dPhi Bin width = " << dPhiBinWidth << endl;

cout << "Fitting the central part of eta" << endl;
ProEta->Fit("pol0","","",-0.2,0.2);
float s = pol0->GetParameter(0); 
cout << "The value we need to normalize for is = " << s << endl;
//cout << "bin(0,0) before = " << hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->GetBinContent(50,37) << endl;
cout << "bin(0,0) before = " << hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->GetBinContent(50,37) << endl;

cout << "Normalizing the Mixed Events" << endl;
//hists_hJetTrackMELeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Scale(EtaBins/s);
hists_hJetTrackMELeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Scale(EtaBins/s);
//cout << "bin(0,0) after = " << hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->GetBinContent(50,37) << endl;
cout << "bin(0,0) after = " << hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->GetBinContent(50,37) << endl;

//hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Draw("LEGO2");
hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Draw("LEGO2");
cout << "Cloning Mixed Event to devide by it" << endl;
//TH2D* ME = hists_hJetTrackMELeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Clone("ME");
TH2D* ME = hists_hJetTrackMELeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Clone("ME");

cout << "Divide Signal by ME ...." <<  endl;
hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Divide(ME);

//cout << "bin(0,0) after = " << hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->GetBinContent(50,37) << endl;
cout << "bin(0,0) after = " << hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->GetBinContent(50,37) << endl;

//hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Draw("Lego2");
hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Draw("Lego2");

//hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("ResultPhi",50, 25)->Draw();

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
//side and check
hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionX("rawEta",1,50)->Draw();
//return 0;    //OLGA
//hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("SBL",21,35 )->Draw();
hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("SBL",21,35 )->Draw();
//hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("SBR",66,80 )->Draw();
hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("SBR",66,80 )->Draw("same");
SBL->Add(SBR);
SBL->Draw();

//hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("Signal",36,65 )->Draw();
//Signal->Add(SBL,-1);
// Signal->Draw();
hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionX("Check",1, 50)->Draw();


TH1D *Signal = (TH1D*)hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("Signal",36, 65);
Signal->Draw();
SBL->Draw("SAME");

Signal->Rebin(4);
Signal->SetMarkerColor(kRed);
Signal->SetMarkerStyle(20);
Signal->Draw();
SBL->SetMarkerStyle(20);

SBL->Rebin(4);
SBL->Draw("SAME");
//hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2>ProjectionY("Signal",36,65 )->Draw();
//hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->Scale(1.0/(3 *NumberDijets * dPhiBinWidth ));
/////Signal->Scale(1.0/(3 *NumberDijets * dPhiBinWidth ));
//hists_hJetTrackSignalBackgroundLeadingCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2->ProjectionY("Signal",36,65 )->Draw();
//Signal->Draw();

                   }                                           
