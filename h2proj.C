// Example displaying a 2D histogram with its two projections.
// Author: Olivier Couet
{
   TCanvas *c1 = new TCanvas("c1", "c1",900,900);
<<<<<<< HEAD
   gStyle->SetOptStat(0);
=======
   //gStyle->SetOptStat(0);
>>>>>>> ef71881c634705e9a0d81e499b30bc1a5dd6f107
   
   // Create the three pads
   TPad *center_pad = new TPad("center_pad", "center_pad",0.0,0.0,0.6,0.6);
   center_pad->Draw();

   right_pad = new TPad("right_pad", "right_pad",0.55,0.0,1.0,0.6);
   right_pad->Draw();

   top_pad = new TPad("top_pad", "top_pad",0.0,0.55,0.6,1.0);
   top_pad->Draw();

   // Create, fill and project a 2D histogram.
<<<<<<< HEAD
   TH2F *h2 = new TH2F("h2","",40,-20,20,40,-20,20);
   Float_t px, py;
   //for (Int_t i = 0; i < 25000; i++) {
   for (Int_t i = 0; i < 40; i++) {
   for (Int_t j = 0; j < 40; j++) {


      //gRandom->Rannor(px,1);
      //h2->Fill(px,5*py);
int k = 0;
	if (i > 15 && i < 25 && j > 15 && j < 25)
	k = i * j;
	else 
	if (i < 15 || i > 25)
	k = 50;

      h2->SetBinContent(i,j, k);
   }
	}
   TH1D * projh2X = h2->ProjectionX();
   TH1D * projh2Y = h2->ProjectionY();

	projh2X->Scale(1. / 40.);
	projh2Y->Scale(1. / 40.);
=======
   TH2F *h2 = new TH2F("h2","",40,-4,4,40,-20,20);
   Float_t px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      //h2->Fill(px,5*py);
      h2->Fill(1,11);
   }
   TH1D * projh2X = h2->ProjectionX("projh2X",32,33);
   TH1D * projh2Y = h2->ProjectionY("projh2Y",25,26);
>>>>>>> ef71881c634705e9a0d81e499b30bc1a5dd6f107

   // Drawing
   center_pad->cd();
   gStyle->SetPalette(1);
<<<<<<< HEAD
   //h2->Draw("COL");
   h2->Draw("LEGO2");
=======
   h2->Draw("COL");
>>>>>>> ef71881c634705e9a0d81e499b30bc1a5dd6f107

   top_pad->cd();
   projh2X->SetFillColor(kBlue+1);
   projh2X->Draw("bar");

   right_pad->cd();
   projh2Y->SetFillColor(kBlue-2);
<<<<<<< HEAD
	projh2Y->GetXaxis()->SetTitle("ProjY");
=======
>>>>>>> ef71881c634705e9a0d81e499b30bc1a5dd6f107
   projh2Y->Draw("hbar");
   
   c1->cd();
   TLatex *t = new TLatex();
   t->SetTextFont(42);
   t->SetTextSize(0.02);
   t->DrawLatex(0.6,0.88,"This example demonstrate how to display");
   t->DrawLatex(0.6,0.85,"a histogram and its two projections.");
}
