
{

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <stdio.h>

#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>


using namespace std;

void FitAHisto();

TH1F *fit1d(TH1F *hist, int i, const char *func, double lolim, double hilim)
{
  hist->SetLineColor(i+1);
  hist->SetMarkerColor(i+1);

  char fitname[256];
  strcpy (fitname,hist->GetName());
  strcat (fitname,"Fitted");
 
  TF1 *fitfunc = new TF1 (fitname,func, lolim, hilim);
  fitfunc->SetLineColor(i+1);
  hist->Fit(fitname,"0EMR"); // "0" = do NOT automatically draw "hist"
  hist->GetFunction(fitname)->ResetBit(TF1::kNotDraw); // make "fitname" visible (= 1<<9)
 
  return hist;

}

void FitAHisto(){
  TH1F *hh=new TH1F("hh","hh",100,-1.,1.);
  TH1F *hh2=new TH1F("hh2","hh2",100,-1.,1.);
  hh->FillRandom("gaus",100000);
  hh2->FillRandom("gaus",90000);
 
  TCanvas *cc =new TCanvas("cc","cc",500,500);
  cc->Divide(2,3);
  cc->cd(1);
  fit1d(hh, 1, "gaus", -1.,1.)->Draw();
  fit1d(hh2, 6, "gaus", -1.,1.)->Draw("same");
  cc->cd(2);
  hh->Draw();
  hh2->Draw("same");
  cc->cd(3);
  fit1d(hh, 1, "gaus", -1.,1.)->Draw();
  cc->cd(4);
  fit1d(hh2, 6, "gaus", -1.,1.)->Draw();

  TH1F *hhh=new TH1F("hhh","hhh",100,-1.,1.);
  TH1F *hhh2=new TH1F("hhh2","hhh2",100,-1.,1.);
  hhh->SetLineColor(2);
  hhh2->SetLineColor(7);
  hhh->FillRandom("gaus",100000);
  hhh2->FillRandom("gaus",90000);


  TF1 *fff = new TF1("fff","gaus",-1.,1.);
  TF1 *fff2 = new TF1("fff2","gaus",-1.,1.);
  fff->SetLineColor(2);
  fff2->SetLineColor(7);

  cc->cd(5);
  hhh->Fit("fff","emr"); // automatically draws "hhh"
  hhh2->Fit("fff2","emr","same"); // automatically draws "hhh2"
  cc->cd(6);
  hhh->Draw();
  hhh2->Draw("same");

}

}
