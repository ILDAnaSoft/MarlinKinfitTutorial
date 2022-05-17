
#include "TH1F.h"
#include "TFile.h"

 
void hipl3(TString hname, TString href, bool log, TString leg1, TString leg2, const char* FILE1, const char* FILE2) {

  TFile* f1 = new TFile(FILE1);
  TFile* f2 = new TFile(FILE2);
  TString hnamefull = "MyZHllqqFit/" + hname;
  TString hreffull = "MyZHllqqFit/" + href;
  TH1F* h1 = (TH1F *) f1->Get(hnamefull.Data());
  TH1F* h2 = (TH1F *) f2->Get(hnamefull.Data());
  TH1F* h3 = (TH1F *) f1->Get(hreffull.Data());
  
  float ymax = h1->GetMaximum();
  if (h2->GetMaximum() > ymax) ymax = h2->GetMaximum();
  if (h3->GetMaximum() > ymax) ymax = h3->GetMaximum();
  if (log) {
    h1->SetMaximum (ymax * 2); 
    h1->SetMinimum (1); 
  } else {
    h1->SetMaximum (ymax * 1.1); 
    h1->SetMinimum (0); 
  }  
  h1->SetLineColor (1);  
  h1->SetLineStyle (1);  
  h2->SetLineColor (2);  
  h2->SetLineStyle (2);  
  h3->SetLineColor (4);  
  h3->SetLineStyle (3); 
  TString htitle = h1->GetTitle();
  htitle = " ; " + htitle;
  h1->SetTitle (htitle); 

  gStyle->SetOptStat(0);
  
  TCanvas* c1 = new TCanvas("c1","hname",600,400);
  
  if (log) gPad->SetLogy();
  
  h1->Draw();
  h2->Draw("same");
  h3->Draw("same");
  
  TString leg3 = "nofit";
  TLegend* leg = new TLegend(0.6,0.65,0.9,0.9);
  leg->SetHeader(hname); // option "C" allows to center the header
  leg->AddEntry(h1,leg1,"L");
  leg->AddEntry(h2,leg2,"L");
  leg->AddEntry(h3,leg3,"L");
  leg->Draw();
 
  TString outnamepdf = TString(hname) + "_" + leg1 + "_" + leg2 + "_" + leg3 + ".pdf";
  c1->Print(outnamepdf);

}
