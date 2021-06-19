//Macro to plot the histograms from root file. Use cvsq.C while starting root initially.

#ifndef __CINT__
#endif
#include <iostream>
#include "TTree.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"

using namespace std;

void plot_macro2()
{
  for (int i=65;i<=80;i++)
  {
    char* s = (Form("scopeRD359%d.root",i));
    cout<<s<<endl;
    //TCanvas* c = new TCanvas("c","Efficiency vs xmod ymod",800,800) ;      
    TFile *inputfile = new TFile(s,"READ");
    TH1 *heg  = (TH1*)inputfile->Get("linnpxvsxmym");
    heg->Draw("zcol");
    heg->SetMinimum(1.0);
    heg->SetMaximum(1.7);
    c1->SaveAs(Form("linnpxvsxmym_%d.png",i));
    inputfile->Close();
  }
}
