#include <iostream>
#include <fstream>
#include <cmath>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>

using namespace std;
double effvsbias()
{
  int n=14;
  float x[n],y[n];
  ifstream infile;
  infile.open("effvsbias.dat");
  for(int i=1;i<=n;i++)
    {
      infile>>x[i]>>y[i];
      cout<<x[i]<<" "<<y[i]<<endl;
    }
  // plot x,y;
  c1= new TCanvas("c1","Efficiency Vs Bias voltage",800,600);
  TGraph *gr1 = new TGraph (n,x,y); 
  gr1->Draw("AC*");
  c1->Draw();
  infile.close();
  return 0;
}
