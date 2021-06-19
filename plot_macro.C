#include <iostream>
#include <fstream>
#include <cmath>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <string> 
//#include <pngwriter.h>

using namespace std;
std:: string ss;
std:: string first ("scopeRD359");
std:: string last (".root");
double plot_macro()
{
  int n;
  double x,y;
  ifstream infile;
  for (int i=65;i<=80;i++)
    {
      std::string s = std::to_string(i);
      ss = first + s + last;
      infile.open(ss);
      //infile->cd("")
      if(infile.is_open()){
	cout<<ss<<endl;
      } else{
	cout<<"file isn't open"<<endl;
      }
      
      //if (infile == true)
      //{
      //cout <<ss<<endl;
      //}
      //pngwriter png(1024, 768, 0.0, "output.png");
      //png.plot(23, 42, 0.5, 1.0, 0.7);
      //png.close();
     
      infile.close();
    }
 return 0;
}
