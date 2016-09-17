#include <TApplication.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <cstdlib>

#include <fstream>
#include <iostream>

using namespace std;

const double nc=3;
const double nf=3;
const double pref=1/(16*M_PI*M_PI);
const double gamma_0=pref*(11*nc-2*nf)/3;
const double gamma_1=pref*pref*(34*nc*nc*nc-10*nc*nc*nf-3*nf*(nc*nc-1))/(3*nc);
double a(double *x,double *p)
{
  double beta=x[0];
  double lambda=p[0];
  double g02=2*nc/beta;
  return 1/lambda*pow(g02*gamma_0,-gamma_1/(2*gamma_0*gamma_0))*exp(-1/(2*gamma_0*g02));
}

int main(int narg,char **arg)
{
  int i;
  double x,y,ey;
  ifstream file;
  TApplication myapp("App",NULL,NULL);
  TCanvas tela;
  TGraphErrors graf;
  TF1 fun("fun",a,3.6,4.2,1);
  
  if(narg<2)
    {
      cerr<<"usa: "<<arg[0]<<" file"<<endl;
      cerr<<"il file è del tipo:"<<endl;
      cerr<<"x1 y1 ey1"<<endl;
      cerr<<"x2 y2 ey2"<<endl;
      cerr<<"etc"<<endl;
      exit(0);
    }
  
  file.open(arg[1]);
  
  if(file.good()==0)
    {
      cout<<"File: "<<arg[1]<<" inesistente!"<<endl;
      exit(0);
    }
  
  i=0;
  while(file>>x>>y>>ey)
    {
      cout<<x<<" "<<y<<" "<<ey<<endl;
      graf.SetPoint(i,x,y);
      graf.SetPointError(i,0,ey);
      i++;
    }
  
  fun.SetParameter(0,1.013);
  //fun.FixParameter(0,0.013);
  
  graf.Draw("AP");
  //fun.Draw("L");
  graf.Fit("fun");
  
  tela.Modified();
  tela.Update();
  
  myapp.Run(true);
  
  return 0;
}
