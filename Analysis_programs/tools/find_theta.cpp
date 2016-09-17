//program to find optimal theta to give to L

#include <iostream>
#include <math.h>

using namespace std;

double MH,ML;
int L;

double sqr(double x)
{return x*x;}
  
double momentum(double th)
{return M_PI/L*th;}

double E(double M,double th)
{return 2*asinh(sqrt(3*sqr(sin(momentum(th)/2))+sqr(sinh(M/2))));}

double Q2_HEAVY_REST(double th)
//{return MH*MH+ML*ML-2*MH*E(ML,th);}
{return sqr(MH-E(ML,th))-3*sqr(momentum(th));}

double find_th_null_Q2(double (*Q2)(double))
{
  double thL=0;
  double thR=1.5;
  double Q2L=Q2(thL);
  double Q2R=Q2L;
  
  while(Q2R*Q2L>=0)
    {
      //cout<<"thR: "<<thR<<", Q2R: "<<Q2R<<endl;
      thR*=1.1;
      Q2R=Q2(thR);
    }
  
  do
    {
      double thM=(thL+thR)/2;
      double Q2M=Q2(thM);
      //cout<<"th: "<<thL<<" "<<thM<<" "<<thR<<endl;
      //cout<<"Q2: "<<Q2L<<" "<<Q2M<<" "<<Q2R<<endl;
      if(Q2R*Q2M>=0)
	{
	  Q2R=Q2M;
	  thR=thM;
	}
      else
	{
	  Q2L=Q2M;
	  thL=thM;
	}
    }
  while(fabs(thR-thL)>1.e-10);

  return thR;
}

int main()
{
  cout<<"MH: ";
  cin>>MH;
  
  cout<<"ML: ";
  cin>>ML;
  
  cout<<"L: ";
  cin>>L;
  
  if(MH<ML)
    {
      cerr<<"Error, MH<ML!"<<endl;
      exit(1);
    }
  
  cout<<"Theta: "<<(MH*MH-ML*ML)/(2*MH)*L/sqrt(3)/M_PI<<endl;
  
  cout<<"Theta Isotropic with lattice disp rel: "<<find_th_null_Q2(Q2_HEAVY_REST)<<endl;
  
  double Heavy_rest=(MH*MH-ML*ML)/(2*MH);
  double flm=M_PI/L;
  cout<<"Heavy rest frame: |p|="<<Heavy_rest/flm<<"*pi/L"<<endl;
  
  double Light_rest=(MH*MH-ML*ML)/(2*ML);
  cout<<"Light rest frame: |p|="<<Light_rest/flm<<"*pi/L"<<endl;
  
  double Breit=(MH*MH-ML*ML)/sqrt(8*(MH*MH+ML*ML));
  cout<<"Breit frame: |p|="<<Breit/flm<<"*pi/L"<<"="<<Breit<<endl;
  
  return 0;
}
