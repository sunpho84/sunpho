#include <TMinuit.h>
#include <TMath.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

using namespace std;

const double dl1b=-0.4,dl2b=4.3,dl3b=2.9,dl4b=4.4;

double G1t(double x)
{
  const int nterm=20;
  const int mul[nterm]={6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24};
  double G1t=0,y;

  for(int nn=0;nn<nterm;nn++)
    {
      y=x*sqrt((double)(nn+1));
      G1t+=4*mul[nn]*TMath::BesselK1(y)/y;
    }

  return G1t;
}      

void volume(double xl,double & g1t,double & g2tm,double & g2tf,double & g2tv)
{
  double coef1=-55.0/18+4*dl1b+8*dl2b/3-2.5*dl3b-2*dl4b;
  double coef2=112.0/9-8*dl1b/3-32*dl2b/3;
  double coef3=-7.0/9+2*dl1b+4*dl2b/3-3*dl4b;
  double g0=2-0.5*M_PI;
  double g1=0.25*M_PI-0.5;
  double g2=0.5-0.125*M_PI;
  double g3=3.0*M_PI/16-0.5;
  double dmn,b0,b2;

  const int nterm=20;
  const int mul[nterm]={6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24};
  int nn;
  double y;

  g1t=g2tm=g2tf=g2tv=0;

  for(nn=0;nn<nterm;nn++)
    {
      y=xl*sqrt((double)(nn+1));
      dmn=mul[nn];
      b2=TMath::BesselK0(y);
      g2tv+=dmn*b2;
      b0=TMath::BesselK1(y);
      b2=(b2+2*b0/y)/y;
      g1t+=4*dmn*b0/y;
      g2tm+=dmn*(coef1*b0+coef2*b2+13*g0*b0/3-(40*g0+32*g1+26*g2)*b2/3)/y;
      g2tf+=2*dmn*(coef3*b0+coef2*b2+(8*g0-13*g1)*b0/6-(40*g0-12*g1-8*g2-13*g3)*b2/3)/y;
    }
}

void corr_vf(double &corrvm,double &corrvf,double ampi,double afpi,double lato)
{
  double dn=(4*M_PI*4*M_PI);
      
  double corr0;
  double ampi2=ampi*ampi;
  double xl=lato*ampi;
  double csi=2*ampi2/(dn*afpi*afpi);

  double g1t,g2tm,g2tf,g2tv;

  volume(xl,g1t,g2tm,g2tf,g2tv);
  
  corr0=-2*(dl4b*g1t-0.5*dl3b*g2tv);
  
  corrvm=1+0.25*csi*g1t-csi*csi*(g2tm-0.25*corr0);
  corrvf=1-csi*g1t+csi*csi*(g2tf-corr0);
}

void no_corr_vf(double &corr_mpi,double &corr_fpi,double ampi,double afpi,double lato)
{
  double xlam=ampi*lato;
  double xil=ampi*ampi/pow(4*M_PI*afpi,2);

  corr_mpi=sqrt(1+ xil*G1t(xlam));
  corr_fpi=1-2*xil*G1t(xlam);
}


