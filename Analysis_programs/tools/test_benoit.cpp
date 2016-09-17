#include "include.h"

int njack=20,T=48,TH=24;

//impulse
double momentum(double th)
{return M_PI/TH*th;}

double E(double M,double th)
{return 2*asinh(sqrt(3*sqr(sin(momentum(th)/2))+sqr(sinh(M/2))));}

jack E(jack &M,double th)
{
  jack o(njack);
  
  for(int ijack=0;ijack<=njack;ijack++) o.data[ijack]=E(M[ijack],th);
  
  return o;
}

jvec load(const char *path)
{
  jvec corr(T,njack);
  
  ifstream in(path);
  
  for(int t=0;t<T;t++)
    {
      int i;
      double y,ey;
      
      in>>i>>y>>ey;

      //load
      corr.data[t].data[njack]=0;
      for(int ijack=0;ijack<njack;ijack++)
	{
	  in>>corr.data[t].data[ijack];
	  corr.data[t].data[njack]+=corr.data[t].data[ijack];
	}

      corr.data[t].data[njack]/=njack;
    }
  
  return corr.simmetrized(1);
}

int main()
{
  jvec corr0=load("correlmas1theta00.dat");
  
  double th[4]={0.0,0.1,0.2,0.3};
  
  ofstream out_teo("out_teo.xmg");
  out_teo<<"@type xydy"<<endl;
  
  ofstream out_fit("out_fit.xmg");
  out_fit<<"@type xydy"<<endl;
  
  ofstream out_bis("out_bis.xmg");
  out_bis<<"@type xydy"<<endl;
  
  jack M[4],deltaM[4];
  for(int ith=0;ith<=3;ith++)
    {
      jvec corri=load(combine("correlmas1theta0%d.dat",ith).c_str());
      
      int tmin=10;
      int tmax=23;
      M[ith]=constant_fit(effective_mass(corri),tmin,tmax,combine("out_%d.xmg",ith).c_str());
      deltaM[ith]=constant_fit(aperiodic_effective_mass(corri/corr0),4,6,combine("delta_%d.xmg",ith).c_str());
      
      double q2=3*sqr(th[ith]*M_PI/TH);
      
      out_teo<<q2<<" "<<E(M[0],th[ith])-M[0]<<endl;
      
      out_fit<<q2+0.00005<<" "<<M[ith]-M[0]<<endl;
      out_bis<<q2+0.00005<<" "<<deltaM[ith]<<endl;
      
      jack Mfit,Z2fit;
      two_pts_migrad_fit(Mfit,Z2fit,corri,tmin,tmin+4);

      cout.precision(16);
      cout<<"M["<<ith<<"]: "<<smart_print(M[ith])<<" "<<smart_print(Mfit)<<endl;
    }
  
  return 0;
}
