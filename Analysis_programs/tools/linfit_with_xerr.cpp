#include "../src/include.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

vector<jack> y;
vector<double> x,dx,dy;
const int njacks=1000;
int ijack_fit;

template <class _T1,class _T2> _T1 fun_fit(_T1 q,_T1 m,_T2 x)
{
  return q+m*x;
}

void ch2_fit(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  double q=p[0];
  double m=p[1];
  
  for(size_t i=0;i<x.size();i++)
    ch+=sqr(y[i][ijack_fit]-fun_fit(q,m,x[i]))/(sqr(dx[i]*m)+sqr(dy[i]));
}

void fit(jack &q,jack &m)
{
  TMinuit minu;
  minu.SetPrintLevel(-1);
  minu.SetFCN(ch2_fit);

  minu.DefineParameter(0,"q",1,0.001,0,0);
  minu.DefineParameter(1,"m",1,0.001,0,0);
  
  for(ijack_fit=0;ijack_fit<=njacks;ijack_fit++)
    {
      minu.Migrad();
      double dum;
      minu.GetParameter(0,q[ijack_fit],dum);
      minu.GetParameter(1,m[ijack_fit],dum);
    }
  
  ijack_fit=njacks;
  double ch2,grad[2],par[2]={q[njacks],m[njacks]};
  minu.Eval(2,grad,ch2,par,3);
  cout<<"ch2: "<<ch2<<endl;
}

int main(int narg,char **arg)
{
  if(narg<2) crash("use %s filein [fileout]",arg[0]);
    
  ifstream in(arg[1]);
  double _x,_dx,_y,_dy;
  int seed=235413264;
  double xmax=0;
  int i=0;
  while(in>>_x>>_dx>>_y>>_dy)
    {
      x.push_back(_x);
      dx.push_back(_dx);
      y.push_back(fill_gauss(_y,_dy,seed,njacks));
      dy.push_back(_dy);
      xmax=std::max(xmax,_x);
      cout<<x[i]<<" "<<dx[i]<<" "<<y[i]<<endl;
      i++;
    }
  
  cout<<"N points: "<<x.size()<<endl;
  
  jvec pars_fit(2,njacks);
  fit(pars_fit[0],pars_fit[1]);
  
  //plot
  if(narg>2)
    {
      ofstream out(arg[2]);
      out<<"@type xy"<<endl<<
	write_poly_with_error(pars_fit,0,xmax*1.1)<<endl
	 <<"&"<<endl
	 <<"@s1 line color 1"<<endl
	 <<"@type xydxdy"<<endl
	 <<"@s2 linestyle 0"<<endl
	 <<"@s2 symbol 1"<<endl
	 <<"@s2 symbol color 2"<<endl
	 <<"@s2 errorbar color 2"<<endl
	;
	
      for(size_t i=0;i<x.size();i++)
	out<<x[i]<<" "<<y[i].med()<<" "<<dx[i]<<" "<<y[i].err()<<endl;
    }
  
  cout<<pars_fit<<endl;
  
  return 0;
}
