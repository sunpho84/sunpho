#include "include.h"

//////////////////////////////////// to find phys point /////////////////////////////

const int nbeta=5;
const double a[nbeta]={0.63401,0.50203,0.418274,0.358883,0.290355};

vector<tuple<int,double,double,double>> data;

void ch2(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;

  cout<<"-------------------"<<endl;
  
  for(auto &it : data)
    {
      int &ib=get<0>(it);
      double &L=get<1>(it);
      double &ynum=get<2>(it);
      double &err=get<3>(it);
      
      double yteo=p[2*ib+0]+p[2*ib+1]*exp(-p[2*nbeta]*L);
      double contr=sqr((yteo-ynum)/err);
      
      cout<<ib<<" "<<1/L<<" "<<yteo<<" "<<p[2*ib+0]<<" "<<p[2*ib+1]<<" "<<err<<" "<<contr<<endl;
      
      ch+=contr;
    }
}

vector<double> fit()
{
  vector<double> out(2*nbeta+1);
  
  TMinuit minu;
  //minu.SetPrintLevel(-1);
  minu.SetFCN(ch2);
  //set tolerance
  double tol=1e-16;
  int iflag;
  //minu.SetMaxIterations(10000000);
  //minu.mnexcm("SET ERR",&tol,1,iflag);
  
  for(int ib=0;ib<nbeta;ib++)
    {
      minu.DefineParameter(2*ib+0,"infi",1,0.001,0,0);
      minu.DefineParameter(2*ib+1,"coef",1,0.001,0,0);
    }
  minu.DefineParameter(2*nbeta,"M",0,0.001,0,0);
  
  minu.Migrad();
  
  double dum;
  for(int ib=0;ib<=2*nbeta;ib++) minu.GetParameter(ib,out[ib],dum);
  
  return out;
}

int main(int narg,char **arg)
{
  int ib;
  double Lfra,me,er;
  ifstream in("data");
  if(!in.good()) crash("unable to open \"data\"");
  while(in>>ib>>Lfra>>me>>er) data.push_back(make_tuple(ib,Lfra*a[ib],me,er));
  
  vector<double> extr=fit();
  
  ofstream out("extr.xmg");
  out<<"@type xydy"<<endl;
  for(int ib=0;ib<nbeta;ib++)
    {
      for(auto &it : data) if(get<0>(it)==ib) out<<get<1>(it)<<" "<<get<2>(it)<<" "<<get<3>(it)<<endl;
      out<<"&"<<endl;
    }
  out<<"@type xy"<<endl;
  for(int ib=0;ib<nbeta;ib++)
    {
      for(double x=5;x<50;x+=1) out<<x<<" "<<extr[2*ib]+extr[2*ib+1]*exp(-extr[2*nbeta]*x)<<endl;
      out<<"&"<<endl;
    }
  
  return 0;
}

