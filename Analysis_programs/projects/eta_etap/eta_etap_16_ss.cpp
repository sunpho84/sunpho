#include "include.h"
#include <TVectorD.h>

//#define TWO_PTS_FIT two_pts_migrad_fit
#define TWO_PTS_FIT two_pts_fit

int njacks=25;
int T=32,TH=T/2;

int tmin_ETA=11,tmax_ETA=14;
int tmin_ETAP=5,tmax_ETAP=6;
int tmin_two_states_fit=4,tmax_two_states_fit=14;

double a=1/1.75;

//fit the mass and the matrix element in SS and SL combo
double *c_two_states_fit,*e_two_states_fit;

double fun_two_states_migrad_fit(double Z2ETA,double META,double Z2ETAP,double METAP,double t)
{
  double c1=Z2ETA*exp(-META*TH)*cosh(META*(TH-t))/META;
  double c2=Z2ETAP*exp(-METAP*TH)*cosh(METAP*(TH-t))/METAP;
  
  return c1+c2;
}

void ch2_two_states_migrad_fit(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  double Z2ETA=p[0];
  double META=p[1];
  double Z2ETAP=p[2];
  double METAP=p[3];
  
  for(int t=tmin_two_states_fit;t<=min(tmax_two_states_fit,TH);t++)
    {
      double num=c_two_states_fit[t];
      double teo=fun_two_states_migrad_fit(Z2ETA,META,Z2ETAP,METAP,t);
      double diff=num-teo;
      double err=e_two_states_fit[t];
      double cont=sqr(diff/err);
      ch+=cont;
      //if(flag==3) cout<<" t="<<t<<", diff=("<<num<<"-"<<teo<<")="<<diff<<" err="<<err<<" cont="<<cont<<endl;
    }
}

jvec load_ss(int isme_so,int isme_si)
{
  //Connected ss
  jvec Css=jvec_load(combine("data/%02d_%02d/SS_conn",isme_so,isme_si).c_str(),T,njacks,0).simmetrized(1);
  ofstream out_Css(combine("plots/%02d_%02d/Css.xmg",isme_so,isme_si).c_str());
  out_Css<<"@type xydy"<<endl;
  out_Css<<Css<<endl;
  cout<<"SS conn: "<<constant_fit(effective_mass(Css),9,15,combine("plots/%02d_%02d/Css_eff.xmg",isme_so,isme_si).c_str())<<endl;
  
  //Disconnected ss
  jvec Dss=jvec_load(combine("data/%02d_%02d/SS_disconn",isme_so,isme_si).c_str(),T,njacks,0).simmetrized(1);
  ofstream out_Dss(combine("plots/%02d_%02d/Dss.xmg",isme_so,isme_si).c_str());
  out_Dss<<"@type xydy"<<endl;
  out_Dss<<Dss<<"&"<<endl;
  constant_fit(effective_mass(Dss),9,15,combine("plots/%02d_%02d/Dss_eff.xmg",isme_so,isme_si).c_str());
  
  //total ss
  jvec ss=Css+Dss;
  ofstream out_ss(combine("plots/%02d_%02d/ss.xmg",isme_so,isme_si).c_str());
  out_ss<<"@type xydy"<<endl;
  out_ss<<ss<<"&"<<endl;
  jack Mss=constant_fit(effective_mass(ss),6,9,combine("plots/%02d_%02d/ss_eff.xmg",isme_so,isme_si).c_str());
  cout<<"SS: "<<smart_print(Mss)<<endl;
  
  jack Z2ss;
  two_pts_fit(Mss,Z2ss,ss,6,9);
  jvec ss_exc=ss;
  for(int t=0;t<=TH;t++) ss_exc[t]-=Z2ss*exp(-Mss*TH)*cosh(Mss*(TH-t))/Mss;
  ofstream out_ss_exc(combine("plots/%02d_%02d/ss_exc.xmg",isme_so,isme_si).c_str());
  out_ss_exc<<"@type xydy"<<endl;
  out_ss_exc<<effective_mass(ss_exc)<<"&"<<endl;
  
  return ss;
}

int main(int narg,char **arg)
{
  jvec Css=jvec_load("data_only_strange/00_00/SS_conn",T,njacks,0).simmetrized(1);
  ofstream out_Css("plots_only_strange/00_00/Css.xmg");
  out_Css<<"@type xydy"<<endl;
  out_Css<<Css<<endl;

  //Disconnected ss
  jvec Dss=jvec_load("data_only_strange/00_00/SS_disconn",T,njacks,0).simmetrized(1);
  ofstream out_Dss("plots_only_strange/00_00/Dss.xmg");
  out_Dss<<"@type xydy"<<endl;
  out_Dss<<Dss<<"&"<<endl;

  ofstream out_ss("plots_only_strange/00_00/ss.xmg");
  out_ss<<"@type xydy"<<endl;
  out_ss<<Css+Dss<<"&"<<endl;
  
  return 0;
  
  jvec ss[16];
  for(int isme_so=0;isme_so<4;isme_so++)
    for(int isme_si=0;isme_si<4;isme_si++)
      ss[isme_so*4+isme_si]=load_ss(8*isme_so,8*isme_si);

  //symmetrize
  for(int isme_so=0;isme_so<4;isme_so++)
      for(int isme_si=0;isme_si<isme_so;isme_si++)
	ss[isme_so*4+isme_si]=ss[isme_si*4+isme_so];

  //cout
  for(int isme_so=0;isme_so<4;isme_so++)
    {
      for(int isme_si=0;isme_si<4;isme_si++) cout<<smart_print(ss[isme_so*4+isme_si][3])<<" ";
      cout<<endl;
    }

  int t0=2,nsme=2,sme[]={2,3};
  gevp_pars_t gevp(nsme,120,T/2,t0);
  for(int isme1=0;isme1<nsme;isme1++)
    for(int isme2=0;isme2<nsme;isme2++)
      gevp.data[isme1*nsme+isme2]=ss[sme[isme1]*4+sme[isme2]];
  
  ofstream gevp_in("plots/gevp_in.xmg");
  gevp_in<<"@type xydy"<<endl;
  for(int isme=0;isme<nsme;isme++)
    {
      gevp_in<<effective_mass(gevp.data[isme])<<endl;
      gevp_in<<"&"<<endl;
    }
  
  gevp.check_norm(4);
  gevp.check_singularity(4);
  gevp.gevp();
  gevp.reorder_eig();
  
  ofstream gevp_out("plots/gevp_out.xmg");
  gevp_out<<"@type xydy"<<endl;
  for(int isme=0;isme<nsme;isme++)
    {
      gevp_out<<effective_mass(gevp.eig_va[isme])<<endl;
      gevp_out<<"&"<<endl;
    }
  
  return 0;
}
