#include "include.h"
#include <TVectorD.h>

//#define TWO_PTS_FIT two_pts_migrad_fit
#define TWO_PTS_FIT two_pts_fit

int njacks=18;
int T=32,TH=T/2;
double disc_norm=1;

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

int main(int narg,char **arg)
{
  if(narg<2) crash("use %s n",arg[0]);
  
  //Connected ll
  jvec Cll=jvec_load(combine("data/%s/LL_conn",arg[1]).c_str(),T,njacks,0).simmetrized(1);
  ofstream out_Cll(combine("plots/%s/Cll.xmg",arg[1]).c_str());
  out_Cll<<"@type xydy"<<endl;
  out_Cll<<Cll<<endl;
  jack M_Pi,Z2_Pi;
  two_pts_migrad_fit(M_Pi,Z2_Pi,Cll,9,15,combine("plots/%s/Cll_eff.xmg",arg[1]).c_str());
  
  //compute covariance matrix of effective mass
  char label[3][30]={"correlator","effective mass","Antonin effective mass"};
  for(int iter=0;iter<3;iter++)
    {
      cout<<"Normalized correlation matrix for: "<<label[iter]<<endl;
      jvec te;
      switch(iter)
	{
	case 0:te=Cll;break;
	case 1:te=effective_mass(Cll);break;
	case 2:te=antonin_effective_mass(Cll);break;
	default:crash("what the hell");
	}
      
      //fill the matrix
      TMatrixD d(5,5);
      for(int t=0;t<5;t++)
	for(int s=0;s<5;s++)
	  d[t][s]=cov(te[t+9],te[s+9]);///sqrt(cov(te[t+9],te[t+9])*cov(te[s+9],te[s+9]));
      
      //compute on and off diagonal part
      double tot=0,par=0;
      for(int t=0;t<5;t++)
	{
	  par+=d[t][t]*d[t][t];
	  for(int s=0;s<5;s++) tot+=d[t][s]*d[t][s];
	}
      cout<<(par/tot)<<endl;
    }

  //Disconnected ll
  jvec Dll=jvec_load(combine("data/%s/LL_disconn",arg[1]).c_str(),T,njacks,0).simmetrized(1)*disc_norm;
  ofstream out_Dll(combine("plots/%s/Dll.xmg",arg[1]).c_str());
  out_Dll<<"@type xydy"<<endl;
  out_Dll<<Dll<<endl;
  ofstream out_Dll_eff(combine("plots/%s/Dll_eff.xmg",arg[1]).c_str());
  out_Dll_eff<<"@type xydy"<<endl;
  out_Dll_eff<<effective_mass(Dll)<<"&"<<endl;
  
  jack M_Pi_D,Z2_Pi_D;
  two_pts_migrad_fit(M_Pi_D,Z2_Pi_D,Dll,9,15,combine("plots/%s/Dll_fit.xmg",arg[1]).c_str());
  cout<<"Disc pi: "<<smart_print(M_Pi_D)<<"&"<<endl;
  
  //Connected ls
  jvec Cls=jvec_load(combine("data/%s/LS_conn",arg[1]).c_str(),T,njacks,0).simmetrized(1);
  ofstream out_Cls(combine("plots/%s/Cls.xmg",arg[1]).c_str());
  out_Cls<<"@type xydy"<<endl;
  out_Cls<<Cls<<"&"<<endl;
  ofstream out_Cls_eff(combine("plots/%s/Cls_eff.xmg",arg[1]).c_str());
  out_Cls_eff<<"@type xydy"<<endl;
  out_Cls_eff<<effective_mass(Cls)<<"&"<<endl;
    
  //Disconnected ls
  jvec Dls=jvec_load(combine("data/%s/LS_disconn",arg[1]).c_str(),T,njacks,0).simmetrized(1)*disc_norm;
  ofstream out_Dls(combine("plots/%s/Dls.xmg",arg[1]).c_str());
  out_Dls<<"@type xydy"<<endl;
  out_Dls<<Dls<<endl;
  ofstream out_Dls_eff(combine("plots/%s/Dls_eff.xmg",arg[1]).c_str());
  out_Dls_eff<<"@type xydy"<<endl;
  out_Dls_eff<<effective_mass(Dls)<<"&"<<endl;

  //Disconnected sl
  jvec Dsl=jvec_load(combine("data/%s/SL_disconn",arg[1]).c_str(),T,njacks,0).simmetrized(1)*disc_norm;
  ofstream out_Dsl(combine("plots/%s/Dsl.xmg",arg[1]).c_str());
  out_Dsl<<"@type xydy"<<endl;
  out_Dsl<<Dsl<<endl;
  ofstream out_Dsl_eff(combine("plots/%s/Dsl_eff.xmg",arg[1]).c_str());
  out_Dsl_eff<<"@type xydy"<<endl;
  out_Dsl_eff<<effective_mass(Dsl)<<"&"<<endl;

  //Connected ss
  jvec Css=jvec_load(combine("data/%s/SS_conn",arg[1]).c_str(),T,njacks,0).simmetrized(1);
  ofstream out_Css(combine("plots/%s/Css.xmg",arg[1]).c_str());
  out_Css<<"@type xydy"<<endl;
  out_Css<<Css<<endl;
  ofstream out_Css_eff(combine("plots/%s/Css_eff.xmg",arg[1]).c_str());
  out_Css_eff<<"@type xydy"<<endl;
  out_Css_eff<<effective_mass(Css)<<"&"<<endl;
  
  //Disconnected ss
  jvec Dss=jvec_load(combine("data/%s/SS_disconn",arg[1]).c_str(),T,njacks,0).simmetrized(1)*disc_norm;
  ofstream out_Dss(combine("plots/%s/Dss.xmg",arg[1]).c_str());
  out_Dss<<"@type xydy"<<endl;
  out_Dss<<Dss<<"&"<<endl;
  
  //total ll
  jvec ll=Cll-2*Dll;
  ofstream out_ll(combine("plots/%s/ll.xmg",arg[1]).c_str());
  out_ll<<"@type xydy"<<endl;
  out_ll<<ll<<"&"<<endl;
  ofstream out_ll_eff(combine("plots/%s/ll_eff.xmg",arg[1]).c_str());
  out_ll_eff<<"@type xydy"<<endl;
  out_ll_eff<<effective_mass(ll)<<"&"<<endl;
  
  //total ls
  jvec ls=-sqrt(2)*Dls;
  ofstream out_ls(combine("plots/%s/ls.xmg",arg[1]).c_str());
  out_ls<<"@type xydy"<<endl;
  out_ls<<ls<<"&"<<endl;
  ofstream out_ls_eff(combine("plots/%s/ls_eff.xmg",arg[1]).c_str());
  out_ls_eff<<"@type xydy"<<endl;
  out_ls_eff<<effective_mass(ls)<<"&"<<endl;
  
  //total sl
  jvec sl=-sqrt(2)*Dsl;
  ofstream out_sl(combine("plots/%s/sl.xmg",arg[1]).c_str());
  out_sl<<"@type xydy"<<endl;
  out_sl<<sl<<"&"<<endl;
  ofstream out_sl_eff(combine("plots/%s/sl_eff.xmg",arg[1]).c_str());
  out_sl_eff<<"@type xydy"<<endl;
  out_sl_eff<<effective_mass(sl)<<"&"<<endl;
  
  //total ss
  jvec ss=Css-Dss;
  ofstream out_ss(combine("plots/%s/ss.xmg",arg[1]).c_str());
  out_ss<<"@type xydy"<<endl;
  out_ss<<ss<<"&"<<endl;
  ofstream out_ss_eff(combine("plots/%s/ss_eff.xmg",arg[1]).c_str());
  out_ss_eff<<"@type xydy"<<endl;
  out_ss_eff<<effective_mass(ss)<<"&"<<endl;
  
  int t0=2;
  gevp_pars_t gevp(2,144,T/2,t0);
  gevp.data[0]=ll;
  gevp.data[1]=gevp.data[2]=ls;
  gevp.data[3]=ss;
  
  cout<<smart_print(ll[t0])<<" "<<smart_print(ls[t0])<<endl;
  cout<<smart_print(ls[t0])<<" "<<smart_print(ss[t0])<<endl;
  
  cout<<smart_print(ll[t0]/ll[t0])<<" "<<smart_print(ls[t0]/sqrt(ll[t0]*ss[t0]))<<endl;
  cout<<smart_print(ls[t0]/sqrt(ll[t0]*ss[t0]))<<" "<<smart_print(ss[t0]/ss[t0])<<endl;
  
  //gevp.check_orthogonality();
  gevp.gevp();
  gevp.reorder_eig();
  
  ofstream gevp_out_file(combine("plots/%s/gevp_out_eff.xmg",arg[1]).c_str());
  gevp_out_file<<"@type xydy"<<endl;
  gevp_out_file<<effective_mass(gevp.eig_va[0])<<endl;
  gevp_out_file<<"&"<<endl;
  gevp_out_file<<effective_mass(gevp.eig_va[1])<<endl;
  gevp_out_file<<"&"<<endl;  
  
  /*
  
  //////////////////////////////////////
  
  ofstream out(combine("plots/%s/out.xmg",arg[1]).c_str());
  out<<"@type xydy"<<endl;
  out<<Cll<<endl;
  out<<"@s0 legend \"Cll\""<<endl;
  out<<"&"<<endl;
  out<<Css<<endl;
  out<<"@s1 legend \"Css\""<<endl;
  out<<"&"<<endl;
  out<<Dll<<endl;
  out<<"@s2 legend \"Dll\""<<endl;
  out<<"&"<<endl;
  out<<Dls<<endl;
  out<<"@s3 legend \"Dls\""<<endl;
  out<<"&"<<endl;
  out<<Dss<<endl;
  out<<"@s4 legend \"Dss\""<<endl;
  out<<"&"<<endl;
  
  //fit the SL
  jack MSL(njacks),ZL(njacks),ZS(njacks);
  two_pts_SL_fit(MSL,ZL,ZS,ls,ss,tmin_ETA,tmax_ETAP,tmin_ETA,tmax_ETAP,combine("plots/%s/sl_fit.xmg",arg[1]).c_str(),combine("plots/%s/ss_fit.xmg",arg[1]).c_str());

  //fit ss with 2pts
  jack META,Z2ETA;
  TWO_PTS_FIT(META,Z2ETA,ss,tmin_ETA,tmax_ETA,combine("plots/%s/ss_fit.xmg",arg[1]).c_str());
  cout<<META<<endl;
  
  //compute res
  jvec res_ss=ss;
  for(int t=0;t<=TH;t++)
    res_ss[t]-=Z2ETA*exp(-META*TH)*cosh(META*(TH-t))/META;
  
  //fit res with 2pts
  jack METAP,Z2ETAP;
  TWO_PTS_FIT(METAP,Z2ETAP,res_ss,tmin_ETAP,tmax_ETAP,combine("plots/%s/ss_fit_bis.xmg",arg[1]).c_str());
  cout<<METAP<<endl;
  
  TMinuit minu;
  //minu.SetPrintLevel(-1);
  minu.SetFCN(ch2_two_states_migrad_fit);
  minu.DefineParameter(0,"Z2ETA",Z2ETA.med(),Z2ETA.err(),0,0);
  minu.DefineParameter(1,"META",META.med(),META.err(),0,0);
  minu.DefineParameter(2,"Z2ETAP",Z2ETAP.med(),Z2ETAP.err(),0,0);
  minu.DefineParameter(3,"METAP",METAP.med(),METAP.err(),0,0);
  
  c_two_states_fit=new double[TH+1];
  e_two_states_fit=new double[TH+1];
  
  for(int iel=0;iel<=TH;iel++)
    e_two_states_fit[iel]=ss[iel].err();
  
  jack Z2ETA_2st(njacks),META_2st(njacks),Z2ETAP_2st(njacks),METAP_2st(njacks);
  jack ch2(njacks);
  for(int ijack=njacks;ijack>=0;ijack--)
    {
      for(int iel=0;iel<=TH;iel++) c_two_states_fit[iel]=ss[iel][ijack];
      minu.Migrad();
      minu.SetPrintLevel(-1);
      double dum;
      minu.GetParameter(0,Z2ETA_2st.data[ijack],dum);
      minu.GetParameter(1,META_2st.data[ijack],dum);
      minu.GetParameter(2,Z2ETAP_2st.data[ijack],dum);
      minu.GetParameter(3,METAP_2st.data[ijack],dum);
      
      double grad[4],par[4]={Z2ETA_2st[ijack],META_2st[ijack],Z2ETAP_2st[ijack],METAP_2st[ijack]};

      minu.Eval(4,grad,ch2.data[ijack],par,(ijack==njacks)?3:2);
    }
  
  cout<<"bare_masses: "<<smart_print(META)<<endl;
  cout<<"bare_masses: "<<smart_print(METAP)<<endl;
  
  cout<<"ch2: "<<smart_print(ch2)<<endl;
  */

  return 0;
}
