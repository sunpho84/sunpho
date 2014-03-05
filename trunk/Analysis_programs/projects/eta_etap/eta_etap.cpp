#include "include.h"

//#define TWO_PTS_FIT two_pts_migrad_fit
#define TWO_PTS_FIT two_pts_fit

int njacks=72;
int T=96,TH=T/2;

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

int main()
{
  jvec Cll=jvec_load("ll/total",T,njacks,0).simmetrized(1);
  jvec Cll_exa=jvec_load("ll/exact",T,njacks,0).simmetrized(1);
  jvec Dll=jvec_load("disc/total",T,njacks,0).simmetrized(1);
  jvec Dll_exa=jvec_load("disc/exact",T,njacks,0).simmetrized(1);
  jvec Dll_ine=jvec_load("disc/inexact",T,njacks,0).simmetrized(1);
  jvec Dll_res=jvec_load("disc/residue",T,njacks,0).simmetrized(1);
  
  jvec CllWP=jvec_load("llWP/total",T,njacks,0).simmetrized(1);
  jvec CllWP_exa=jvec_load("llWP/exact",T,njacks,0).simmetrized(1);
    
  jvec Dls=jvec_load("disc/total",T,njacks,2).simmetrized(1);
  jvec Dls_ine=jvec_load("disc/inexact",T,njacks,2).simmetrized(1);
  jvec Dls_exa=jvec_load("disc/exact",T,njacks,2).simmetrized(1);
  jvec Dls_res=jvec_load("disc/residue",T,njacks,2).simmetrized(1);
  
  jvec Css=jvec_load("ss/total",T,njacks,0).simmetrized(1);
  jvec Css_exa=jvec_load("ss/exact",T,njacks,0).simmetrized(1);
  jvec Dss=jvec_load("disc/total",T,njacks,4).simmetrized(1);  
  jvec Dss_exa=jvec_load("disc/exact",T,njacks,4).simmetrized(1);  
  
  jvec ll_ine=Cll-2*Dll_ine;
  jvec ll=Cll-2*Dll;
  jvec ls=-sqrt(2)*Dls;
  jvec ls_exa=-sqrt(2)*Dls_exa;
  jvec ss=Css-Dss;
  
  //////////////////////////////////////
  
  ofstream out_Dss("plots/Dss.xmg");
  out_Dss<<"@type xydy"<<endl;
  out_Dss<<Dss_exa<<endl;
  out_Dss<<"&"<<endl;
  out_Dss<<Dss<<endl;
  out_Dss<<"&"<<endl;
  
  ofstream out_Dss_eff("plots/Dss_eff.xmg");
  out_Dss_eff<<"@type xydy"<<endl;
  out_Dss_eff<<effective_mass(Dss)<<endl;
  out_Dss_eff<<"&"<<endl;
  
  ofstream out_Dls("plots/Dls.xmg");
  out_Dls<<"@type xydy"<<endl;
  out_Dls<<Dls<<endl;
  
  ofstream out_Dls_jack("plots/Dls_jack.xmg");
  out_Dls_jack<<"@type xydy"<<endl;
  for(int ijack=0;ijack<njacks;ijack++)
    {
      for(int t=0;t<=48;t++)
	out_Dls_jack<<Dls[t][ijack]<<endl;
      out_Dls_jack<<"&"<<endl;
    }
  
  ofstream out_Dls_exa("plots/Dls_exa.xmg");
  out_Dls_exa<<"@type xydy"<<endl;
  out_Dls_exa<<Dls_exa<<endl;

  ofstream out_Dls_ine("plots/Dls_ine.xmg");
  out_Dls_ine<<"@type xydy"<<endl;
  out_Dls_ine<<Dls_ine<<endl;
  
  ofstream out_Dls_res("plots/Dls_res.xmg");
  out_Dls_res<<"@type xydy"<<endl;
  out_Dls_res<<"&"<<endl;
  out_Dls_res<<Dls_res/Dls<<endl;
  
  //inexact total ll effective mass
  ofstream out_ll_ine("plots/ll_ine.xmg");
  jvec te=ll_ine;
  for(int t=0;t<TH;t+=4) te[t]=ll[t];
  out_ll_ine<<"@type xydy"<<endl;
  out_ll_ine<<aperiodic_effective_mass(ll_ine)<<endl;
  out_ll_ine<<"&"<<endl;
  out_ll_ine<<aperiodic_effective_mass(te)<<endl;

  //relative residue of Disconnected ll
  ofstream out_Dll_res("plots/Dll_res.xmg");
  out_Dll_res<<"@type xydy"<<endl;
  for(int t=0;t<TH;t++) if(Dll_exa[t][njacks]!=0) out_Dll_res<<t<<" "<<Dll_res[t]/Dll[t]<<endl;
  
  //inexact of Disconnected ll
  ofstream out_Dll_ine("plots/Dll_ine.xmg");
  out_Dll_ine<<"@type xydy"<<endl;
  for(int t=0;t<TH;t++) if(Dll_exa[t][njacks]!=0) out_Dll_ine<<t<<" "<<Dll_ine[t]<<endl;
  
  //notes plot for Disconnected ll
  ofstream notes_Dll("plots/Dll.xmg");
  notes_Dll<<"@type xydy"<<endl;
  for(int t=0;t<TH;t++) if(Dll_exa[t][njacks]!=0) notes_Dll<<t<<" "<<Dll_exa[t]<<endl;
  notes_Dll<<"@s0 legend \"Without AMA\""<<endl;
  notes_Dll<<"&"<<endl;
  for(int t=0;t<TH;t++) if(Dll_exa[t][njacks]!=0) notes_Dll<<t+0.6<<" "<<Dll_exa[t]+Dll_res[t]<<endl;
  notes_Dll<<"@s1 legend \"Inexact (same sources)\""<<endl;
  notes_Dll<<"&"<<endl;
  for(int t=0;t<TH;t++) if(Dll_exa[t][njacks]!=0) notes_Dll<<t+1.2<<" "<<Dll_ine[t]<<endl;
  notes_Dll<<"@s2 legend \"Inexact (all sources)\""<<endl;
  notes_Dll<<"&"<<endl;
  for(int t=0;t<TH;t++) if(Dll_exa[t][njacks]!=0) notes_Dll<<t+1.8<<" "<<Dll[t]<<endl;
  notes_Dll<<"@s3 legend \"AMA improved\""<<endl;
  notes_Dll<<"&"<<endl;
  for(int iset=0;iset<4;iset++)
    {
      int col[5]={15,11,6,2,2};
      int type[5]={2,3,4,1,1};
      int size[5]={2,2,2,2,1};
      notes_Dll<<"@s"<<iset<<" line type 0"<<endl;
      notes_Dll<<"@s"<<iset<<" symbol color "<<col[iset]<<endl;
      notes_Dll<<"@s"<<iset<<" errorbar color "<<col[iset]<<endl;
      notes_Dll<<"@s"<<iset<<" symbol "<<type[iset]<<endl;
      notes_Dll<<"@s"<<iset<<" symbol linewidth "<<size[iset]<<endl;
      notes_Dll<<"@s"<<iset<<" errorbar linewidth "<<size[iset]<<endl;
      notes_Dll<<"@s"<<iset<<" errorbar riser linewidth "<<size[iset]<<endl;
    }
  
  //total ll: exact and inexact
  ofstream out_ll("plots/ll.xmg");
  out_ll<<"@type xydy"<<endl;
  out_ll<<ll_ine<<endl;
  out_ll<<"&"<<endl;
  for(int t=0;t<TH;t++) if(Dll_exa[t][njacks]!=0) out_ll<<t<<" "<<ll[t]<<endl;
  out_ll<<"&"<<endl;
  
  //total ss
  ofstream out_ss("plots/ss.xmg");
  out_ss<<"@type xydy"<<endl;
  out_ss<<" "<<ss<<endl;
  out_ss<<"&"<<endl;
  out_ss<<" "<<Css<<endl;
  
  //cancellation between ll and 2*Dll_ine
  ofstream out_ll_cancellation("plots/ll_cancellation.xmg");
  out_ll_cancellation<<"@type xydy"<<endl;
  out_ll_cancellation<<effective_mass(Cll)<<endl;
  out_ll_cancellation<<"&"<<endl;
  out_ll_cancellation<<effective_mass(2*Dll_ine)<<endl;
  out_ll_cancellation<<"&"<<endl;
  
  ofstream ls_aper("plots/ls.xmg");
  ls_aper<<"@type xydy"<<endl;
  ls_aper<<aperiodic_effective_mass(ls)<<endl;
  ls_aper<<"&"<<endl;
  ls_aper<<aperiodic_effective_mass(ls_exa)<<endl;
  
  ofstream ls_and_ss_eff_mass("plots/ls_and_ss_eff_mass.xmg");
  ls_and_ss_eff_mass<<"@type xydy"<<endl;
  ls_and_ss_eff_mass<<effective_mass(ls)<<endl;
  ls_and_ss_eff_mass<<"&"<<endl;
  ls_and_ss_eff_mass<<effective_mass(ss)<<endl;
  ls_and_ss_eff_mass<<"&"<<endl;
  ls_and_ss_eff_mass<<"@type xy"<<endl;
  ls_and_ss_eff_mass<<"0 "<<0.547*a<<endl;
  ls_and_ss_eff_mass<<"48 "<<0.547*a<<endl;
  
  
  ofstream out("plots/out.xmg");
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
  
  //comparison between ll and llWP
  ofstream outWP("plots/Cll_CllWP.xmg");
  outWP<<"@type xydy"<<endl;
  outWP<<effective_mass(Cll)<<endl;
  outWP<<"&"<<endl;
  outWP<<effective_mass(CllWP)<<endl;
  outWP<<"&"<<endl;
  outWP<<effective_mass(CllWP_exa)<<endl;

  //fit the pion
  jack MPION(njacks);
  ofstream out_pion("plots/pion.xmg");
  jvec pion_effmass=effective_mass(Cll);
  for(int t=13;t+5<48;t+=5)
    {
      jack inte_MPION=constant_fit(pion_effmass,t,t+5);
      cout<<t<<" "<<t+5<<" "<<smart_print(inte_MPION)<<endl;
      out_pion<<write_constant_with_error(inte_MPION,t+0.1,t+5-0.1);
      if(t==28) MPION=inte_MPION;
    }
  out_pion<<"@type xydy\n"<<pion_effmass<<endl;
  out_pion<<"&"<<endl;
  out_pion<<"@type xydy\n"<<effective_mass(Cll_exa)<<endl;
  
  ofstream out_ss_conn("plots/ss_conn.xmg");
  out_ss_conn<<"@type xydy"<<endl;
  out_ss_conn<<effective_mass(Css)<<endl;
  out_ss_conn<<"&"<<endl;
  out_ss_conn<<effective_mass(Css_exa)<<endl;

  jack ss_conn_mass=constant_fit(effective_mass(Css),10,47,"plots/ss_conn_fit.xmg");
  cout<<"SS conn mass: "<<smart_print(ss_conn_mass)<<" = "<<smart_print(ss_conn_mass/a)<<" GeV"<<endl;
  
  //fit the SL
  jack MSL(njacks),ZL(njacks),ZS(njacks);
  two_pts_SL_fit(MSL,ZL,ZS,ls,ss,tmin_ETA,tmax_ETAP,tmin_ETA,tmax_ETAP,"plots/sl_fit.xmg","plots/ss_fit.xmg");

  //fit ss with 2pts
  jack META,Z2ETA;
  TWO_PTS_FIT(META,Z2ETA,ss,tmin_ETA,tmax_ETA,"plots/ss_fit.xmg");
  cout<<META<<endl;
  
  //compute res
  jvec res_ss=ss;
  for(int t=0;t<=TH;t++)
    res_ss[t]-=Z2ETA*exp(-META*TH)*cosh(META*(TH-t))/META;
  
  //fit res with 2pts
  jack METAP,Z2ETAP;
  TWO_PTS_FIT(METAP,Z2ETAP,res_ss,tmin_ETAP,tmax_ETAP,"plots/ss_fit_bis.xmg");
  cout<<METAP<<endl;
  
  cout<<"Eta: "<<META/MPION*0.135-0.547<<endl;
  cout<<"Etap: "<<METAP/MPION*0.135-0.958<<endl;
  
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
  
  cout<<"META_2st: "<<smart_print(META*1.75/*MPION*0.135*/)<<"expected 0.547"<<endl;
  cout<<"METAP_2st: "<<smart_print(METAP*1.75/*MPION*0.135*/)<<" expected 0.958"<<endl;
  
  cout<<"META_2st: "<<smart_print(META_2st/MPION*0.135)<<" epected: 0.547"<<endl;
  cout<<"METAP_2st: "<<smart_print(METAP_2st/MPION*0.135)<<" expected -0.958"<<endl;
  cout<<"ch2: "<<smart_print(ch2)<<endl;
  
  return 0;
}
