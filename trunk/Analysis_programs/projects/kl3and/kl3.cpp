#include "include.h"

const int T=96,L=48,TH=L,njacks=83;
const double a=1/1.75,mom_i=0.589278*2*M_PI/L;
const double mu=0.00078,ms=0.0362;
const int atw=1;

const int nonminimal=1;
const int gl_stsep=0;
const int gl_ntsep=6;

jvec total_f_S_move_R1(6,njacks);
jvec total_f_V_move_R1(6,njacks);

int icombo(int tsep,int iv)
{return iv+3*tsep;}

jvec partially_simmetric(jvec in,int tsep)
{
  jvec out(T,njacks);
  for(int t=0;t<=tsep;t++) out[tsep-t]=in[t];
  for(int t=1;t<T-tsep;t++) out[T-t]=in[t+tsep];  

  return out;
}

jvec simmetric(jvec in)
{
  jvec out(T,njacks);
  for(int t=0;t<T;t++) out[(T-t)%T]=in[t];

  return out;
}

jvec load(const char *path,int icombo)
{
  jvec out(T,njacks);
  out.load(combine("data/total_%s",path).c_str(),icombo);
  
  return out;
}

void simple_plot(const char *path,jvec in,int tsep=-1,ios::openmode mode=ios::out)
{
  ofstream out(combine("plots/%s.xmg",path).c_str(),mode);
  out<<"@type xydy"<<endl;
  if(tsep!=-1) out<<"#tsep"<<tsep<<endl;
  out<<in<<"&"<<endl;
}

jvec theta(int t0)
{
  jvec out(T,njacks);
  for(int t=0;t<t0;t++) out[t]=0;
  out[t0]=0.5;
  for(int t=t0+1;t<T;t++) out[t]=1;
  
  return out;
}

jack fit_ratio(const char *path,jvec ratio,int *tint,int tsep,ios::openmode mode=ios::out)
{
  jack f=constant_fit(ratio.subset(0,tsep),tint[0],tint[1],combine("plots/%s_%02d.xmg",path,tsep).c_str());
  ofstream f_out(combine("plots/%s.xmg",path).c_str(),mode);
  f_out<<tsep<<" "<<f<<endl;
  f_out.close();
  
  return f;
}

class global_fit_S_t
{
public:
  int ijack_fit;
  int tmin_two_pts_fit,tmin_three_pts_fit[gl_ntsep],tmax_three_pts_fit[gl_ntsep],tsep[gl_ntsep];
  jack M_Pi,E_Pi,M_K,F;
  jack ZWW_Pi_rest,ZWW_Pi_move,ZWW_K_rest,ZWP_Pi_rest,ZWP_K_rest;
  jvec pion_rest,pion_move,kaon_rest,kp_S0_move_d[gl_ntsep][2];
  jvec point_pion_rest,point_kaon_rest;
};

global_fit_S_t gs;

void ch2_global_S_migrad_fit(int &npar,double *fuf,double &ch,double *p,int flag)
{
  int nd=0;
  ch=0;
  double F=p[0];
  double M_Pi=p[1];
  double E_Pi=p[2];
  double M_K=p[3];
  double ZWW_Pi_rest=p[4];
  double ZWW_Pi_move=p[5];
  double ZWW_K_rest=p[6];
  double ZWP_Pi_rest=p[7];
  double ZWP_K_rest=p[8];
  
  TH_two_pts_fit=TH;
  
  //add the two points pion rest contribution
  //for(int t=gs.tmin_two_pts_fit;t<TH;t++)
  for(int t=10;t<=35;t++)
    {
      double n=gs.pion_rest[t][gs.ijack_fit];
      double p=fun_two_pts_migrad_fit(ZWW_Pi_rest,M_Pi,t);
      double e=gs.pion_rest[t].err();
      
      ch+=sqr((n-p)/e);
      //cout<<" t: "<<t<<" n: "<<n<<" p: "<<p<<" e: "<<e<<endl;
      nd++;
    }
    
  //add the two points pion moving contribution
  //for(int t=gs.tmin_two_pts_fit;t<TH;t++)
  for(int t=6;t<=30;t++)
    {
      double n=gs.pion_move[t][gs.ijack_fit];
      double p=fun_two_pts_migrad_fit(ZWW_Pi_move,E_Pi,t);
      double e=gs.pion_move[t].err();
      
      ch+=sqr((n-p)/e);
      //cout<<" t: "<<t<<" n: "<<n<<" p: "<<p<<" e: "<<e<<endl;
      nd++;
    }
    
  //add the two points kaon at rest contribution
  //for(int t=gs.tmin_two_pts_fit;t<TH;t++)
  for(int t=10;t<=34;t++)
    {
      double n=gs.kaon_rest[t][gs.ijack_fit];
      double p=fun_two_pts_migrad_fit(ZWW_K_rest,M_K,t);
      double e=gs.kaon_rest[t].err();
      
      ch+=sqr((n-p)/e);
      //cout<<" t: "<<t<<" n: "<<n<<" p: "<<p<<" e: "<<e<<endl;
      nd++;
    }
    
  //add the three point function contribution
  for(int isep=gl_stsep;isep<gl_ntsep;isep++)
    for(int id=0;id<2;id++)
      for(int t=gs.tmin_three_pts_fit[isep];t<=gs.tmax_three_pts_fit[isep];t++)
	{
	  double n=gs.kp_S0_move_d[isep][id][t][gs.ijack_fit];
	  double E_Pi_loc=E_Pi;//sqrt(M_Pi*M_Pi+3*mom_i*mom_i); //works much better using the fitted results
	  double p=F*sqrt(ZWW_K_rest*ZWW_Pi_move)*exp(-M_K*t)*exp(-E_Pi_loc*(gs.tsep[isep]-t))/(2*E_Pi_loc*2*M_K);
	  p/=(ms-mu)/(M_K*M_K-M_Pi*M_Pi);
	  double e=gs.kp_S0_move_d[isep][id][t].err();
	  
	  ch+=sqr((n-p)/e);
	  //cout<<" t: "<<t<<" n: "<<n<<" p: "<<p<<" e: "<<e<<endl;
	  nd++;
	}
  
  if(nonminimal)
    {
      //add two points pion at rest with point source
      for(int t=10;t<=40;t++)
	{
	  double n=gs.point_pion_rest[t][gs.ijack_fit];
	  double p=fun_two_pts_migrad_fit(ZWP_Pi_rest,M_Pi,t);
	  double e=gs.point_pion_rest[t].err();
	  
	  ch+=sqr((n-p)/e);
	  //cout<<" t: "<<t<<" n: "<<n<<" p: "<<p<<" e: "<<e<<endl;
	  nd++;
	}
      
      //add two kaon at rest with point source
      for(int t=12;t<=42;t++)
	{
	  double n=gs.point_kaon_rest[t][gs.ijack_fit];
	  double p=fun_two_pts_migrad_fit(ZWP_K_rest,M_K,t);
	  double e=gs.point_kaon_rest[t].err();
	  
	  ch+=sqr((n-p)/e);
	  //cout<<" t: "<<t<<" n: "<<n<<" p: "<<p<<" e: "<<e<<endl;
	  nd++;
	}
    }
}

void global_S_migrad_fit()
{
  TMinuit minu;
  if(gs.ijack_fit!=0) minu.SetPrintLevel(-1);
  
  minu.SetFCN(ch2_global_S_migrad_fit);

  minu.DefineParameter(0,"F",gs.F.med(),gs.F.err(),0,0);
  minu.DefineParameter(1,"M_Pi",gs.M_Pi.med(),gs.M_Pi.err(),0,0);
  minu.DefineParameter(2,"E_Pi",gs.E_Pi.med(),gs.E_Pi.err(),0,0);
  minu.DefineParameter(3,"M_K",gs.M_K.med(),gs.M_K.err(),0,0);
  minu.DefineParameter(4,"ZWW_Pi_rest",gs.ZWW_Pi_rest.med(),gs.ZWW_Pi_rest.err(),0,0);
  minu.DefineParameter(5,"ZWW_Pi_move",gs.ZWW_Pi_move.med(),gs.ZWW_Pi_move.err(),0,0);
  minu.DefineParameter(6,"ZWW_K_rest",gs.ZWW_K_rest.med(),gs.ZWW_K_rest.err(),0,0);
  minu.DefineParameter(7,"ZWP_Pi_rest",gs.ZWP_Pi_rest.med(),gs.ZWP_Pi_rest.err(),0,0);
  minu.DefineParameter(8,"ZWP_K_rest",gs.ZWP_K_rest.med(),gs.ZWP_K_rest.err(),0,0);
  if(nonminimal==0) minu.FixParameter(7);
  if(nonminimal==0) minu.FixParameter(8);
  
  for(gs.ijack_fit=0;gs.ijack_fit<=njacks;gs.ijack_fit++)
    {
      minu.Migrad();  /////////////HET BACK
      double dum;
      minu.GetParameter(0,gs.F[gs.ijack_fit],dum);
      minu.GetParameter(1,gs.M_Pi[gs.ijack_fit],dum);
      minu.GetParameter(2,gs.E_Pi[gs.ijack_fit],dum);
      minu.GetParameter(3,gs.M_K[gs.ijack_fit],dum);
      minu.GetParameter(4,gs.ZWW_Pi_rest[gs.ijack_fit],dum);
      minu.GetParameter(5,gs.ZWW_Pi_move[gs.ijack_fit],dum);
      minu.GetParameter(6,gs.ZWW_K_rest[gs.ijack_fit],dum);
      minu.GetParameter(7,gs.ZWP_Pi_rest[gs.ijack_fit],dum);
      minu.GetParameter(8,gs.ZWP_K_rest[gs.ijack_fit],dum);
    }
}

jack study_single_tsep(int tsep,ios::openmode mode=ios::out)
{
  int it=tsep/4;
  
  //fit pion mass
  jvec pion_rest=load("pion-00WW",0).simmetrized(1),pion_rest_alt=load("pion-00WP",0).simmetrized(1);
  jvec pion_rest_alt_bf=load("pion-00WP",0);
  simple_plot("pion_rest_effmass",effective_mass(pion_rest));
  simple_plot("pion_rest",pion_rest);
  simple_plot("pion_rest_alt",pion_rest_alt_bf);
  jack M_Pi=constant_fit(effective_mass(pion_rest_alt),15,T/2-1,"plots/pion_rest_alt_effmass.xmg");
  cout<<"Most precise pion mass: "<<smart_print(M_Pi)<<endl;
  jack ZWW_Pi_rest(njacks),dum_rest(njacks);
  two_pts_fit(dum_rest,ZWW_Pi_rest,pion_rest,15,T/2-1);
  jack ZWP_Pi_rest(njacks),dum_rest_alt(njacks);
  two_pts_fit(dum_rest_alt,ZWP_Pi_rest,pion_rest_alt,15,T/2-1);
  jack ZW_Pi_rest=sqrt(ZWW_Pi_rest),ZP_Pi_rest=ZWP_Pi_rest/ZW_Pi_rest;
  jvec pion_rest_reco(T,njacks);
  for(int t=0;t<T;t++) pion_rest_reco[t]=ZW_Pi_rest*ZW_Pi_rest*exp(-M_Pi*t)/(2*M_Pi);
  simple_plot("pion_rest_reco",pion_rest_reco);
  
  //fit pion energy
  jvec pion_move=load("pion-01WW",0).simmetrized(1),pion_move_alt=load("pion-01WP",0).simmetrized(1);
  simple_plot("pion_move_effmass",effective_mass(pion_move));
  simple_plot("pion_move",pion_move);
  jack E_Pi_fit=constant_fit(effective_mass(pion_move_alt),15,T/2-1,"plots/pion_move_alt.xmg");
  jack ZWW_Pi_move(njacks),dum_move(njacks);
  two_pts_fit(dum_move,ZWW_Pi_move,pion_move,15,T/2-2);
  jack ZWP_Pi_move(njacks),dum_move_alt(njacks);
  two_pts_fit(dum_move_alt,ZWP_Pi_move,pion_move_alt,15,T/2-2);
  jack ZW_Pi_move=sqrt(ZWW_Pi_move),ZP_Pi_move=ZWP_Pi_move/ZW_Pi_move;
  
  //fit kaon mass
  jvec kaon_rest=load("kaon-00WW",0).simmetrized(1),kaon_rest_alt=load("kaon-00WP",0).simmetrized(1);
  simple_plot("kaon_rest_effmass",effective_mass(kaon_rest));
  simple_plot("kaon_rest",kaon_rest);
  jack M_K=constant_fit(effective_mass(kaon_rest_alt),14,T/2-1,"plots/kaon_rest_alt.xmg");  
  jack ZWW_K_rest(njacks),dum_K_rest(njacks);
  two_pts_fit(dum_K_rest,ZWW_K_rest,kaon_rest,10,T/2-2);
  jack ZWP_K_rest(njacks),dum_K_rest_alt(njacks);
  two_pts_fit(dum_K_rest_alt,ZWP_K_rest,kaon_rest_alt,14,T/2-2);
  jack ZW_K_rest=sqrt(ZWW_K_rest),ZP_K_rest=ZWP_K_rest/ZW_K_rest;
  cout<<"M_K: "<<smart_print(M_K)<<endl;
  
  //compare Z of pions, moving and at rest
  cout<<"ZP_Pi at rest: "<<smart_print(ZP_Pi_rest)<<endl;
  cout<<"ZP_Pi moving:  "<<smart_print(ZP_Pi_move)<<endl;
  cout<<"ZP_K at rest:  "<<smart_print(ZP_K_rest)<<endl;
  cout<<"ZK_Pi at rest:  "<<smart_print(ZW_Pi_rest)<<endl;
  cout<<"ZW_Pi at rest: "<<smart_print(ZW_Pi_rest)<<endl;
  cout<<"ZW_K at rest:  "<<smart_print(ZW_K_rest)<<endl;
  
  //plot the ratios of WW/WP
  ofstream out("/tmp/WWWP.xmg");
  out<<"@type xydy"<<endl;
  out<<pion_rest/pion_rest_alt<<"&"<<endl;
  out<<kaon_rest/kaon_rest_alt<<"&"<<endl;
  
  //verify dispertion relation
  jack E_Pi=sqrt(M_Pi*M_Pi+3*mom_i*mom_i);
  cout<<"E_Pi: "<<smart_print(E_Pi)<<" = "<<smart_print(E_Pi_fit)<<endl;
  
  //compute kinematic factors
  jack P0=M_K+E_Pi,Q0=M_K-E_Pi;
  jack PK(njacks),QK(njacks);
  PK=+mom_i;
  QK=-mom_i;
  jack delta=P0*QK-Q0*PK;
  
  ////////////////////////////////// renormalization //////////////////////////////////
  
  //load insertion of V0 between pions at rest
  jvec pp_V0_rest=(load("zpa-00",icombo(tsep,1))-atw*simmetric(load("zpa-00",icombo(T-tsep,1))))/(atw+1);
  simple_plot(combine("pp_V0_rest_%02d",tsep).c_str(),pp_V0_rest,tsep,mode);
  jvec pp_reno_rest=pp_V0_rest;
  for(int t=0;t<=tsep;t++) pp_reno_rest[t]/=-pion_rest[T/4]*exp(-M_Pi*(tsep-T/4));
  for(int t=tsep+1;t<T;t++) pp_reno_rest[t]/=pion_rest[T/4]*exp(-M_Pi*(T-tsep-T/4));
  simple_plot("pp_reno_rest",pp_reno_rest,tsep,mode);

  //load insertion of V0 between kaons at rest (on s)
  jvec kk_V0_rest_a=(load("zka-00",icombo(tsep,1))-atw*simmetric(load("zka-00",icombo(T-tsep,1))))/(atw+1);
  jvec kk_reno_rest_a=kk_V0_rest_a;
  for(int t=0;t<=tsep;t++) kk_reno_rest_a[t]/=-kaon_rest[T/4]*exp(-M_K*(tsep-T/4));  
  for(int t=tsep+1;t<T;t++) kk_reno_rest_a[t]/=kaon_rest[T/4]*exp(-M_K*(T-tsep-T/4));
  simple_plot("kk_reno_rest_a",kk_reno_rest_a,tsep,mode);

  //load insertion of V0 between kaons at rest (on l)
  jvec kk_V0_rest_b=(load("zkb-00",icombo(tsep,1))-atw*simmetric(load("zkb-00",icombo(T-tsep,1))))/(atw+1);
  jvec kk_reno_rest_b=kk_V0_rest_b;
  for(int t=0;t<=tsep;t++) kk_reno_rest_b[t]/=-kaon_rest[T/4]*exp(-M_K*(tsep-T/4));  
  for(int t=tsep+1;t<T;t++) kk_reno_rest_b[t]/=kaon_rest[T/4]*exp(-M_K*(T-tsep-T/4));
  simple_plot("kk_reno_rest_b",kk_reno_rest_b,tsep,mode);
  
  //load insertion of V0 between pions moving
  jvec pp_V0_move=(load("zpa-11",icombo(tsep,1))-atw*simmetric(load("zpa-11",icombo(T-tsep,1))))/(atw+1);
  jvec pp_reno_move=pp_V0_move;
  for(int t=0;t<=tsep;t++) pp_reno_move[t]/=-pion_move[T/4]*exp(-E_Pi*(tsep-T/4));
  for(int t=tsep+1;t<T;t++) pp_reno_move[t]/=pion_move[T/4]*exp(-E_Pi*(T-tsep-T/4));
  simple_plot("pp_reno_move",pp_reno_move,tsep,mode);
  
  //verify the ratio between pion correlators at rest and in motion
  jvec pp_V0_move_fr_rest=(pp_V0_move/pion_move[tsep])/(pp_V0_rest/pion_rest[tsep]);
  //this should go as the ratio of P0...
  simple_plot("pp_V0_move_fr_rest",pp_V0_move_fr_rest,tsep,mode);
  
  //define tilded
  jvec pion_rest_tilde(T,njacks);
  jvec pion_move_tilde(T,njacks);
  jvec kaon_rest_tilde(T,njacks);
  for(int t=0;t<=T/2;t++)
    {
      pion_rest_tilde[t]=pion_rest[t]-pion_rest[T/2]*exp(-M_Pi*(T/2-t))/2;
      pion_move_tilde[t]=pion_move[t]-pion_move[T/2]*exp(-E_Pi*(T/2-t))/2;
      kaon_rest_tilde[t]=kaon_rest[t]-kaon_rest[T/2]*exp(-M_K*(T/2-t))/2;
    }
  for(int t=T/2;t<T;t++)
    {
      pion_rest_tilde[t]=pion_rest[T/2]*exp(-M_Pi*(t-T/2))/2;
      pion_move_tilde[t]=pion_move[T/2]*exp(-E_Pi*(t-T/2))/2;
      kaon_rest_tilde[t]=kaon_rest[T/2]*exp(-M_K*(t-T/2))/2;
    }
  
  simple_plot("pion_rest_tilde",pion_rest_tilde);

  //normalization
  jvec norm_a=pp_V0_rest*kk_V0_rest_a;
  simple_plot("norm_a",norm_a,tsep,mode);
  jvec norm_b=pp_V0_rest*kk_V0_rest_b;
  simple_plot("norm_b",norm_b,tsep,mode);
  
  //Zv
  jvec Zv_pion_rest_corr=-pion_rest_tilde[tsep]/pp_V0_rest;
  simple_plot("Zv_pion_rest",-1/Zv_pion_rest_corr.subset(0,tsep),tsep,mode);
  jvec Zv_pion_move_corr=-pion_move_tilde[tsep]/pp_V0_move;
  simple_plot("Zv_pion_move",Zv_pion_move_corr.subset(0,tsep),tsep,mode);
  jvec Zv_kaon_rest_a_corr=-kaon_rest_tilde[tsep]/kk_V0_rest_a;
  simple_plot("Zv_kaon_rest_a",Zv_kaon_rest_a_corr.subset(0,tsep),tsep,mode);
  jvec Zv_kaon_rest_b_corr=-kaon_rest_tilde[tsep]/kk_V0_rest_b;
  simple_plot("Zv_kaon_rest_b",Zv_kaon_rest_b_corr.subset(0,tsep),tsep,mode);
  
  //extract Zv from pion and kaon_a at rest
  jack Zv_pion=constant_fit(Zv_pion_rest_corr,tsep/2-tsep/4,tsep/2+tsep/4);
  jack Zv_pion_move=constant_fit(Zv_pion_move_corr,tsep/2-tsep/4,tsep/2+tsep/4);
  jack Zv_kaon=constant_fit(Zv_kaon_rest_a_corr,tsep/2-tsep/4,tsep/2+tsep/4);
  jack Zv_kaon_alt=constant_fit(Zv_kaon_rest_b_corr,tsep/2-tsep/4,tsep/2+tsep/4);
  jack Zv_mixed=constant_fit(sqrt(Zv_pion_rest_corr*Zv_kaon_rest_a_corr),tsep/2-tsep/4,tsep/2+tsep/4);
  cerr<<"ZV: "<<tsep<<" pi_rest: "<<smart_print(Zv_pion)<<" pi_move(unused): "<<smart_print(Zv_pion_move)<<
    " k_rest(s): "<<smart_print(Zv_kaon)<<" k_rest(l)(unused): "<<smart_print(Zv_kaon_alt)<<endl;
  cerr<<"ZV: "<<tsep<<" "<<smart_print(Zv_pion)<<" "<<smart_print(Zv_kaon)<<endl;
  cerr<<"ZV_mixed: "<<tsep<<" "<<smart_print(Zv_mixed)<<" "<<smart_print(sqrt(Zv_kaon*Zv_pion))<<
    "in the middle: "<<sqrt(Zv_pion_rest_corr*Zv_kaon_rest_a_corr)[tsep/2]<<endl;
  cerr<<"diff: "<<smart_print((Zv_pion_rest_corr-Zv_kaon_rest_a_corr)[tsep/2])<<endl;
  //////////////////////////////////// kl3 //////////////////////////////////
  
  //load kl3 with insertions
  jvec kp_VK_move=(load("kl3-01",icombo(tsep,0))+atw*simmetric(load("kl3-01",icombo(T-tsep,0))))/(atw+1);
  jvec kp_V0_move=(load("kl3-01",icombo(tsep,1))-atw*simmetric(load("kl3-01",icombo(T-tsep,1))))/(atw+1);
  jvec kp_S0_move=(load("kl3-01",icombo(tsep,2))+atw*simmetric(load("kl3-01",icombo(T-tsep,2))))/(atw+1);
  jvec kp_VK_rest=(load("kl3-00",icombo(tsep,0))+atw*simmetric(load("kl3-00",icombo(T-tsep,0))))/(atw+1);
  jvec kp_V0_rest=(load("kl3-00",icombo(tsep,1))-atw*simmetric(load("kl3-00",icombo(T-tsep,1))))/(atw+1);

  //20,24,28,32,36,40
  static bool re_fl=false;
  if(re_fl==false)
    for(int isep=gl_stsep;isep<gl_ntsep;isep++)
      {
	gs.kp_S0_move_d[isep][0]=load("kl3-01",icombo((isep+5)*4,2));
	gs.kp_S0_move_d[isep][1]=simmetric(load("kl3-01",icombo(T-(isep+5)*4,2)));
      }
  re_fl=true;
  
  //check that VK is 0
  simple_plot("kp_VK_move",kp_VK_move,tsep,mode);
  simple_plot("kp_VK_rest",kp_VK_rest,tsep,mode);

  //compare moving and rest kl3
  simple_plot("kp_V0_move",kp_V0_move,tsep,mode);
  simple_plot("kp_V0_rest",kp_V0_rest,tsep,mode);
  
  //solve the correlator for V and for S
  jvec kp_V_rest=kp_V0_move/(M_Pi+M_K);
  jvec kp_V_move=(kp_V0_move*QK-Q0*kp_VK_move)/delta;
  jvec kp_F_move=(kp_V0_move*PK-P0*kp_VK_move)/delta;
  jvec kp_S_move=kp_S0_move*(ms-mu)/(M_K*M_K-M_Pi*M_Pi);
  
  ////////////////////// * ANDREAS comparison * ///////////////////////////
  
  //build ratios
  jvec kp_V0_rest_R1=sqrt(4*M_Pi*M_K*kp_V0_rest*partially_simmetric(kp_V0_rest,tsep)/
				  (pion_rest_tilde[tsep]*kaon_rest_tilde[tsep]));
  jvec kp_V0_move_R1=sqrt(4*E_Pi*M_K*kp_V0_move*partially_simmetric(kp_V0_move,tsep)/
				   (pion_move_tilde[tsep]*kaon_rest_tilde[tsep]));
  jvec kp_VK_move_R1=sqrt(4*E_Pi*M_K*kp_VK_move*partially_simmetric(kp_VK_move,tsep)/
				 (pion_move_tilde[tsep]*kaon_rest_tilde[tsep]));
  jvec kp_V0_rest_R2=2*sqrt(M_Pi*M_K)*sqrt(kp_V0_rest*partially_simmetric(kp_V0_rest,tsep)/(pp_V0_rest*kk_V0_rest_a));
  jvec kp_V0_move_R2=2*sqrt(E_Pi*M_K)*sqrt(kp_V0_move*partially_simmetric(kp_V0_move,tsep)/(pp_V0_move*kk_V0_rest_a));
  jvec kp_VK_move_R2=2*sqrt(E_Pi*M_K)*sqrt(kp_VK_move*partially_simmetric(kp_VK_move,tsep)/(pp_V0_move*kk_V0_rest_a));
  
  jvec kp_VK_move_R4_BOB=kp_VK_move;
  for(int t=0;t<tsep;t++)
    kp_VK_move_R4_BOB[t]/=ZW_Pi_rest*ZW_K_rest*exp(-M_K*t)*exp(-E_Pi*(tsep-t))/(2*E_Pi*2*M_K);
  
  simple_plot("kp_V0_rest_R1",kp_V0_rest_R1.subset(0,tsep),tsep,mode);
  simple_plot("kp_V0_move_R1",kp_V0_move_R1.subset(0,tsep),tsep,mode);
  simple_plot("kp_VK_move_R1",kp_VK_move_R1.subset(0,tsep),tsep,mode);
  simple_plot("kp_V0_rest_R2",kp_V0_rest_R2.subset(0,tsep),tsep,mode);
  simple_plot("kp_V0_move_R2",kp_V0_move_R2.subset(0,tsep),tsep,mode);
  simple_plot("kp_VK_move_R2",kp_VK_move_R2.subset(0,tsep),tsep,mode);

  simple_plot(combine("kp_VK_move_R4_BOB_tsep_%02d",tsep).c_str(),kp_VK_move_R4_BOB.subset(0,tsep),tsep,mode);
  
  //fit V move ratios             0     4     8     12    16     20     24     28     32      36     40
  int f_V0_move_R1_tint[11][2]={{0,0},{1,3},{4,6},{4,8},{5,11},{6,14},{6,20},{4,24},{14,17},{8,28},{10,35}};
  int f_VK_move_R1_tint[11][2]={{0,0},{1,3},{4,6},{4,8},{5,11},{6,14},{6,20},{4,24},{14,17},{8,28},{10,35}};
  jack f_V0_move_R1=fit_ratio("f_V0_move_R1",kp_V0_move_R1,f_V0_move_R1_tint[it],tsep,mode);
  jack f_VK_move_R1=fit_ratio("f_VK_move_R1",kp_VK_move_R1,f_VK_move_R1_tint[it],tsep,mode);
  int f_V0_move_R2_tint[11][2]={{0,0},{1,3},{4,6},{4,8},{5,11},{6,14},{6,20},{4,24},{14,17},{8,28},{10,35}};
  int f_VK_move_R2_tint[11][2]={{0,0},{1,3},{4,6},{4,8},{5,11},{6,14},{6,20},{4,24},{14,17},{8,28},{10,35}};
  jack f_V0_move_R2=fit_ratio("f_V0_move_R2",kp_V0_move_R2,f_V0_move_R2_tint[it],tsep,mode);
  jack f_VK_move_R2=fit_ratio("f_VK_move_R2",kp_VK_move_R2,f_VK_move_R2_tint[it],tsep,mode);
  //cout<<"V0 tsep "<<tsep<<" at timeslice 4: "<<kp_V0_move_R1[4]<<endl;
  //cout<<"VK tsep "<<tsep<<" at timeslice 4: "<<kp_VK_move_R1[4]<<endl;

  /////////////////////////////////////////////////////////////////////////

  //build ratios
  jvec kp_V_rest_R1=Zv_mixed*sqrt(4*M_Pi*M_K*kp_V_rest*partially_simmetric(kp_V_rest,tsep)/
				 (pion_rest_tilde[tsep]*kaon_rest_tilde[tsep]));
  jvec kp_V_move_R1=Zv_mixed*sqrt(4*E_Pi*M_K*kp_V_move*partially_simmetric(kp_V_move,tsep)/
				 (pion_move_tilde[tsep]*kaon_rest_tilde[tsep]));
  jvec kp_F_move_R1=-Zv_mixed*sqrt(4*E_Pi*M_K*kp_F_move*partially_simmetric(kp_F_move,tsep)/
				 (pion_move_tilde[tsep]*kaon_rest_tilde[tsep]));
  
  //jvec kp_S_move_R1=kp_S0_move*partially_simmetric(kp_S0_move,tsep)/(pion_move_tilde[tsep]*kaon_rest_tilde[tsep]);
  jvec kp_S_move_R1=sqrt(4*E_Pi*M_K*kp_S_move*partially_simmetric(kp_S_move,tsep)/
			 (pion_move_tilde[tsep]*kaon_rest_tilde[tsep]));
  
  jvec kp_V_rest_R2=sqrt(4*M_Pi*M_K*kp_V_rest*partially_simmetric(kp_V_rest,tsep)/(pp_V0_rest*kk_V0_rest_a));
  jvec kp_V_move_R2=sqrt(4*E_Pi*M_K*kp_V_move*partially_simmetric(kp_V_move,tsep)/(pp_V0_move*kk_V0_rest_a));
  jvec kp_V_rest_R3=-4*sqrt(M_Pi*M_K)*kp_V_rest/pion_rest_tilde[tsep];
  jvec kp_V_move_R3=-4*sqrt(E_Pi*M_K)*kp_V_move/pion_move_tilde[tsep];
  for(int t=0;t<tsep;t++)
    {
      kp_V_rest_R3[t]*=sqrt(kaon_rest_tilde[tsep-t]*pion_rest_tilde[t]*pion_rest_tilde[tsep]);
      kp_V_rest_R3[t]/=sqrt(pion_rest_tilde[tsep-t]*kaon_rest_tilde[t]*kaon_rest_tilde[tsep]);
      kp_V_move_R3[t]*=sqrt(kaon_rest_tilde[tsep-t]*pion_move_tilde[t]*pion_move_tilde[tsep]);
      kp_V_move_R3[t]/=sqrt(pion_move_tilde[tsep-t]*kaon_rest_tilde[t]*kaon_rest_tilde[tsep]);
    }
  jvec kp_V_rest_R4=-kp_V_rest*Zv_mixed;
  jvec kp_V_move_R4=-kp_V_move*Zv_mixed;
  jvec kp_S_move_R4=kp_S_move;
  for(int t=0;t<tsep;t++)
    {
      kp_V_rest_R4[t]/=ZW_Pi_rest*ZW_K_rest*exp(-M_K*t)*exp(-M_Pi*(tsep-t))/(2*M_Pi*2*M_K);
      kp_V_move_R4[t]/=ZW_Pi_move*ZW_K_rest*exp(-M_K*t)*exp(-E_Pi*(tsep-t))/(2*E_Pi*2*M_K);
      kp_S_move_R4[t]/=ZW_Pi_move*ZW_K_rest*exp(-M_K*t)*exp(-E_Pi*(tsep-t))/(2*E_Pi*2*M_K);
    }
  simple_plot("kp_V_rest_R1",kp_V_rest_R1.subset(0,tsep),tsep,mode);
  simple_plot("kp_V_rest_R2",kp_V_rest_R2.subset(0,tsep),tsep,mode);
  simple_plot("kp_V_rest_R3",kp_V_rest_R3.subset(0,tsep),tsep,mode);
  simple_plot("kp_V_rest_R4",kp_V_rest_R4.subset(0,tsep),tsep,mode);
  
  /////////////////////////////////////////// plots /////////////////////////////////////  
  
  //fit V rest ratios            0     4     8     12    16     20     24     28      32      36      40
  int f_V_rest_R1_tint[11][2]={{0,0},{1,3},{4,6},{4,8},{5,11},{6,14},{6,20},{4,24},{14,17},{8,28},{10,35}};
  int f_V_rest_R2_tint[11][2]={{0,0},{1,3},{4,6},{4,8},{5,11},{5,14},{6,20},{4,24},{14,17},{8,28},{10,35}};
  int f_V_rest_R3_tint[11][2]={{0,0},{1,3},{4,6},{4,8},{5,11},{5,14},{6,20},{4,24},{14,17},{8,28},{10,35}};
  int f_V_rest_R4_tint[11][2]={{0,0},{1,3},{4,6},{4,8},{5,11},{5,14},{6,20},{4,24},{14,17},{8,28},{10,35}};
  jack f_V_rest_R1=fit_ratio("f_V_rest_R1",kp_V_rest_R1,f_V_rest_R1_tint[it],tsep,mode);
  jack f_V_rest_R2=fit_ratio("f_V_rest_R2",kp_V_rest_R2,f_V_rest_R2_tint[it],tsep,mode);
  jack f_V_rest_R3=fit_ratio("f_V_rest_R3",kp_V_rest_R3,f_V_rest_R3_tint[it],tsep,mode);
  jack f_V_rest_R4=fit_ratio("f_V_rest_R4",kp_V_rest_R4,f_V_rest_R4_tint[it],tsep,mode);

  //fit V move ratios            0     4     8     12    16    20     24     28      32      36      40
  int f_V_move_R1_tint[11][2]={{0,0},{1,3},{4,6},{7,8},{7,9},{9,11},{9,15},{8,22},{17,22},{12,31},{13,33}};
  int f_V_move_R2_tint[11][2]={{0,0},{1,3},{4,6},{7,8},{8,10},{9,13},{9,18},{8,22},{17,22},{12,31},{13,33}};
  int f_V_move_R3_tint[11][2]={{0,0},{1,3},{4,6},{7,8},{8,10},{9,13},{9,18},{8,22},{17,22},{12,31},{13,33}};
  int f_V_move_R4_tint[11][2]={{0,0},{1,3},{4,6},{7,8},{8,10},{9,13},{9,18},{8,22},{17,22},{12,31},{13,33}};
  jack f_V_move_R1=fit_ratio("f_V_move_R1",kp_V_move_R1,f_V_move_R1_tint[it],tsep,mode);
  jack f_V_move_R2=fit_ratio("f_V_move_R2",kp_V_move_R2,f_V_move_R2_tint[it],tsep,mode);
  jack f_V_move_R3=fit_ratio("f_V_move_R3",kp_V_move_R3,f_V_move_R3_tint[it],tsep,mode);
  jack f_V_move_R4=fit_ratio("f_V_move_R4",kp_V_move_R4,f_V_move_R4_tint[it],tsep,mode);
  
  //fit V move ratios            0     4     8     12    16    20     24     28      32      36      40
  int f_F_move_R1_tint[11][2]={{0,0},{1,3},{4,6},{7,8},{7,9},{9,11},{9,15},{8,22},{17,22},{12,31},{13,33}};
  jack f_F_move_R1=fit_ratio("f_F_move_R1",kp_F_move_R1,f_F_move_R1_tint[it],tsep,mode);
  
  //fit S move ratios            0     4     8     12    16    20     24     28     32       36      40
  int f_S_move_R1_tint[11][2]={{0,0},{1,3},{4,6},{7,8},{8,8},{8,12},{7,17},{9,19},{11,21},{12,24},{13,31}};
  //int f_S_move_R2_tint[11][2]={{0,0},{1,3},{4,6},{7,8},{8,10},{9,13},{8,19},{8,23},{10,23},{12,24},{13,33}};
  //int f_S_move_R3_tint[11][2]={{0,0},{1,3},{4,6},{7,8},{8,10},{9,13},{8,19},{8,23},{10,23},{12,24},{13,33}};
  int f_S_move_R4_tint[11][2]={{0,0},{1,3},{4,6},{7,8},{8,10},{9,13},{8,19},{8,23},{10,23},{12,24},{13,33}};
  jack f_S_move_R1=fit_ratio("f_S_move_R1",kp_S_move_R1,f_S_move_R1_tint[it],tsep,mode);
  //jack f_S_move_R2=fit_ratio("f_S_move_R2",kp_S_move_R2,f_S_move_R2_tint[it],tsep,mode);
  //jack f_S_move_R3=fit_ratio("f_S_move_R3",kp_S_move_R3,f_S_move_R3_tint[it],tsep,mode);
  jack f_S_move_R4=fit_ratio("f_S_move_R4",kp_S_move_R4,f_S_move_R4_tint[it],tsep,mode);
  
  jvec kp_SV_diff_move_R1=kp_S_move_R1-kp_V_move_R1;
  simple_plot("kp_SV_diff_move_R1",kp_SV_diff_move_R1.subset(0,tsep),tsep,mode);
  jvec kp_SV_ratio_move_R1=kp_S_move_R1/kp_V_move_R1-1;
  simple_plot("kp_SV_ratio_move_R1",kp_SV_ratio_move_R1.subset(0,tsep),tsep,mode);
  
  jack combo_R1=(f_S_move_R1/sqr(f_S_move_R1.err())+f_V_move_R1/sqr(f_V_move_R1.err()))/
    (1/sqr(f_S_move_R1.err())+1/sqr(f_V_move_R1.err()));
  //jack combo_R2=(f_S_move_R2/sqr(f_S_move_R2.err())+f_V_move_R2/sqr(f_V_move_R2.err()))/
  //(1/sqr(f_S_move_R2.err())+1/sqr(f_V_move_R2.err()));
  
  cout<<"S only: "<<smart_print(f_S_move_R1)<<" Vector only: "<<smart_print(f_V_move_R1)<<endl;
  
  cout<<"combo_R1: "<<smart_print(combo_R1)<<endl;
  //cout<<"combo_R2: "<<smart_print(combo_R2)<<endl;
  
  //////////////////////////// DAVID comparison //////////////////////////
  
  gs.F=f_S_move_R1;
  gs.M_Pi=M_Pi;
  gs.E_Pi=E_Pi;
  gs.M_K=M_K;
  gs.ZWW_Pi_rest=ZWW_Pi_rest;
  gs.ZWW_Pi_move=ZWW_Pi_move;
  gs.ZWP_Pi_rest=ZWP_Pi_rest;
  gs.ZWW_K_rest=ZWW_K_rest;
  gs.ZWP_K_rest=ZWP_K_rest;
  for(int isep=gl_stsep;isep<gl_ntsep;isep++) gs.tsep[isep]=(isep+5)*4;
  gs.tmin_two_pts_fit=15;  
  
  //    tsep               0     4     8     12    16     20     24     28     32     36     40 
  int DAVID_tint[11][2]={{0,0},{1,3},{4,6},{4,8},{5,11},{7,15},{6,20},{8,22},{6,27},{6,30},{9,33}};

  for(int isep=gl_stsep;isep<gl_ntsep;isep++)
    {
      gs.tmin_three_pts_fit[isep]=DAVID_tint[isep+5][0];
      gs.tmax_three_pts_fit[isep]=DAVID_tint[isep+5][1];
    }

  gs.pion_rest=pion_rest;
  gs.point_pion_rest=pion_rest_alt;
  gs.pion_move=pion_move;
  gs.kaon_rest=kaon_rest;
  gs.point_kaon_rest=kaon_rest_alt;

  global_S_migrad_fit();

  jack ch2_total_OUR(njacks),ch2_total_DAVID(njacks);
  ch2_total_OUR=ch2_total_DAVID=0;
  int ndof=0;
  
  cout<<"ZpWW00: "<<(sqrt(gs.ZWW_Pi_rest))<<endl;
  cout<<"ZpWP00: "<<(sqrt(gs.ZWP_Pi_rest))<<endl;
  cout<<"ZkWP00: "<<(sqrt(gs.ZWP_K_rest))<<endl;
  cout<<"ZpWW01: "<<(sqrt(gs.ZWW_Pi_move))<<endl;
  cout<<"ZkWW00: "<<(sqrt(gs.ZWW_K_rest))<<endl;
  cout<<"m_K: "<<(gs.M_K)<<endl;
  cout<<"m_pi: "<<(gs.M_Pi)<<endl;
  cout<<"Ep: "<<(gs.E_Pi)<<endl;
  cout<<"f_p: "<<(gs.F)<<endl;
  cout<<"---"<<endl;
  
  jack ZpWW00=fill_gauss(28554.1,56.3,0,njacks);
  jack ZpWW01=fill_gauss(26346.2,70.3,0,njacks);
  jack ZkWW00=fill_gauss(30733.9,50.1,0,njacks);
  jack m_K=fill_gauss(0.288650,0.000171,0,njacks);
  jack m_pi=fill_gauss(0.0804754,0.0001113,0,njacks);
  jack Ep=fill_gauss(0.156852,0.000327,0,njacks);
  jack f_p=fill_gauss(0.967734,0.003127,0,njacks);

  //add the three points contribution
  for(int isep=gl_stsep;isep<gl_ntsep;isep++)
    for(int id=0;id<2;id++)
      {
	cout<<"---three points time ordering "<<id<<" time separation "<<isep<<" contribution---"<<endl;
	for(int t=DAVID_tint[isep+5][0];t<=DAVID_tint[isep+5][1];t++)
	{
	  jack three_pts_reco_DAVID=(m_K*m_K-m_pi*m_pi)/(ms-mu)*f_p*ZpWW01*ZkWW00*exp(-m_K*t)*exp(-Ep*(tsep-t))/(2*Ep*2*m_K);
	  jack three_pts_reco_OUR=(gs.F*sqrt(gs.ZWW_K_rest*gs.ZWW_Pi_move)*exp(-gs.M_K*t)*exp(-gs.E_Pi*(gs.tsep[isep]-t))/(2*gs.E_Pi*2*gs.M_K)/((ms-mu)/(gs.M_K*gs.M_K-gs.M_Pi*gs.M_Pi)));
	  
	  jack diff_DAVID=three_pts_reco_DAVID-gs.kp_S0_move_d[isep][id][t].med();
	  jack diff_OUR=three_pts_reco_OUR-gs.kp_S0_move_d[isep][id][t];
	  jack cchi_DAVID=sqr(diff_DAVID/gs.kp_S0_move_d[isep][id][t].err());
	  jack cchi_OUR=sqr(diff_OUR/gs.kp_S0_move_d[isep][id][t].err());
	  
	  ch2_total_OUR+=cchi_OUR;
	  ch2_total_DAVID+=cchi_DAVID;
	  
	  cout<<"t: "<<t<<" "<<(cchi_OUR.med())<<"=[("<<three_pts_reco_OUR.med()<<"-"<<gs.kp_S0_move_d[isep][id][t].med()<<")/"<<gs.kp_S0_move_d[isep][id][t].err()<<"]^2"<<endl;
	  ndof++;
	}
    }
  
  //add the two points pion at rest contribution
  cout<<"---two points pion at rest contribution---"<<endl;
  for(int t=10;t<=35;t++)
    {
      jack n=gs.pion_rest[t];
      jack p_OUR(njacks),p_DAVID(njacks);
      for(int ijack=0;ijack<=njacks;ijack++)
	{
	  p_OUR[ijack]=fun_two_pts_migrad_fit(gs.ZWW_Pi_rest[ijack],gs.M_Pi[ijack],t);
	  p_DAVID[ijack]=fun_two_pts_migrad_fit(ZpWW00[ijack]*ZpWW00[ijack],m_pi[ijack],t);
	}
      double e=gs.pion_rest[t].err();
      jack cchi_OUR=sqr((n-p_OUR)/e);
      ch2_total_OUR+=cchi_OUR;
      ch2_total_DAVID+=sqr((n-p_DAVID)/e);
      
      cout<<"t: "<<t<<" "<<cchi_OUR.med()<<"=[("<<p_OUR.med()<<"-"<<n.med()<<")/"<<e<<"]^2"<<endl;
      ndof++;
    }
    
  //add the two points pion moving contribution
  cout<<"---two points pion moving contribution---"<<endl;
  for(int t=6;t<=30;t++)
    {
      jack n=gs.pion_move[t];
      jack p_OUR(njacks),p_DAVID(njacks);
      for(int ijack=0;ijack<=njacks;ijack++)
	{
	  p_OUR[ijack]=fun_two_pts_migrad_fit(gs.ZWW_Pi_move[ijack],gs.E_Pi[ijack],t);
	  p_DAVID[ijack]=fun_two_pts_migrad_fit(ZpWW01[ijack]*ZpWW01[ijack],Ep[ijack],t);
	}
      double e=gs.pion_move[t].err();
      
      jack cchi_OUR=sqr((n-p_OUR)/e);
      ch2_total_OUR+=cchi_OUR;
      ch2_total_DAVID+=sqr((n-p_DAVID)/e);

      cout<<"t: "<<t<<" "<<cchi_OUR.med()<<"=[("<<p_OUR.med()<<"-"<<n.med()<<")/"<<e<<"]^2"<<endl;
      ndof++;
    }
    
  //add the two points kaon at rest contribution
  cout<<"---two points kaon at rest contribution---"<<endl;
  jvec test=load("kaon-00WW",0);
  for(int t=10;t<=34;t++)
    {
      jack n=gs.kaon_rest[t];
      jack p_OUR(njacks),p_DAVID(njacks);
      for(int ijack=0;ijack<=njacks;ijack++)
	{
	  p_OUR[ijack]=fun_two_pts_migrad_fit(gs.ZWW_K_rest[ijack],gs.M_K[ijack],t);
	  p_DAVID[ijack]=fun_two_pts_migrad_fit(ZkWW00[ijack]*ZkWW00[ijack],m_K[ijack],t);
	}
      double e=gs.kaon_rest[t].err();
      
      //jack cchi_OUR=sqr((n-p_OUR)/e);
      jack cchi_OUR=sqr((test[t]-p_OUR)/test[t].err());
      ch2_total_OUR+=cchi_OUR;
      ch2_total_DAVID+=sqr((n-p_DAVID)/e);
      
      cout<<"t: "<<t<<" "<<cchi_OUR.med()<<"=[("<<p_OUR.med()<<"-"<<n.med()<<")/"<<e<<"]^2"<<endl;
      cout<<"t: "<<t<<" "<<cchi_OUR.med()<<"=[("<<p_OUR.med()<<"-"<<test[t].med()<<")/"<<test[t].err()<<"]^2"<<endl;
      ndof++;
    }

  if(nonminimal)
    {
      //add the two points point pion at rest contribution
      cout<<"---two points point pion at rest contribution---"<<endl;
      for(int t=10;t<=40;t++)
	{
	  jack n=gs.point_pion_rest[t];
	  jack p_OUR(njacks),p_DAVID(njacks);
	  for(int ijack=0;ijack<=njacks;ijack++)
	    {
	      p_OUR[ijack]=fun_two_pts_migrad_fit(gs.ZWP_Pi_rest[ijack],gs.M_Pi[ijack],t);
	      p_DAVID[ijack]=fun_two_pts_migrad_fit(ZpWW00[ijack]*ZpWW00[ijack],m_pi[ijack],t);
	    }
	  double e=gs.point_pion_rest[t].err();
	  jack cchi_OUR=sqr((n-p_OUR)/e);
	  ch2_total_OUR+=cchi_OUR;
	  //ch2_total_DAVID+=sqr((n-p_DAVID)/e);
	  
	  cout<<"t: "<<t<<" "<<cchi_OUR.med()<<"=[("<<p_OUR.med()<<"-"<<n.med()<<")/"<<e<<"]^2"<<endl;
	  ndof++;
	}

      //add the two points point kaon at rest contribution
      cout<<"---two points point kaon at rest contribution---"<<endl;
      for(int t=12;t<=42;t++)
	{
	  jack n=gs.point_kaon_rest[t];
	  jack p_OUR(njacks),p_DAVID(njacks);
	  for(int ijack=0;ijack<=njacks;ijack++)
	    {
	      p_OUR[ijack]=fun_two_pts_migrad_fit(gs.ZWP_K_rest[ijack],gs.M_K[ijack],t);
	      p_DAVID[ijack]=fun_two_pts_migrad_fit(ZpWW00[ijack]*ZpWW00[ijack],m_pi[ijack],t);
	    }
	  double e=gs.point_kaon_rest[t].err();
	  jack cchi_OUR=sqr((n-p_OUR)/e);
	  ch2_total_OUR+=cchi_OUR;
	  //ch2_total_DAVID+=sqr((n-p_DAVID)/e);
	  
	  cout<<"t: "<<t<<" "<<cchi_OUR.med()<<"=[("<<p_OUR.med()<<"-"<<n.med()<<")/"<<e<<"]^2"<<endl;
	  ndof++;
	}
    }
    
  cout<<"Npoints: "<<ndof<<endl;
  cout<<"Total ch2_OUR("<<tsep<<"): "<<smart_print(ch2_total_OUR)<<"/"<<ndof-7<<"="<<smart_print(ch2_total_OUR/(ndof-7))<<endl;
  cout<<"Total ch2_DAVID: "<<ch2_total_DAVID/(ndof-7)<<endl;  
  
  cout<<"F: "<<smart_print(gs.F)<<endl;
  cout<<"M_Pi: "<<smart_print(gs.M_Pi)<<" "<<smart_print(M_Pi)<<endl;
  cout<<"E_Pi: "<<smart_print(gs.E_Pi)<<" "<<smart_print(E_Pi)<<endl;
  cout<<"M_K: "<<smart_print(gs.M_K)<<" "<<smart_print(M_K)<<endl;
  cout<<"ZWW_Pi_rest: "<<smart_print(gs.ZWW_Pi_rest)<<" "<<smart_print(ZWW_Pi_rest)<<endl;
  cout<<"ZWW_Pi_move: "<<smart_print(gs.ZWW_Pi_move)<<" "<<smart_print(ZWW_Pi_move)<<endl;
  cout<<"ZWW_K_rest: "<<smart_print(gs.ZWW_K_rest)<<" "<<smart_print(ZWW_K_rest)<<endl;
  cout<<"ZWP_Pi_rest: "<<smart_print(gs.ZWP_Pi_rest)<<" "<<smart_print(ZWP_Pi_rest)<<endl;
  cout<<"ZWP_K_rest: "<<smart_print(gs.ZWP_K_rest)<<" "<<smart_print(ZWP_K_rest)<<endl;
  
  if(tsep>16) total_f_S_move_R1[tsep/4-5]=f_S_move_R1;
  if(tsep>16) total_f_V_move_R1[tsep/4-5]=f_V_move_R1;
  
  return (tsep==16)?f_S_move_R1:combo_R1;
}

int main()
{
  std::vector<jack> results;
  results.push_back(study_single_tsep(4));
  for(int tsep=8;tsep<=40;tsep+=4) results.push_back(study_single_tsep(tsep,ios::out|ios::app));
  
  jack tot(njacks);
  tot=0;
  double tw=0;
  for(int i=4;i<results.size();i++)
    {
      cout<<i*4+4<<" "<<smart_print(results[i])<<endl;
      double w=1/sqr(results[i].err());
      tot+=results[i]*w;
      tw+=w;
    }
  tot/=tw;
  
  cout<<smart_print(tot)<<" "<<1/sqrt(tw)<<endl;
  
  cout<<"---av S"<<endl;
  cout<<total_f_S_move_R1<<endl;
  cout<<constant_fit(total_f_S_move_R1,0,6)<<endl;
  
  cout<<"---av V"<<endl;
  cout<<total_f_V_move_R1<<endl;
  cout<<constant_fit(total_f_V_move_R1,0,6)<<endl;
  
  cout<<constant_fit(paste(total_f_S_move_R1,total_f_V_move_R1),0,12,"/tmp/unc.xmg")<<endl;
  cout<<correlated_constant_fit(paste(total_f_S_move_R1,total_f_V_move_R1),0,12,"/tmp/corr.xmg")<<endl;
  
  return 0;
}
