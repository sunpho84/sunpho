#include "include.h"

const int T=96,L=48,njacks=81;
const double a=1/1.75,mom_i=0.589278*2*M_PI/L;
const double mu=0.00078,ms=0.0362;

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

void study_single_tsep(int tsep,ios::openmode mode=ios::out)
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
  jvec pp_V0_rest=load("zpa-00",icombo(tsep,1));
  simple_plot("pp_V0_rest",pp_V0_rest,tsep,mode);
  jvec pp_reno_rest=pp_V0_rest;
  for(int t=0;t<=tsep;t++) pp_reno_rest[t]/=-pion_rest[T/4]*exp(-M_Pi*(tsep-T/4));
  for(int t=tsep+1;t<T;t++) pp_reno_rest[t]/=pion_rest[T/4]*exp(-M_Pi*(T-tsep-T/4));
  simple_plot("pp_reno_rest",pp_reno_rest,tsep,mode);

  //load insertion of V0 between kaons at rest (on s)
  jvec kk_V0_rest_a=load("zka-00",icombo(tsep,1));
  jvec kk_reno_rest_a=kk_V0_rest_a;
  for(int t=0;t<=tsep;t++) kk_reno_rest_a[t]/=-kaon_rest[T/4]*exp(-M_K*(tsep-T/4));  
  for(int t=tsep+1;t<T;t++) kk_reno_rest_a[t]/=kaon_rest[T/4]*exp(-M_K*(T-tsep-T/4));
  simple_plot("kk_reno_rest_a",kk_reno_rest_a,tsep,mode);

  //load insertion of V0 between kaons at rest (on l)
  jvec kk_V0_rest_b=load("zkb-00",icombo(tsep,1));
  jvec kk_reno_rest_b=kk_V0_rest_b;
  for(int t=0;t<=tsep;t++) kk_reno_rest_b[t]/=-kaon_rest[T/4]*exp(-M_K*(tsep-T/4));  
  for(int t=tsep+1;t<T;t++) kk_reno_rest_b[t]/=kaon_rest[T/4]*exp(-M_K*(T-tsep-T/4));
  simple_plot("kk_reno_rest_b",kk_reno_rest_b,tsep,mode);
  
  //load insertion of V0 between pions moving
  jvec pp_V0_move=load("zpa-11",icombo(tsep,1));
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
  jvec Zv_pion_rest=-pion_rest_tilde[tsep]/pp_V0_rest;
  simple_plot("Zv_pion_rest",-1/Zv_pion_rest.subset(0,tsep),tsep,mode);
  jvec Zv_pion_move=-pion_move_tilde[tsep]/pp_V0_move;
  simple_plot("Zv_pion_move",Zv_pion_move.subset(0,tsep),tsep,mode);
  jvec Zv_kaon_rest_a=-kaon_rest_tilde[tsep]/kk_V0_rest_a;
  simple_plot("Zv_kaon_rest_a",Zv_kaon_rest_a.subset(0,tsep),tsep,mode);
  jvec Zv_kaon_rest_b=-kaon_rest_tilde[tsep]/kk_V0_rest_b;
  simple_plot("Zv_kaon_rest_b",Zv_kaon_rest_b.subset(0,tsep),tsep,mode);
  
  //extract Zv from pion and kaon_a at rest
  jack Zv_pion=constant_fit(Zv_pion_rest,tsep/2-tsep/4,tsep/2+tsep/4);
  jack Zv_kaon=constant_fit(Zv_kaon_rest_a,tsep/2-tsep/4,tsep/2+tsep/4);
  jack Zv_mixed=constant_fit(sqrt(Zv_pion_rest*Zv_kaon_rest_a),tsep/2-tsep/4,tsep/2+tsep/4);
  cerr<<"ZV: "<<tsep<<" "<<smart_print(Zv_pion)<<" "<<smart_print(Zv_kaon)<<endl;
  cerr<<"ZV_mixed: "<<tsep<<" "<<smart_print(Zv_mixed)<<" "<<smart_print(sqrt(Zv_kaon*Zv_pion))<<
    "in the middle: "<<sqrt(Zv_pion_rest*Zv_kaon_rest_a)[tsep/2]<<endl;
  
  //////////////////////////////////// kl3 //////////////////////////////////
  
  //load kl3 with insertions
  jvec kp_VK_move=(load("kl3-01",icombo(tsep,0))+simmetric(load("kl3-01",icombo(T-tsep,0))))/2;
  jvec kp_V0_move=(load("kl3-01",icombo(tsep,1))-simmetric(load("kl3-01",icombo(T-tsep,1))))/2;
  jvec kp_S0_move=load("kl3-01",icombo(tsep,2));
  jvec kp_VK_rest=(load("kl3-00",icombo(tsep,0))+simmetric(load("kl3-00",icombo(T-tsep,0))))/2;
  jvec kp_V0_rest=(load("kl3-00",icombo(tsep,1))-simmetric(load("kl3-00",icombo(T-tsep,1))))/2;
  
  //check that VK is 0
  simple_plot("kp_VK_move",kp_VK_move,tsep,mode);
  simple_plot("kp_VK_rest",kp_VK_rest,tsep,mode);
  
  //compare moving and rest kl3
  simple_plot("kp_V0_move",kp_V0_move,tsep,mode);
  simple_plot("kp_V0_rest",kp_V0_rest,tsep,mode);
  
  //solve the correlator for V and for S
  jvec kp_V_rest=kp_V0_move/(M_Pi+M_K);
  jvec kp_V_move=(kp_V0_move*QK-Q0*kp_VK_move)/delta;
  jvec kp_S_move=kp_S0_move*(ms-mu)/(M_K*M_K-M_Pi*M_Pi);
  
  ////////////////////// * ANDREAS comparison * ///////////////////////////
  
  //build ratios
  jvec kp_V0_rest_R1=Zv_mixed*sqrt(4*M_Pi*M_K*kp_V0_rest*partially_simmetric(kp_V0_rest,tsep)/
				  (pion_rest_tilde[tsep]*kaon_rest_tilde[tsep]));
  jvec kp_V0_move_R1=Zv_mixed*sqrt(4*E_Pi*M_K*kp_V0_move*partially_simmetric(kp_V0_move,tsep)/
				 (pion_move_tilde[tsep]*kaon_rest_tilde[tsep]));
  jvec kp_VK_move_R1=Zv_mixed*sqrt(4*E_Pi*M_K*kp_VK_move*partially_simmetric(kp_VK_move,tsep)/
				 (pion_move_tilde[tsep]*kaon_rest_tilde[tsep]));
  
  
  kp_VK_move_R1=kp_VK_move*partially_simmetric(kp_VK_move,tsep)/
    (pion_move_tilde[tsep]*kaon_rest_tilde[tsep]);
  
  simple_plot("kp_V0_rest_R1",kp_V0_rest_R1.subset(0,tsep),tsep,mode);
  simple_plot("kp_V0_move_R1",kp_V0_move_R1.subset(0,tsep),tsep,mode);
  simple_plot("kp_VK_move_R1",kp_VK_move_R1.subset(0,tsep),tsep,mode);
  
  //fit V move ratios             0     4     8     12    16     20     24     28     32      36     40
  int f_VK_move_R1_tint[11][2]={{0,0},{1,3},{4,6},{4,8},{5,11},{6,14},{6,20},{4,24},{14,17},{8,28},{10,35}};
  jack f_VK_move_R1=fit_ratio("f_VK_move_R1",kp_VK_move_R1,f_VK_move_R1_tint[it],tsep,mode);
  cout<<"VK tsep "<<tsep<<" at timeslice 4: "<<kp_VK_move_R1[4]<<endl;

  /////////////////////////////////////////////////////////////////////////

  //build ratios
  jvec kp_V_rest_R1=Zv_mixed*sqrt(4*M_Pi*M_K*kp_V_rest*partially_simmetric(kp_V_rest,tsep)/
				 (pion_rest_tilde[tsep]*kaon_rest_tilde[tsep]));
  jvec kp_V_move_R1=Zv_mixed*sqrt(4*E_Pi*M_K*kp_V_move*partially_simmetric(kp_V_move,tsep)/
				 (pion_move_tilde[tsep]*kaon_rest_tilde[tsep]));
  cout<<"S0 tsep "<<tsep<<" at timeslice 4: "<<kp_S0_move[4]<<endl;
  //jvec kp_S_move_R1=kp_S0_move*partially_simmetric(kp_S0_move,tsep)/(pion_move_tilde[tsep]*kaon_rest_tilde[tsep]);
  jvec kp_S_move_R1=sqrt(4*E_Pi*M_K*kp_S_move*partially_simmetric(kp_S_move,tsep)/
			 (pion_move_tilde[tsep]*kaon_rest_tilde[tsep]));
  cout<<"S0 tsep "<<tsep<<" at timeslice 4: "<<kp_S_move_R1[4]<<endl;
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
  int f_V_rest_R1_tint[11][2]={{0,0},{1,3},{4,6},{4,8},{5,11},{5,14},{6,20},{4,24},{14,17},{8,28},{10,35}};
  int f_V_rest_R2_tint[11][2]={{0,0},{1,3},{4,6},{4,8},{5,11},{5,14},{6,20},{4,24},{14,17},{8,28},{10,35}};
  int f_V_rest_R3_tint[11][2]={{0,0},{1,3},{4,6},{4,8},{5,11},{5,14},{6,20},{4,24},{14,17},{8,28},{10,35}};
  int f_V_rest_R4_tint[11][2]={{0,0},{1,3},{4,6},{4,8},{5,11},{5,14},{6,20},{4,24},{14,17},{8,28},{10,35}};
  jack f_V_rest_R1=fit_ratio("f_V_rest_R1",kp_V_rest_R1,f_V_rest_R1_tint[it],tsep,mode);
  jack f_V_rest_R2=fit_ratio("f_V_rest_R2",kp_V_rest_R2,f_V_rest_R2_tint[it],tsep,mode);
  jack f_V_rest_R3=fit_ratio("f_V_rest_R3",kp_V_rest_R3,f_V_rest_R3_tint[it],tsep,mode);
  jack f_V_rest_R4=fit_ratio("f_V_rest_R4",kp_V_rest_R4,f_V_rest_R4_tint[it],tsep,mode);

  //fit V move ratios            0     4     8     12    16     20     24     28      32      36      40
  int f_V_move_R1_tint[11][2]={{0,0},{1,3},{4,6},{7,8},{7,9},{8,12},{9,15},{8,22},{17,22},{12,31},{13,33}};
  int f_V_move_R2_tint[11][2]={{0,0},{1,3},{4,6},{7,8},{8,10},{9,13},{9,18},{8,22},{17,22},{12,31},{13,33}};
  int f_V_move_R3_tint[11][2]={{0,0},{1,3},{4,6},{7,8},{8,10},{9,13},{9,18},{8,22},{17,22},{12,31},{13,33}};
  int f_V_move_R4_tint[11][2]={{0,0},{1,3},{4,6},{7,8},{8,10},{9,13},{9,18},{8,22},{17,22},{12,31},{13,33}};
  jack f_V_move_R1=fit_ratio("f_V_move_R1",kp_V_move_R1,f_V_move_R1_tint[it],tsep,mode);
  jack f_V_move_R2=fit_ratio("f_V_move_R2",kp_V_move_R2,f_V_move_R2_tint[it],tsep,mode);
  jack f_V_move_R3=fit_ratio("f_V_move_R3",kp_V_move_R3,f_V_move_R3_tint[it],tsep,mode);
  jack f_V_move_R4=fit_ratio("f_V_move_R4",kp_V_move_R4,f_V_move_R4_tint[it],tsep,mode);
  
  //fit S move ratios            0     4     8     12    16    20     24     28     32       36      40
  int f_S_move_R1_tint[11][2]={{0,0},{1,3},{4,6},{7,8},{8,9},{8,12},{7,17},{9,19},{11,21},{12,24},{13,31}};
  //int f_S_move_R2_tint[11][2]={{0,0},{1,3},{4,6},{7,8},{8,10},{9,13},{8,19},{8,23},{10,23},{12,24},{13,33}};
  //int f_S_move_R3_tint[11][2]={{0,0},{1,3},{4,6},{7,8},{8,10},{9,13},{8,19},{8,23},{10,23},{12,24},{13,33}};
  int f_S_move_R4_tint[11][2]={{0,0},{1,3},{4,6},{7,8},{8,10},{9,13},{8,19},{8,23},{10,23},{12,24},{13,33}};
  jack f_S_move_R1=fit_ratio("f_S_move_R1",kp_S_move_R1,f_S_move_R1_tint[it],tsep,mode);
  //jack f_S_move_R2=fit_ratio("f_S_move_R2",kp_S_move_R2,f_S_move_R2_tint[it],tsep,mode);
  //jack f_S_move_R3=fit_ratio("f_S_move_R3",kp_S_move_R3,f_S_move_R3_tint[it],tsep,mode);
  jack f_S_move_R4=fit_ratio("f_S_move_R4",kp_S_move_R4,f_S_move_R4_tint[it],tsep,mode);
}

int main()
{
  study_single_tsep(4);
  for(int tsep=8;tsep<=40;tsep+=4) study_single_tsep(tsep,ios::out|ios::app);
  
  return 0;
}
