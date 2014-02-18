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

void simple_plot(const char *path,jvec in,ios::openmode mode=ios::out)
{
  ofstream out(combine("plots/%s.xmg",path).c_str(),mode);
  out<<"@type xydy"<<endl;
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

void study_single_tsep(int tsep,ios::openmode mode=ios::out)
{
  //fit pion mass
  jvec pion_rest=load("pion-00WW",0).simmetrized(1),pion_rest_alt=load("pion-00WP",0).simmetrized(1);
  simple_plot("pion_rest_effmass",effective_mass(pion_rest));
  simple_plot("pion_rest",pion_rest);
  jack M_Pi=constant_fit(effective_mass(pion_rest_alt),15,T/2-1,"plots/pion_rest_alt.xmg");
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
  jack M_K=constant_fit(effective_mass(kaon_rest_alt),15,T/2-1,"plots/kaon_rest_alt.xmg");  
  jack ZWW_K_rest(njacks),dum_K_rest(njacks);
  two_pts_fit(dum_K_rest,ZWW_K_rest,kaon_rest,15,T/2-2);
  jack ZWP_K_rest(njacks),dum_K_rest_alt(njacks);
  two_pts_fit(dum_K_rest_alt,ZWP_K_rest,kaon_rest_alt,15,T/2-2);
  jack ZW_K_rest=sqrt(ZWW_K_rest),ZP_K_rest=ZWP_K_rest/ZW_K_rest;
  
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
  jvec pp_reno_rest=pp_V0_rest;
  for(int t=0;t<=tsep;t++) pp_reno_rest[t]/=-pion_rest[T/4]*exp(-M_Pi*(tsep-T/4));
  for(int t=tsep+1;t<T;t++) pp_reno_rest[t]/=pion_rest[T/4]*exp(-M_Pi*(T-tsep-T/4));
  simple_plot("pp_reno_rest",pp_reno_rest,mode);

  //load insertion of V0 between kaons at rest (on s)
  jvec kk_V0_rest_a=load("zka-00",icombo(tsep,1));
  jvec kk_reno_rest_a=kk_V0_rest_a;
  for(int t=0;t<=tsep;t++) kk_reno_rest_a[t]/=-kaon_rest[T/4]*exp(-M_K*(tsep-T/4));  
  for(int t=tsep+1;t<T;t++) kk_reno_rest_a[t]/=kaon_rest[T/4]*exp(-M_K*(T-tsep-T/4));
  simple_plot("kk_reno_rest_a",kk_reno_rest_a,mode);

  //load insertion of V0 between kaons at rest (on l)
  jvec kk_V0_rest_b=load("zkb-00",icombo(tsep,1));
  jvec kk_reno_rest_b=kk_V0_rest_b;
  for(int t=0;t<=tsep;t++) kk_reno_rest_b[t]/=-kaon_rest[T/4]*exp(-M_K*(tsep-T/4));  
  for(int t=tsep+1;t<T;t++) kk_reno_rest_b[t]/=kaon_rest[T/4]*exp(-M_K*(T-tsep-T/4));
  simple_plot("kk_reno_rest_b",kk_reno_rest_b,mode);
  
  //load insertion of V0 between pions moving
  jvec pp_V0_move=load("zpa-11",icombo(tsep,1));
  jvec pp_reno_move=pp_V0_move;
  for(int t=0;t<=tsep;t++) pp_reno_move[t]/=-pion_move[T/4]*exp(-E_Pi*(tsep-T/4));
  for(int t=tsep+1;t<T;t++) pp_reno_move[t]/=pion_move[T/4]*exp(-E_Pi*(T-tsep-T/4));
  simple_plot("pp_reno_move",pp_reno_move,mode);
  
  //verify the ratio between pion correlators at rest and in motion
  jvec pp_V0_move_fr_rest=(pp_V0_move/pion_move[tsep])/(pp_V0_rest/pion_rest[tsep]);
  //this should go as the ratio of P0...
  simple_plot("pp_V0_move_fr_rest",pp_V0_move_fr_rest,mode);
  
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
  simple_plot("norm_a",norm_a,mode);
  jvec norm_b=pp_V0_rest*kk_V0_rest_b;
  simple_plot("norm_b",norm_b,mode);
  
  //Zv
  jvec Zv_pion_rest=-pion_rest_tilde[tsep]/pp_V0_rest;
  simple_plot("Zv_pion_rest",Zv_pion_rest.subset(0,tsep),mode);
  jvec Zv_pion_move=-pion_move_tilde[tsep]/pp_V0_move;
  simple_plot("Zv_pion_move",Zv_pion_move.subset(0,tsep),mode);
  jvec Zv_kaon_rest_a=-kaon_rest_tilde[tsep]/kk_V0_rest_a;
  simple_plot("Zv_kaon_rest_a",Zv_kaon_rest_a.subset(0,tsep),mode);
  jvec Zv_kaon_rest_b=-kaon_rest_tilde[tsep]/kk_V0_rest_b;
  simple_plot("Zv_kaon_rest_b",Zv_kaon_rest_b.subset(0,tsep),mode);
  
  //extract Zv from pion and kaon_a at rest
  jack Zv_pion=constant_fit(Zv_pion_rest,tsep/2-tsep/4,tsep/2+tsep/4);
  jack Zv_kaon=constant_fit(Zv_kaon_rest_a,tsep/2-tsep/4,tsep/2+tsep/4);
  cerr<<"ZV: "<<tsep<<" "<<smart_print(Zv_pion)<<" "<<smart_print(Zv_kaon)<<endl;
  
  //////////////////////////////////// kl3 //////////////////////////////////
  
  //load kl3 with insertions
  jvec kp_VK_move=(load("kl3-01",icombo(tsep,0))+0*simmetric(load("kl3-01",icombo(T-tsep,0))))/2*2;
  jvec kp_V0_move=(load("kl3-01",icombo(tsep,1))-simmetric(load("kl3-01",icombo(T-tsep,1))))/2;
  jvec kp_S0_move=load("kl3-01",icombo(tsep,2));
  jvec kp_VK_rest=(load("kl3-00",icombo(tsep,0))+simmetric(load("kl3-00",icombo(T-tsep,0))))/2;
  jvec kp_V0_rest=(load("kl3-00",icombo(tsep,1))-simmetric(load("kl3-00",icombo(T-tsep,1))))/2;
  
  //check that VK is 0
  simple_plot("kp_VK_move",kp_VK_move,mode);
  simple_plot("kp_VK_rest",kp_VK_rest,mode);
  
  //compare moving and rest kl3
  simple_plot("kp_V0_move",kp_V0_move,mode);
  simple_plot("kp_V0_rest",kp_V0_rest,mode);
  
  //plot the plateaux for V0 at rest
  jvec kp_V0_rest_R1=Zv_kaon*sqrt(4*M_Pi*M_K*kp_V0_rest*partially_simmetric(kp_V0_rest,tsep)/
				  (pion_rest_tilde[tsep]*kaon_rest_tilde[tsep]));
  jvec kp_V0_rest_R2=sqrt(4*M_Pi*M_K*kp_V0_rest*partially_simmetric(kp_V0_rest,tsep)/(pp_V0_rest*kk_V0_rest_a));
  jvec kp_V0_rest_R3=-4*sqrt(M_Pi*M_K)*kp_V0_rest/pion_rest_tilde[tsep];
  jvec kp_V0_rest_R4=-kp_V0_rest*Zv_kaon;
  for(int t=0;t<tsep;t++) kp_V0_rest_R3[t]*=sqrt(kaon_rest_tilde[tsep-t]*pion_rest_tilde[t]*pion_rest_tilde[tsep]);
  for(int t=0;t<tsep;t++) kp_V0_rest_R3[t]/=sqrt(pion_rest_tilde[tsep-t]*kaon_rest_tilde[t]*kaon_rest_tilde[tsep]);
  for(int t=0;t<tsep;t++) kp_V0_rest_R4[t]/=ZW_Pi_rest*ZW_K_rest*exp(-M_K*t)*exp(-M_Pi*(tsep-t))/(2*M_Pi*2*M_K);
  simple_plot("kp_V0_rest_R1",kp_V0_rest_R1.subset(0,tsep)/(M_Pi+M_K),mode);
  simple_plot("kp_V0_rest_R2",kp_V0_rest_R2.subset(0,tsep)/(M_Pi+M_K),mode);
  simple_plot("kp_V0_rest_R3",kp_V0_rest_R3.subset(0,tsep)/(M_Pi+M_K),mode);
  simple_plot("kp_V0_rest_R4",kp_V0_rest_R4.subset(0,tsep)/(M_Pi+M_K),mode);
  
  //plot the plateaux for V0 in motion
  jvec kp_V0_move_R1=Zv_kaon*sqrt(4*E_Pi*M_K*kp_V0_move*partially_simmetric(kp_V0_move,tsep)/
				  (pion_move_tilde[tsep]*kaon_rest_tilde[tsep]));
  jvec kp_V0_move_R2=sqrt(4*E_Pi*M_K*kp_V0_move*partially_simmetric(kp_V0_move,tsep)/(pp_V0_move*kk_V0_rest_a));
  jvec kp_V0_move_R4=-kp_V0_move*Zv_kaon;
  for(int t=0;t<tsep;t++) kp_V0_move_R4[t]/=ZW_Pi_move*ZW_K_rest*exp(-M_K*t)*exp(-E_Pi*(tsep-t))/(2*E_Pi*2*M_K);
  simple_plot("kp_V0_move_R1",kp_V0_move_R1.subset(0,tsep)/(E_Pi+M_K),mode);
  simple_plot("kp_V0_move_R2",kp_V0_move_R2.subset(0,tsep)/(E_Pi+M_K),mode);
  simple_plot("kp_V0_move_R4",kp_V0_move_R4.subset(0,tsep)/(E_Pi+M_K),mode);

  //plot the plateaux for VK in motion
  jvec kp_VK_move_R1=Zv_kaon*sqrt(4*E_Pi*M_K*kp_VK_move*partially_simmetric(kp_VK_move,tsep)/
				  (pion_move_tilde[tsep]*kaon_rest_tilde[tsep]));
  jvec kp_VK_move_R2=sqrt(4*E_Pi*M_K*kp_VK_move*partially_simmetric(kp_VK_move,tsep)/(pp_V0_move*kk_V0_rest_a));
  jvec kp_VK_move_R4=-kp_VK_move*Zv_kaon;
  for(int t=0;t<tsep;t++) kp_VK_move_R4[t]/=ZW_Pi_move*ZW_K_rest*exp(-M_K*t)*exp(-E_Pi*(tsep-t))/(2*E_Pi*2*M_K);
  simple_plot("kp_VK_move_R1",kp_VK_move_R1.subset(0,tsep)/mom_i,mode);
  simple_plot("kp_VK_move_R2",kp_VK_move_R2.subset(0,tsep)/mom_i,mode);
  simple_plot("kp_VK_move_R4",kp_VK_move_R4.subset(0,tsep)/mom_i,mode);
  
  //solve the ratios
  jvec f_R1_corr=(kp_V0_move_R1*QK-Q0*kp_VK_move_R1)/delta;
  jvec f_R2_corr=(kp_V0_move_R2*QK-Q0*kp_VK_move_R2)/delta;
  jvec f_R4_corr=(kp_V0_move_R4*QK-Q0*kp_VK_move_R4)/delta;
  simple_plot("kp_f_move_R1_solve",f_R1_corr.subset(0,tsep),mode);
  simple_plot("kp_f_move_R2_solve",f_R2_corr.subset(0,tsep),mode);
  simple_plot("kp_f_move_R4_solve",f_R4_corr.subset(0,tsep),mode);
  
  //plot the plateaux for S0 in motion
  jvec kp_S0_move_R4=kp_S0_move;
  for(int t=0;t<tsep;t++) kp_S0_move_R4[t]/=ZW_Pi_move*ZW_K_rest*exp(-M_K*t)*exp(-E_Pi*(tsep-t))/(2*E_Pi*2*M_K);
  simple_plot("kp_S0_move_R4",kp_S0_move_R4.subset(0,tsep)*(M_K*M_K-M_Pi*M_Pi)/(ms-mu),mode);
  
  //fit solve ratio 4            0     4     8     12    16     20     24     28      32      36      40
  int f_R4_tint[11][2]={{0,0},{1,3},{3,5},{5,7},{6,10},{9,11},{6,18},{12,16},{14,19},{12,30},{13,33}};
  jack f_R4=constant_fit(f_R4_corr.subset(0,tsep),
			 f_R4_tint[tsep/4][0],f_R4_tint[tsep/4][1],combine("plots/f_R4_%02d.xmg",tsep).c_str());
  ofstream f_R4_out("plots/f_R4.xmg",mode);
  f_R4_out<<tsep<<" "<<f_R4<<endl;
  f_R4_out.close();
  
  //fit VK ratio 1            0     4     8     12    16     20     24     28      32      36      40
  int f0_R1_VK_tint[11][2]={{0,0},{1,3},{3,5},{5,7},{6,10},{9,11},{6,18},{12,16},{14,19},{12,30},{13,33}};
  jack f0_R1_VK=constant_fit(kp_VK_move_R1.subset(0,tsep)/mom_i,
	     f0_R1_VK_tint[tsep/4][0],f0_R1_VK_tint[tsep/4][1],combine("plots/kp_VK_move_R1_%02d.xmg",tsep).c_str());
  ofstream f0_R1_VK_out("plots/f0_R1_VK.xmg",mode);
  f0_R1_VK_out<<tsep<<" "<<f0_R1_VK<<endl;
  f0_R1_VK_out.close();

  //fit VK ratio 2            0     4     8     12    16     20     24     28      32      36      40
  int f0_R2_VK_tint[11][2]={{0,0},{1,3},{4,6},{4,6},{5,10},{8,11},{9,18},{6,18},{17,22},{12,30},{13,33}};
  jack f0_R2_VK=constant_fit(kp_VK_move_R2.subset(0,tsep)/mom_i,
	     f0_R2_VK_tint[tsep/4][0],f0_R2_VK_tint[tsep/4][1],combine("plots/kp_VK_move_R2_%02d.xmg",tsep).c_str());
  ofstream f0_R2_VK_out("plots/f0_R2_VK.xmg",mode);
  f0_R2_VK_out<<tsep<<" "<<f0_R2_VK<<endl;
  f0_R2_VK_out.close();

  //fit VK ratio 4            0     4     8     12    16     20     24     28      32      36      40
  int f0_R4_VK_tint[11][2]={{0,0},{1,3},{4,6},{4,6},{5,10},{8,11},{9,18},{6,18},{17,22},{12,30},{13,33}};
  jack f0_R4_VK=constant_fit(kp_VK_move_R4.subset(0,tsep)/mom_i,
	     f0_R4_VK_tint[tsep/4][0],f0_R4_VK_tint[tsep/4][1],combine("plots/kp_VK_move_R4_%02d.xmg",tsep).c_str());
  ofstream f0_R4_VK_out("plots/f0_R4_VK.xmg",mode);
  f0_R4_VK_out<<tsep<<" "<<f0_R4_VK<<endl;
  f0_R4_VK_out.close();

  //fit V0 ratio 1            0     4     8     12    16     20     24     28      32      36      40
  int f0_R1_V0_tint[11][2]={{0,0},{1,3},{4,6},{4,8},{5,11},{5,14},{6,20},{4,24},{14,17},{8,28},{10,35}};
  jack f0_R1_V0=constant_fit(kp_V0_move_R1.subset(0,tsep)/(E_Pi+M_K),
	     f0_R1_V0_tint[tsep/4][0],f0_R1_V0_tint[tsep/4][1],combine("plots/kp_V0_move_R1_%02d.xmg",tsep).c_str());
  ofstream f0_R1_V0_out("plots/f0_R1_V0.xmg",mode);
  f0_R1_V0_out<<tsep<<" "<<f0_R1_V0<<endl;
  f0_R1_V0_out.close();

  //fit V0 ratio 4            0     4     8     12    16     20     24     28      32      36      40
  int f0_R4_V0_tint[11][2]={{0,0},{1,3},{4,6},{7,8},{8,10},{9,13},{9,18},{8,22},{17,22},{12,31},{13,33}};
  jack f0_R4_V0=constant_fit(kp_V0_move_R4.subset(0,tsep)/(E_Pi+M_K),
	     f0_R4_V0_tint[tsep/4][0],f0_R4_V0_tint[tsep/4][1],combine("plots/kp_V0_move_R4_%02d.xmg",tsep).c_str());
  ofstream f0_R4_V0_out("plots/f0_R4_V0.xmg",mode);
  f0_R4_V0_out<<tsep<<" "<<f0_R4_V0<<endl;
  f0_R4_V0_out.close();
  
  ///////////////////////////
  
  //fit S0 ratio 4            0     4     8     12    16     20     24     28      32      36      40
  int f0_R4_S0_tint[11][2]={{0,0},{1,3},{4,6},{7,8},{8,10},{9,13},{8,19},{8,23},{10,23},{12,24},{13,33}};
  jack f0_R4_S0=constant_fit(kp_S0_move_R4.subset(0,tsep)*(ms-mu)/(M_K*M_K-M_Pi*M_Pi),
	     f0_R4_S0_tint[tsep/4][0],f0_R4_S0_tint[tsep/4][1],combine("plots/kp_S0_move_R4_%02d.xmg",tsep).c_str());
  ofstream f0_R4_S0_out("plots/f0_R4_S0.xmg",mode);
  f0_R4_S0_out<<tsep<<" "<<f0_R4_S0<<endl;
  f0_R4_S0_out.close();
}

int main()
{
  study_single_tsep(4);
  for(int tsep=8;tsep<=40;tsep+=4) study_single_tsep(tsep,ios::out|ios::app);
  
  return 0;
}
