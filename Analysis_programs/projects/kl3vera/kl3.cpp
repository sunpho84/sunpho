#include "include.h"

#include <vector>
#include <omp.h>

#include "loc.hpp"

//mkdir vera
//sshfs edinbgq:/work/dp008/dp008/DWF/2+1f/48nt96/IWASAKI/ls12/M1.8/b2.37/ms0.019/ml0.002/ukrun_kl3_dec_2015/results vera

using namespace std;

const int T=96,L=48;
const int ntot_confs=82;
const int clust_size=4;
const int njacks=ntot_confs/clust_size;
const int nconfs=clust_size*njacks;
const int corr_size=2*T*16*16;
const size_t read_size=sizeof(double)*corr_size;
const size_t padding_size=44;
const int nsep=3;
const int tseps[nsep]={27,35,39};
const int nmel=3;
//--
const int nfl_combo_3pts=2;
const char m_tag_3pts[nfl_combo_3pts][20]={"s_0.02144","l_0.002144"};
const char fly_tag_3pts[nfl_combo_3pts][100]={"0_0_0","0.65123416_0.65123416_0.65123416"};
const int i3pts_KK=0,i3pts_KPi=1;
const char templ_3pts_path[]="vera/%d/meson/meson_m%s_Z2Wall_LL_0.0_t%d_twist_%s_ml_0.002144_Z2Wall_LL_0.0_t%d_twist_0_0_0_seq_ms_0.02144_Z2Wall_LL_0.0_t%d_twist_0_0_0.%d.bin";
const int iS0=0,iVi=1,iV0=2;
//--
const int nfl_combo_2pts=3;
const char m1_tag_2pts[nfl_combo_2pts][20]={"l_0.002144","l_0.002144","s_0.02144"};
const char m2_tag_2pts[nfl_combo_2pts][20]={"s_0.02144","l_0.002144","s_0.02144"};
const char fly_tag_2pts[nfl_combo_2pts][100]={"0_0_0","0_0_0"};
const int i2pts_K=0,i2pts_Pi_rest=1,i2pts_SS_rest=2;
const char templ_2pts_path[]="vera/%d/meson/meson_m%s_Z2Wall_LL_0.0_t%d_twist_0_0_0_m%s_Z2Wall_LL_0.0_t%d_twist_0_0_0.%d.bin";
//--
const int nmel_tot=nfl_combo_3pts*nmel*nsep;
const int ntpts_tot=nfl_combo_2pts;
const char loc_tpts_path[]="tpts";
const char loc_mels_path[]="mels";
const double theta=0.65123416;
const double pi=theta*M_PI/L;
const double p2=3*sqr(pi);
const double a=1/2.77;
jack P0,Q0,delta,Q2;
double PK=pi,QK=-pi;

//data
vector<jvec> mels(nmel_tot,jvec(T,njacks));
vector<jvec> tpts(ntpts_tot,jvec(T,njacks));

int big_endian;

void change_endianness(char *in,size_t n)
{
  for(size_t i=0;i<n/2;i++)
    std::swap(in[i],in[n-1-i]);
}

void change_endianness(vector<jvec> &v)
{for(size_t i=0;i<v.size();i++) for(int t=0;t<T;t++) for(int c=0;c<=njacks;c++) change_endianness((char*)&(v[i][t][c]),sizeof(double));}

int ind_data(int t,int ig_so,int ig_si,int ri)
{return ri+2*(ig_so+16*(ig_si+16*t));}

int ind_mel(int ifl_combo,int isep,int imel)
{return imel+nmel*(isep+nsep*ifl_combo);}

void add_3pts(vector<jvec> &out,int ifl_combo,int isep,int imel,int iconf,double *in,int ig_si,int ri)
{
  const int ig_so=15;
  for(int t=0;t<T;t++)
    out[ind_mel(ifl_combo,isep,imel)][t][iconf/clust_size]+=in[ind_data(t,ig_so,ig_si,ri)];
}

void add_2pts(vector<jvec> &out,int ifl_combo,int iconf,double *in)
{
  const int ig_so=15,ig_si=15,ri=0;
  for(int t=0;t<T;t++)
      out[ifl_combo][t][iconf/clust_size]+=in[ind_data(t,ig_so,ig_si,ri)];
}

void read_data(double *out,const char *path)
{
  //cerr<<"Opening: "<<path<<endl;
  //open
  ifstream in(path);
  if(!in.good()) crash("opening %s",path);
  
  //read
  in.ignore(padding_size);
  in.read((char*)out,read_size);
  in.close();
}

void read_3pts(vector<jvec> &out,const char *path,int ifl_combo,int isep,int iconf)
{
  double buf[read_size];
  read_data(buf,path);
  
  //change endianness
  if(big_endian) for(int i=0;i<corr_size;i++) change_endianness((char*)&(buf[i]),sizeof(double));
  
  //add
  const int nmel_in=5;
  const int ig_si[nmel_in]={0,1,2,4,8};
  const int ri[nmel_in]={0,1,1,1,0};
  const int imel_out[5]={0,1,1,1,2};
  for(int imel_in=0;imel_in<nmel_in;imel_in++)
    add_3pts(out,ifl_combo,isep,imel_out[imel_in],iconf,buf,ig_si[imel_in],ri[imel_in]);
}

void read_2pts(vector<jvec> &out,const char *path,int ifl_combo,int iconf)
{
  double buf[read_size];
  read_data(buf,path);
  
  //change endianness
  if(big_endian) for(int i=0;i<corr_size;i++) change_endianness((char*)&(buf[i]),sizeof(double));
  add_2pts(out,ifl_combo,iconf,buf);
}

void load_from_remote()
{
  //put to zero
  for(size_t i=0;i<mels.size();i++) mels[i]=0;
  for(size_t i=0;i<tpts.size();i++) tpts[i]=0;
  
  //load
#pragma omp parallel for
  for(int iconf=0;iconf<nconfs;iconf++)
    {
      cout<<omp_get_thread_num()<<" "<<iconf<<"/"<<nconfs<<endl;
      
      int conf_name=4000+40*iconf;
      
	for(int t0=0;t0<T;t0+=2)
	  {
	    for(int isep=0;isep<nsep;isep++)
	      for(int ifl_combo=0;ifl_combo<nfl_combo_3pts;ifl_combo++)
		read_3pts(mels,combine(templ_3pts_path,conf_name,m_tag_3pts[ifl_combo],t0,fly_tag_3pts[ifl_combo],t0,(t0+tseps[isep])%T,conf_name).c_str(),ifl_combo,isep,iconf);
	    for(int ifl_combo=0;ifl_combo<nfl_combo_2pts;ifl_combo++)
	      read_2pts(tpts,combine(templ_2pts_path,conf_name,m1_tag_2pts[ifl_combo],t0,m2_tag_2pts[ifl_combo],t0,conf_name).c_str(),ifl_combo,iconf);
	  }
    }
  
  //clusterize
  for(size_t i=0;i<mels.size();i++)
    {
      mels[i].clusterize(clust_size);
      mels[i]/=L*L*L*T/2;
    }
  for(size_t i=0;i<tpts.size();i++)
    {
      tpts[i].clusterize(clust_size);
      tpts[i]/=L*L*L*T/2;
    }
}

void save_locally()
{
  FILE *tpts_out=open_file(loc_tpts_path,"w");
  if(!big_endian) change_endianness(tpts);
  for(size_t i=0;i<tpts.size();i++) tpts[i].write_to_binfile(tpts_out);
  fclose(tpts_out);
  
  FILE *mels_out=open_file(loc_mels_path,"w");
  if(!big_endian) change_endianness(mels);
  for(size_t i=0;i<mels.size();i++) mels[i].write_to_binfile(mels_out);
  fclose(mels_out);
}

void load_locally()
{
  cerr<<"Using local path"<<endl;
  
  FILE *tpts_out=open_file(loc_tpts_path,"r");
  for(size_t i=0;i<tpts.size();i++) tpts[i].load(tpts_out);
  if(big_endian) change_endianness(tpts);
  fclose(tpts_out);
  
  FILE *mels_out=open_file(loc_mels_path,"r");
  for(size_t i=0;i<mels.size();i++) mels[i].load(mels_out);
  if(big_endian) change_endianness(mels);
  fclose(mels_out);
}

jvec tilded_corr(jvec corr,jack aM)
{
  jvec corr_tilde(T,njacks);
  for(int t=0;t<=T/2;t++) corr_tilde[t]=corr[t]-corr[T/2]*exp(-aM*(T/2-t))/2;
  for(int t=T/2;t<T;t++)  corr_tilde[t]=corr[T/2]*exp(-aM*(t-T/2))/2;
  
  return corr_tilde;
}

void print_in_01(ofstream &fout,jvec out)
{
  fout<<"@type xydy"<<endl;
  for(int i=0;i<out.nel;i++) fout<<(double)i/(out.nel-1)<<" "<<out[i]<<endl;
  fout<<"&"<<endl;
}

#ifdef ANALYSIS

void print_in_01(string path,jvec out)
{
  ofstream fout(path);
  print_in_01(fout,out);
}

void study_ratios(jack &fp,jack &fm,jvec V0_corr,jvec VK_corr,jvec norm,int isep,double tmin_0_fr_tsep,double tmax_0_fr_tsep,
		  double tmin_K_fr_tsep,double tmax_K_fr_tsep,int irat)
{
  jvec V0_rat=V0_corr/norm;
  jvec VK_rat=VK_corr/norm;
  
  ofstream V0_out(combine("plots/KPi_V0_sep%d_R%d.xmg",isep+1,irat));
  ofstream VK_out(combine("plots/KPi_VK_sep%d_R%d.xmg",isep+1,irat));
  
  int tmin_0=tseps[isep]*tmin_0_fr_tsep;
  int tmax_0=tseps[isep]*tmax_0_fr_tsep;
  int tmin_K=tseps[isep]*tmin_K_fr_tsep;
  int tmax_K=tseps[isep]*tmax_K_fr_tsep;
  
  jack V0=constant_fit(V0_rat,tmin_0,tmax_0);
  jack VK=constant_fit(VK_rat,tmin_K,tmax_K);
  
  V0_out<<write_constant_with_error(V0,tmin_0_fr_tsep,tmax_0_fr_tsep);
  VK_out<<write_constant_with_error(VK,tmin_K_fr_tsep,tmax_K_fr_tsep);
  
  print_in_01(V0_out,V0_rat);
  print_in_01(VK_out,VK_rat);
  
  fp=(V0*QK-Q0*VK)/delta;
  fm=(V0*PK-P0*VK)/delta;
}

#endif

int main()
{
  check_endianess();
  
  if(file_exists(loc_mels_path)&&file_exists(loc_tpts_path)) load_locally();
  else
    {
      load_from_remote();
      save_locally();
    }
  
#ifdef ANALYSIS

  //take resting Pi
  jvec K_rest=tpts[i2pts_K];
  jvec Pi_rest=tpts[i2pts_Pi_rest];
  jvec SS_rest=tpts[i2pts_SS_rest];
  
  int tmin_K=19;
  jack aM_K(njacks);
  jack Z2_K(njacks);
  two_pts_fit(aM_K, Z2_K,K_rest.simmetrized(1),tmin_K,T/2+1,"plots/K_meff.xmg");
  jack Z_K=sqrt(Z2_K);
  
  int tmin_Pi=19;
  jack aM_Pi(njacks);
  jack Z2_Pi(njacks);
  two_pts_fit(aM_Pi, Z2_Pi,Pi_rest.simmetrized(1),tmin_Pi,T/2+1,"plots/Pi_meff.xmg");
  jack Z_Pi=sqrt(Z2_Pi);
  
  int tmin_SS=19;
  jack aM_SS(njacks);
  jack Z2_SS(njacks);
  two_pts_fit(aM_SS, Z2_SS,SS_rest.simmetrized(1),tmin_SS,T/2+1,"plots/SS_meff.xmg");
  
  jack aE_Pi=sqrt(sqr(aM_Pi)+p2);
  P0=aM_K+aE_Pi;
  Q0=aM_K-aE_Pi;
  Q2=sqr(Q0)-p2;
  delta=P0*QK-Q0*PK;
  
  cout<<"aM_Pi: "<<smart_print(aM_Pi)<<endl;
  cout<<"aE_Pi: "<<smart_print(aE_Pi)<<endl;
  cout<<"aM_K: "<<smart_print(aM_K)<<endl;
  cout<<"aM_SS: "<<smart_print(aM_SS)<<endl;
  cout<<"---------"<<endl;
  cout<<"M_Pi: "<<smart_print(aM_Pi/a)<<endl;
  cout<<"E_Pi: "<<smart_print(aE_Pi/a)<<endl;
  cout<<"M_K: "<<smart_print(aM_K/a)<<endl;
  cout<<"M_SS: "<<smart_print(aM_SS/a)<<endl;
  cout<<"M_SS_chpt: "<<smart_print(sqrt(2*sqr(aM_K)-sqr(aM_Pi))/a)<<endl;
  cout<<"----"<<endl;
  cout<<"Q2: "<<smart_print(Q2)<<endl;
  
  //build the tilded correlation functions
  jvec Pi_rest_tilde=tilded_corr(Pi_rest,aM_Pi);
  jvec K_rest_tilde=tilded_corr(K_rest,aM_K);
  jvec Pi_move_tilde(T,njacks);
  for(int t=0;t<T;t++) Pi_move_tilde[t]=Pi_rest_tilde[t]*exp((aM_Pi-aE_Pi)*t)*aM_Pi/aE_Pi;
  
  //study effective mass
  ofstream stu_eff_V0("plots/3pts_meff_V0.xmg");
  ofstream stu_eff_VK("plots/3pts_meff_VK.xmg");
  stu_eff_V0<<write_line_with_error(aM_K-aE_Pi,aM_K*0,0,1)<<"@type xy\n0.48 0\n0.48 "<<(aM_K-aE_Pi).med()*1.5<<"\n0.52 "<<(aM_K-aE_Pi).med()*1.5<<"\n0.52 0\n&\n";
  stu_eff_VK<<write_line_with_error(aM_K-aE_Pi,aM_K*0,0,1)<<"@type xy\n0.48 0\n0.48 "<<(aM_K-aE_Pi).med()*1.5<<"\n0.52 "<<(aM_K-aE_Pi).med()*1.5<<"\n0.52 0\n&\n";
  
  jvec Zv(nsep,njacks);
  for(int isep=0;isep<nsep;isep++)
    {
      int tsep=tseps[isep];
      
      //renormalization constant
      jvec KK_V0_rest=mels[ind_mel(i3pts_KK,isep,iV0)];
      jvec Zv_K_corr=(-K_rest_tilde[tsep]/KK_V0_rest).subset(0,tsep+1);
      print_in_01(combine("plots/Zv_sep%d.xmg",isep),Zv_K_corr);
      cout<<"Zv["<<isep<<"]: "<<smart_print(Zv_K_corr[tsep/2])<<endl;
      Zv[isep]=Zv_K_corr[tsep/2];
      
      //take corr functions
      jvec KPi_V0_move=(-Zv[isep]*mels[ind_mel(i3pts_KPi,isep,iV0)]).subset(0,tsep+1);
      jvec KPi_VK_move=(-Zv[isep]*mels[ind_mel(i3pts_KPi,isep,iVi)]/3).subset(0,tsep+1);
      
      //define normalization
      jvec KPi_norm_R1(tsep+1,njacks);
      jvec KPi_norm_R2(tsep+1,njacks);
      jvec KPi_norm_R3(tsep+1,njacks);
      for(int t=0;t<=tsep;t++)
	{
	  KPi_norm_R1[t]=sqrt((Pi_move_tilde[tsep]*K_rest_tilde[tsep])/(4*aE_Pi*aM_K));
	  KPi_norm_R2[t]=Z_Pi*Z_K*exp(-t*aE_Pi-(tsep-t)*aM_K)/(2*aM_K*2*aE_Pi);
	  KPi_norm_R3[t]=Pi_move_tilde[t]*K_rest_tilde[tsep-t]/(Z_Pi*Z_K);
	}
      
      //compute ratios
      jvec fp(3,njacks),fm(3,njacks);
      study_ratios(fp[0],fm[0],sqrt(KPi_V0_move*KPi_V0_move.inverted()),sqrt(KPi_VK_move*KPi_VK_move.inverted()),KPi_norm_R1,isep,0.49,0.51,0.49,0.51,1);
      study_ratios(fp[1],fm[1],KPi_V0_move,KPi_VK_move,KPi_norm_R2,isep,0.4,0.47,0.4,0.47,2);
      study_ratios(fp[2],fm[2],KPi_V0_move,KPi_VK_move,KPi_norm_R3,isep,0.4,0.50,0.38,0.45,3);
      
      //print 3pts eff mass
      print_in_01(stu_eff_V0,-aperiodic_effective_mass(KPi_V0_move));
      print_in_01(stu_eff_VK,-aperiodic_effective_mass(KPi_VK_move));
      
      for(int irat=0;irat<3;irat++) cout<<"f+ rat"<<irat+1<<": "<<smart_print(fp[irat])<<endl;
      //cout<<"f-: "<<smart_print(fm)<<endl;
      
    }
  
#endif
  
  return 0;
}
