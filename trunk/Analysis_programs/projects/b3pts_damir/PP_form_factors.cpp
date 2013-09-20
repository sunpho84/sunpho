#include "common.cpp"

//#define USE_PROJ

int im_spec,im_S0,im_S1;

void read_input_header(const char *path)
{
  FILE *input_file=open_file(path,"r");

  read_formatted_from_file_expecting((char*)&im_spec,input_file,"%d","im_spec");
  read_formatted_from_file_expecting((char*)&im_S0,input_file,"%d","im_S0");
  read_formatted_from_file_expecting((char*)&im_S1,input_file,"%d","im_S1");

  fclose(input_file);
}

void solve_fP_f0(jack &fP,jack &fM,jack &f0,jack &V0,jack &VK,jack &P2,jack &P0,double PK,jack &Q2,jack &Q0,double QK,jack &M_K,jack &M_H)
{
  //a global factor 3 is dropped
  jack  delta=P0*QK-Q0*PK;
  jack deltaP=V0*QK-Q0*VK;
  jack deltaM=P0*VK-V0*PK;
  
  fP=deltaP/delta;
  fM=deltaM/delta;
  
  f0=fP+fM*Q2/(M_H*M_H-M_K*M_K);
}

int main()
{
  read_set_pars("../data_pars");
  read_input_header("analysis_pars");
  int nth_S0=5;
  //load the masses
  jack M_K(njack),M_H(njack);
  M_K.load("masses_Z/M_K_P5",0);
  M_H.load("masses_Z/M_H_P5",0);
  
  //load the matrix elements
  jvec V0(nth_S0,njack),VK(nth_S0,njack),TK(nth_S0,njack),S0(nth_S0,njack);
  V0.load("matrix_elements/P5_V0_P5",0);
  VK.load("matrix_elements/P5_VK_P5",0);
  TK.load("matrix_elements/P5_TK_P5",0);
  S0.load("matrix_elements/P5_S0_P5",0);
  
  //loads divita ratios
  jvec rat1(nth_S0,njack),rat2(nth_S0,njack);
  rat1.load("matrix_elements/divita_ratio_1",0);
  rat2.load("matrix_elements/divita_ratio_2",0);
  
  //include renormalization
  TK*=Zt[ibeta];
  VK*=Zv[ibeta];
  V0*=Zv[ibeta];
  //S0/=Zp_fr_Zs[ibeta]; //not required: OS case
  
  //include mass diff to achieve S0 quadri-div
  S0*=mass[nm_S0-1+im_S1]-mass[im_S0];
  cout<<"masses: "<<mass[nm_S0-1+im_S1]<<" - "<<mass[im_S0]<<endl;
  
  //open outputs
  ofstream fT_out("form_factors/plots/fT_PP.xmg");
  ofstream fP_out("form_factors/plots/fP_PP.xmg");
  ofstream fM_out("form_factors/plots/fM_PP.xmg");
  ofstream fPd_out("form_factors/plots/fPd_PP.xmg");
  ofstream f0d_out("form_factors/plots/f0d_PP.xmg");
  ofstream f0_out("form_factors/plots/f0_PP.xmg");
  ofstream f0s_out("form_factors/plots/f0s_PP.xmg");
  ofstream fT_fr_fP_out("form_factors/plots/fT_fr_fP_PP.xmg");
  ofstream f0_fr_fP_out("form_factors/plots/f0_fr_fP_PP.xmg");
  fT_out<<"@type xydxdy"<<endl;
  fP_out<<"@type xydxdy"<<endl;
  fM_out<<"@type xydxdy"<<endl;
  fPd_out<<"@type xydxdy"<<endl;
  f0d_out<<"@type xydxdy"<<endl;
  f0_out<<"@type xydxdy"<<endl;
  f0s_out<<"@type xydxdy"<<endl;
  fT_fr_fP_out<<"@type xydxdy"<<endl;
  f0_fr_fP_out<<"@type xydxdy"<<endl;

  jvec Q2(nth_S0,njack),EP(nth_S0,njack);
  jack P2,P0,Q0;
  jvec f0(nth_S0,njack),fP(nth_S0,njack),fT(nth_S0,njack),fM(nth_S0,njack);
  jvec f0s(nth_S0,njack),fPd(nth_S0,njack),f0d(nth_S0,njack);
  for(int ith=0;ith<nth_S0;ith++)
    {
      //compute cinematic factors
      double PK,QK;
      compute_momentums(P2,P0,PK, Q2[ith],Q0,QK, M_K,th_S0[ith], M_H);
      EP[ith]=cont_en(M_K,th_S0[ith]);
      
      //compute f0 using V0 and VK
#ifdef USE_PROJ
      f0[ith]=(Q0*V0[ith]-3*QK*VK[ith])/(M_H*M_H-M_K*M_K);
#else
      if(ith!=0) solve_fP_f0(fP[ith],fM[ith],f0[ith],V0[ith],VK[ith],P2,P0,PK,Q2[ith],Q0,QK,M_K,M_H);
      else
	{
	  fP[ith]=0;
	  f0[ith]=V0[ith]/(M_K+M_H);
	}
#endif
      f0_out<<Q2[ith].med()<<" "<<f0[ith].med()<<" "<<Q2[ith].err()<<" "<<f0[ith].err()<<endl;
      
      //for test of f0(Q2_max,mi+1)/f0(Q2_max,mi)
      if(ith==0) f0[ith].write_to_binfile("form_factors/f0_Q2_max");
      
      //compute f0 using S0
      f0s[ith]=-S0[ith]/(M_H*M_H-M_K*M_K);
      f0s_out<<Q2[ith].med()<<" "<<f0s[ith].med()<<" "<<Q2[ith].err()<<" "<<f0s[ith].err()<<endl;
      
      //fT and fP defined only at non null momentum
      if(ith!=0)
	{
	  fT[ith]=TK[ith]*(M_K+M_H)/(2*M_H*PK);
	  //factor 3 on K simplified
#ifdef USE_PROJ
	  fP[ith]=(V0[ith]-VK[ith]*(M_H-cont_en(M_K,th_S0[ith]))/QK)/(2*M_H);
#endif
	  fT_out<<Q2[ith].med()<<" "<<fT[ith].med()<<" "<<Q2[ith].err()<<" "<<fT[ith].err()<<endl;
	  fP_out<<Q2[ith].med()<<" "<<fP[ith].med()<<" "<<Q2[ith].err()<<" "<<fP[ith].err()<<endl;
	  fM_out<<Q2[ith].med()<<" "<<fM[ith].med()<<" "<<Q2[ith].err()<<" "<<fM[ith].err()<<endl;
	  
	  f0_fr_fP_out<<Q2[ith].med()<<" "<<(f0[ith]/fP[ith]).med()<<" "<<Q2[ith].err()<<" "<<(f0[ith]/fP[ith]).err()<<endl;
	  fT_fr_fP_out<<Q2[ith].med()<<" "<<(fT[ith]/fP[ith]).med()<<" "<<Q2[ith].err()<<" "<<(fT[ith]/fP[ith]).err()<<endl;
      
	  //////////////////// use divita ratios ////////////
	  
	  double pp=momentum(th_S0[ith]);
	  jack xi_reco=fM[ith]/fP[ith];
	  //cout<<" check rat2: "<<rat2[ith]<<" "<<(1-xi_reco)*pp/(M_H+EP[ith]+(M_H-EP[ith])*xi_reco)<<endl;
	  
	  jack xi=(rat2[ith]*(EP[ith]+M_H)-pp)/(rat2[ith]*(EP[ith]-M_H)-pp);
	  cout<<" xi check: "<<xi<<" "<<xi_reco<<endl;
	  fPd[ith]=f0[0]*rat1[ith]*(M_H+M_K)/(M_H+EP[ith]+(M_H-EP[ith])*xi);
	  f0d[ith]=fPd[ith]*(1+xi*Q2[ith]/(M_H*M_H-M_K*M_K));
	  
	  fPd_out<<Q2[ith].med()<<" "<<fPd[ith].med()<<" "<<Q2[ith].err()<<" "<<fPd[ith].err()<<endl;
	  f0d_out<<Q2[ith].med()<<" "<<f0d[ith].med()<<" "<<Q2[ith].err()<<" "<<f0d[ith].err()<<endl;
	}
      
      f0d[0]=f0[0];
    }
  
  char results_path[]="form_factors/EP_Q2_fP_f0_fT_f0s";
  EP.write_to_binfile(results_path);
  Q2.append_to_binfile(results_path);
  fP.append_to_binfile(results_path);
  f0.append_to_binfile(results_path);
  fT.append_to_binfile(results_path);
  f0s.append_to_binfile(results_path);
  
  return 0;
}
