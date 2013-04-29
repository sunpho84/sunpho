#include "common.cpp"

#define USE_PROJ

int im_spec,im_S0,im_S1;

void read_input_header(const char *path)
{
  FILE *input_file=open_file(path,"r");

  read_formatted_from_file_expecting((char*)&im_spec,input_file,"%d","im_spec");
  read_formatted_from_file_expecting((char*)&im_S0,input_file,"%d","im_S0");
  read_formatted_from_file_expecting((char*)&im_S1,input_file,"%d","im_S1");

  fclose(input_file);
}

void solve_fP_f0(jack &fP,jack &f0,jack &V0,jack &VK,jack &P2,jack &P0,double PK,jack &Q2,jack &Q0,double QK,jack &M_K,jack &M_H)
{
  //a global factor 3 is dropped
  jack  delta=P0*QK-Q0*PK;
  jack deltaP=V0*QK-Q0*VK;
  jack deltaM=P0*VK-V0*PK;
  
  fP=deltaP/delta;
  jack fM=deltaM/delta;
  
  f0=fP+fM*Q2/(M_H*M_H-M_K*M_K);
}

int main()
{
  read_set_pars("../data_pars");
  read_input_header("analysis_pars");
  int nth_S0=5; //ATTENTIONNNNNNNNNNNNNNNNNNNNNNNNNN only fist 5 theta
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
  
  //include renormalization
  TK*=Zt[ibeta];
  VK*=Zv[ibeta];
  V0*=Zv[ibeta];
  S0/=Zp_fr_Zs[ibeta];
  
  //include mass diff to achieve S0 quadri-div
  S0*=mass[nm_S0+im_S1]-mass[im_S0];
  //cout<<"masses: "<<mass[nm_S0+im_S1]<<" - "<<mass[im_S0]<<endl;
  
  //open outputs
  ofstream fT_out("form_factors/plots/fT_PP.xmg");
  ofstream fP_out("form_factors/plots/fP_PP.xmg");
  ofstream f0_out("form_factors/plots/f0_PP.xmg");
  ofstream f0s_out("form_factors/plots/f0s_PP.xmg");
  ofstream fT_fr_fP_out("form_factors/plots/fT_fr_fP_PP.xmg");
  ofstream f0_fr_fP_out("form_factors/plots/f0_fr_fP_PP.xmg");
  fT_out<<"@type xydxdy"<<endl;
  fP_out<<"@type xydxdy"<<endl;
  f0_out<<"@type xydxdy"<<endl;
  f0s_out<<"@type xydxdy"<<endl;
  fT_fr_fP_out<<"@type xydxdy"<<endl;
  f0_fr_fP_out<<"@type xydxdy"<<endl;

  for(int ith=0;ith<nth_S0;ith++)
    {
      //compute cinematic factors
      jack P2,P0,Q2,Q0;
      double PK,QK;
      compute_momentums(P2,P0,PK, Q2,Q0,QK, M_K,th_S0[ith], M_H);
      
      //compute f0 using V0 and VK
#ifdef USE_PROJ
      jack f0=(Q0*V0[ith]-3*QK*VK[ith])/(M_H*M_H-M_K*M_K);
#else
      jack f0(njack),fP(njack);
      if(ith!=0) solve_fP_f0(fP,f0,V0[ith],VK[ith],P2,P0,PK,Q2,Q0,QK,M_K,M_H);
      else
	{
	  fP=0;
	  f0=V0[ith]/(M_K+M_H);
	}
#endif
      f0_out<<Q2.med()<<" "<<f0.med()<<" "<<Q2.err()<<" "<<f0.err()<<endl;
      
      //for test of f0(Q2_max,mi+1)/f0(Q2_max,mi)
      if(ith==0) f0.write_to_binfile("f0_Q2_max");
      
      //compute f0 using S0
      jack f0s=-S0[ith]/(M_H*M_H-M_K*M_K);
      f0s_out<<Q2.med()<<" "<<f0s.med()<<" "<<Q2.err()<<" "<<f0s.err()<<endl;
      
      //fT and fP defined only at non null momentum
      if(ith!=0)
	{
	  jack fT=TK[ith]*(M_K+M_H)/(2*M_H*PK);
	  //factor 3 on K simplified
#ifdef USE_PROJ
	  jack fP=(V0[ith]-VK[ith]*(M_H-cont_en(M_K,th_S0[ith]))/QK)/(2*M_H);
#endif
	  fT_out<<Q2.med()<<" "<<fT.med()<<" "<<Q2.err()<<" "<<fT.err()<<endl;
	  fP_out<<Q2.med()<<" "<<fP.med()<<" "<<Q2.err()<<" "<<fP.err()<<endl;
	  
	  f0_fr_fP_out<<Q2.med()<<" "<<(f0/fP).med()<<" "<<Q2.err()<<" "<<(f0/fP).err()<<endl;
	  fT_fr_fP_out<<Q2.med()<<" "<<(fT/fP).med()<<" "<<Q2.err()<<" "<<(fT/fP).err()<<endl;	  
	}
    }
  
  return 0;
}
