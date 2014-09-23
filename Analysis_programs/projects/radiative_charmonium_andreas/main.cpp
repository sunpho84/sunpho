#include "include.h"

#include "find_theta_internal.cpp"

int njacks;
int T,L,tsep;
int tminPl,tmaxPl;
int tminPs,tmaxPs;
int tminVl,tmaxVl;
int tminVs,tmaxVs;
int tminBl,tmaxBl;
int tminBs,tmaxBs;
double th[3][2];

double a=1/1.73;

void read_analysis_pars()
{
  debug_load=0;
  debug_fit=0;
  
  FILE *input_file=open_file("analysis_pars","r");
  read_formatted_from_file_expecting((char*)(&T),input_file,"%d","T");
  read_formatted_from_file_expecting((char*)(&L),input_file,"%d","L");
  read_formatted_from_file_expecting((char*)(&tsep),input_file,"%d","tsep");
  read_formatted_from_file_expecting((char*)(&njacks),input_file,"%d","njacks");
  read_formatted_from_file_expecting((char*)(&tminPl),input_file,"%d","tintPl");
  read_formatted_from_file((char*)(&tmaxPl),input_file,"%d","tintPl");
  read_formatted_from_file_expecting((char*)(&tminPs),input_file,"%d","tintPs");
  read_formatted_from_file((char*)(&tmaxPs),input_file,"%d","tintPs");
  read_formatted_from_file_expecting((char*)(&tminVl),input_file,"%d","tintVl");
  read_formatted_from_file((char*)(&tmaxVl),input_file,"%d","tintVl");
  read_formatted_from_file_expecting((char*)(&tminVs),input_file,"%d","tintVs");
  read_formatted_from_file((char*)(&tmaxVs),input_file,"%d","tintVs");
  read_formatted_from_file_expecting((char*)(&tminBl),input_file,"%d","tintBl");
  read_formatted_from_file((char*)(&tmaxBl),input_file,"%d","tintBl");
  read_formatted_from_file_expecting((char*)(&tminBs),input_file,"%d","tintBs");
  read_formatted_from_file((char*)(&tmaxBs),input_file,"%d","tintBs");
  expect_string_from_file(input_file,"thetas");
  for(int im=0;im<3;im++)
    for(int ith=0;ith<2;ith++)
      read_formatted_from_file((char*)(&th[im][ith]),input_file,"%lg","tintBs");
  
  fclose(input_file);
}
 
void analysis(int im)
{
  const char path[]="out.dat";
  
  int ncorr=21;
  jvec Pth0_sm_lo=jvec_load(path,T,njacks,ncorr*im+0).simmetrized(1);
  jvec Pth0_sm_sm=jvec_load(path,T,njacks,ncorr*im+1).simmetrized(1);
  jvec Pth1_sm_lo=jvec_load(path,T,njacks,ncorr*im+2).simmetrized(1);
  jvec Pth1_sm_sm=jvec_load(path,T,njacks,ncorr*im+3).simmetrized(1);
  jvec Pth2_sm_lo=jvec_load(path,T,njacks,ncorr*im+4).simmetrized(1);
  jvec Pth2_sm_sm=jvec_load(path,T,njacks,ncorr*im+5).simmetrized(1);
  jvec Vth0_sm_lo=jvec_load(path,T,njacks,ncorr*im+6).simmetrized(1);
  jvec Vth0_sm_sm=jvec_load(path,T,njacks,ncorr*im+7).simmetrized(1);
  jvec Vth1_sm_lo=jvec_load(path,T,njacks,ncorr*im+8).simmetrized(1);
  jvec Vth1_sm_sm=jvec_load(path,T,njacks,ncorr*im+9).simmetrized(1);
  jvec Vth2_sm_lo=jvec_load(path,T,njacks,ncorr*im+10).simmetrized(1);
  jvec Vth2_sm_sm=jvec_load(path,T,njacks,ncorr*im+11).simmetrized(1);
  jvec Bth0_sm_lo=jvec_load(path,T,njacks,ncorr*im+12).simmetrized(1);
  jvec Bth0_sm_sm=jvec_load(path,T,njacks,ncorr*im+13).simmetrized(1);
  jvec Bth1_sm_lo=jvec_load(path,T,njacks,ncorr*im+14).simmetrized(1);
  jvec Bth1_sm_sm=jvec_load(path,T,njacks,ncorr*im+15).simmetrized(1);
  jvec Bth2_sm_lo=jvec_load(path,T,njacks,ncorr*im+16).simmetrized(1);
  jvec Bth2_sm_sm=jvec_load(path,T,njacks,ncorr*im+17).simmetrized(1);
  
  jvec P5VJVK=jvec_load(path,T,njacks,ncorr*im+18).subset(0,tsep);
  jvec P5V0P5=jvec_load(path,T,njacks,ncorr*im+19).subset(0,tsep);
  jvec P5VJBK=jvec_load(path,T,njacks,ncorr*im+20).subset(0,tsep);
  
  jack dum(njacks);
  jack EP0(njacks),ZPS0(njacks),ZPL0(njacks);
  jack EP1(njacks),ZPS1(njacks),ZPL1(njacks);
  jack EP2(njacks),ZPS2(njacks),ZPL2(njacks);
  jack EV0(njacks),ZVS0(njacks),ZVL0(njacks);
  jack EV1(njacks),ZVS1(njacks),ZVL1(njacks);
  jack EV2(njacks),ZVS2(njacks),ZVL2(njacks);
  jack EB0(njacks),ZBS0(njacks),ZBL0(njacks);
  jack EB1(njacks),ZBS1(njacks),ZBL1(njacks);
  jack EB2(njacks),ZBS2(njacks),ZBL2(njacks);
  two_pts_SL_fit(EP0,ZPL0,ZPS0,Pth0_sm_lo,Pth0_sm_sm,tminPl,tmaxPl,tminPs,tmaxPs,
		 combine("plots/%d/PSL_th0.xmg",im).c_str());
  two_pts_SL_fit(EP1,ZPL1,ZPS1,Pth1_sm_lo,Pth1_sm_sm,tminPl,tmaxPl,tminPs,tmaxPs,
		 combine("plots/%d/PSL_th1.xmg",im).c_str());
  two_pts_SL_fit(EP2,ZPL2,ZPS2,Pth2_sm_lo,Pth2_sm_sm,tminPl,tmaxPl,tminPs,tmaxPs,
		 combine("plots/%d/PSL_th2.xmg",im).c_str());
  two_pts_SL_fit(EV0,ZVL0,ZVS0,Vth0_sm_lo,Vth0_sm_sm,tminVl,tmaxVl,tminVs,tmaxVs,
		 combine("plots/%d/VSL_th0.xmg",im).c_str());
  two_pts_SL_fit(EV1,ZVL1,ZVS1,Vth1_sm_lo,Vth1_sm_sm,tminVl,tmaxVl,tminVs,tmaxVs,
		 combine("plots/%d/VSL_th1.xmg",im).c_str());
  two_pts_SL_fit(EV2,ZVL2,ZVS2,Vth2_sm_lo,Vth2_sm_sm,tminVl,tmaxVl,tminVs,tmaxVs,
		 combine("plots/%d/VSL_th2.xmg",im).c_str());
  two_pts_SL_fit(EB0,ZBL0,ZBS0,Bth0_sm_lo,Bth0_sm_sm,tminBl,tmaxBl,tminBs,tmaxBs,
		 combine("plots/%d/BSL_th0.xmg",im).c_str());
  two_pts_SL_fit(EB1,ZBL1,ZBS1,Bth1_sm_lo,Bth1_sm_sm,tminBl,tmaxBl,tminBs,tmaxBs,
		 combine("plots/%d/BSL_th1.xmg",im).c_str());
  two_pts_SL_fit(EB2,ZBL2,ZBS2,Bth2_sm_lo,Bth2_sm_sm,tminBl,tmaxBl,tminBs,tmaxBs,
		 combine("plots/%d/BSL_th2.xmg",im).c_str());
  cout<<"EP0: "<<smart_print(EP0)<<", EV0: "<<smart_print(EV0)<<", EB0: "<<smart_print(EB0)<<endl;
  cout<<"EP1: "<<smart_print(EP1)<<", EV1: "<<smart_print(EV1)<<", EB1: "<<smart_print(EB1)<<endl;
  cout<<"EP2: "<<smart_print(EP2)<<", EV2: "<<smart_print(EV2)<<", EB2: "<<smart_print(EB2)<<endl;
  
  //find target theta
  jack targ_p1V_2rf=find_p1_2rf(EP0,EV0);
  jack targ_p1B_2rf=find_p1_2rf(EP0,EB0);  
  cout<<"-target theta:"<<endl;
  cout<<" in J/Psi rest frame: "<<smart_print(targ_p1V_2rf*L/M_PI)<<endl;
  cout<<" in B rest frame:     "<<smart_print(targ_p1B_2rf*L/M_PI)<<endl;
  
  //extract momentum from P and V
  jack exP1=sqrt((EP1*EP1-EP0*EP0)/3);
  jack exP2=sqrt((EP2*EP2-EP0*EP0)/3);
  jack exV1=sqrt((EV1*EV1-EV0*EV0)/3);
  jack exV2=sqrt((EV2*EV2-EV0*EV0)/3);
  jack exB1=sqrt((EB1*EB1-EB0*EB0)/3);
  jack exB2=sqrt((EB2*EB2-EB0*EB0)/3);
  
  //momentum as expected from theta
  double thP1=th[im][0]*M_PI/L;
  double thP2=th[im][1]*M_PI/L;
  
  cout<<"-extr_mom:"<<endl;
  cout<<" mom1: P="<<smart_print(exP1)<<" V="<<smart_print(exV1)<<", B="<<smart_print(exB1)<<endl;
  cout<<" mom2: P="<<smart_print(exP2)<<" V="<<smart_print(exV2)<<", B="<<smart_print(exB2)<<endl;
  cout<<"-targ_mom:"<<endl;
  cout<<" in J/Psi rest frame: "<<smart_print(targ_p1V_2rf)<<endl;
  cout<<" in B     rest frame: "<<smart_print(targ_p1B_2rf)<<endl;
  cout<<"-Q2:"<<endl;
  cout<<" in J/Psi rest frame: "<<smart_print((sqr(EV0-EP1)-3*sqr(thP1)))<<endl;
  cout<<" in B rest frame:     "<<smart_print((sqr(EB0-EP2)-3*sqr(thP2)))<<endl;
  
  //define dT for V and B
  jvec dTV_nu(tsep,njacks),dTV_sa(tsep,njacks);
  jvec dTB_nu(tsep,njacks),dTB_sa(tsep,njacks);
  for(int t=0;t<tsep;t++)
    {
      dTV_nu[t]=Vth0_sm_lo[t]*Pth1_sm_lo[tsep-t]/ZPL1/ZVL0;
      dTV_sa[t]=ZPS1*ZVS0*exp(-EV0*t-EP1*(tsep-t))/(2*EP1*2*EV0);
      dTB_nu[t]=Bth0_sm_lo[t]*Pth1_sm_lo[tsep-t]/ZPL1/ZBL0;
      dTB_sa[t]=ZPS1*ZBS0*exp(-EB0*t-EP2*(tsep-t))/(2*EP2*2*EB0);
    }
  dTV_nu.print_to_file(combine("plots/%d/dTV_nu.xmg",im).c_str());
  dTV_sa.print_to_file(combine("plots/%d/dTV_sa.xmg",im).c_str());
  dTB_nu.print_to_file(combine("plots/%d/dTB_nu.xmg",im).c_str());
  dTB_sa.print_to_file(combine("plots/%d/dTB_sa.xmg",im).c_str());
  
  //Zv
  P5V0P5.print_to_file(combine("plots/%d/P5V0P5.xmg",im).c_str());
  jvec Zv_corr=1/P5V0P5;
  for(int t=0;t<tsep;t++) Zv_corr[t]*=EP1*ZPS0*ZPS1*exp(-EP0*t-EP1*(tsep-t))/(2*EP1*EP0);
  Zv_corr.print_to_file(combine("plots/%d/Zv.xmg",im).c_str());
  jack Zv=-Zv_corr[tsep/2];
  cout<<"Zv: "<<smart_print(Zv)<<endl;
  
  //define ratios
  jvec RV_nu_corr=P5VJVK/dTV_nu;
  jvec RV_sa_corr=P5VJVK/dTV_sa;
  jvec RB_nu_corr=P5VJBK/dTB_nu;
  jvec RB_sa_corr=P5VJBK/dTB_sa;
  jack RV_nu=constant_fit(RV_nu_corr,tsep/2-tsep/4,tsep/2+tsep/4,combine("plots/%d/RV_nu.xmg",im).c_str());
  jack RV_sa=constant_fit(RV_sa_corr,tsep/2-tsep/4,tsep/2+tsep/4,combine("plots/%d/RV_sa.xmg",im).c_str());
  jack RB_nu=constant_fit(RB_nu_corr,tsep/2-tsep/4,tsep/2+tsep/4,combine("plots/%d/RB_nu.xmg",im).c_str());
  jack RB_sa=constant_fit(RB_sa_corr,tsep/2-tsep/4,tsep/2+tsep/4,combine("plots/%d/RB_sa.xmg",im).c_str());
  
  //compute V
  cout<<"-V form factor"<<endl;
  jack V_nu=RV_nu/exP1*(EV0+EP0)/(2*EV0)*Zv;
  jack V_sa=RV_sa/exP1*(EV0+EP0)/(2*EV0)*Zv;
  cout<<" V_nu: "<<smart_print(V_nu)<<endl;
  cout<<" V_sa: "<<smart_print(V_sa)<<endl;
  
  //compute F
  cout<<"-F form factor"<<endl;
  jack F_nu=RB_nu/EB0*Zv;
  jack F_sa=RB_sa/EB0*Zv;
  cout<<" F_nu: "<<smart_print(F_nu)<<endl;
  cout<<" F_sa: "<<smart_print(F_sa)<<endl;
}

int main(int narg,char **arg)
{
  read_analysis_pars();
  
  for(int im=0;im<3;im++)
    {
      cout<<"====================================="<<endl;
      analysis(im);
    }
  
  return 0;
}
