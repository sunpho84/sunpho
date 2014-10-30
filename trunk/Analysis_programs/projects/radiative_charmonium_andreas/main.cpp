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

int ib;
double a[2]={1/1.73,1/2.3};

double *c_two_pts_ASL_fit[3],*e_two_pts_ASL_fit[3];
int TH_two_pts_ASL_fit;
int tmin_two_pts_ASL_fit;
int tmax_two_pts_ASL_fit;

double fun_two_pts_ASL_fit(double Z1,double Z2,double M,double t,double(*f)(double))
{return Z1*Z2*exp(-M*TH_two_pts_ASL_fit)*f(M*(TH_two_pts_ASL_fit-t))/M;}

void ch2_two_pts_ASL_fit(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  double M=p[0];
  double ZA=p[1];
  double ZS=p[2];
  double ZL=p[3];

  for(int t=tmin_two_pts_ASL_fit;t<min(tmax_two_pts_ASL_fit,TH_two_pts_ASL_fit);t++)
    {
      double teo_ASL[3];
      teo_ASL[0]=fun_two_pts_ASL_fit(ZA,ZS,M,t,sinh);
      teo_ASL[1]=fun_two_pts_ASL_fit(ZS,ZS,M,t,cosh);
      teo_ASL[2]=fun_two_pts_ASL_fit(ZL,ZS,M,t,cosh);
      
      for(int isl=0;isl<3;isl++)
	{
	  double num=c_two_pts_ASL_fit[isl][t];
	  double teo=teo_ASL[isl];
	  double diff=num-teo;
	  double err=e_two_pts_ASL_fit[isl][t];
	  double cont=sqr(diff/err);
	  ch+=cont;
	  //cout<<t<<" "<<isl<<" (("<<teo<<"-"<<num<<")"<<"/"<<err<<")^2="<<cont<<endl;
	}
    }
}

void two_pts_ASL_fit(jack &M,jack &ZA,jack &ZS,jack &ZL,jvec corrSA,jvec corrSS,jvec corrSL,int tmin,int tmax,const char *path)
{
  int TH=TH_two_pts_ASL_fit=corrSL.nel-1;
  
  //perform the fit and set an initial estimate
  jack MSA,MSS,MSL,ZSA,ZSS,ZSL;
  two_pts_fit(MSA,ZSA,corrSA,tmin,tmax,NULL,NULL,TH,-1);
  two_pts_fit(MSS,ZSS,corrSS,tmin,tmax);
  two_pts_fit(MSL,ZSL,corrSL,tmin,tmax);
  
  //get estimates for ZS and ZL
  M=MSL;
  ZS=sqrt(ZSS);
  ZL=ZSL/ZS;
  ZA=ZSA/ZS;
  
  //define minimzer
  TMinuit minu;
  minu.SetPrintLevel(-1);
  minu.SetFCN(ch2_two_pts_ASL_fit);
  
  //copy parameters
  int njack=M.njack;
  tmin_two_pts_ASL_fit=tmin;
  tmax_two_pts_ASL_fit=tmax;
  
  //define temporary structures
  c_two_pts_ASL_fit[0]=new double[TH+1];
  c_two_pts_ASL_fit[1]=new double[TH+1];
  c_two_pts_ASL_fit[2]=new double[TH+1];
  e_two_pts_ASL_fit[0]=new double[TH+1];
  e_two_pts_ASL_fit[1]=new double[TH+1];
  e_two_pts_ASL_fit[2]=new double[TH+1];
  
  //copy errors
  for(int iel=0;iel<=TH;iel++)
    {
      e_two_pts_ASL_fit[0][iel]=corrSA[iel].err();
      e_two_pts_ASL_fit[1][iel]=corrSS[iel].err();
      e_two_pts_ASL_fit[2][iel]=corrSL[iel].err();
    }
  
  //loop over jacknife
  for(int ijack_fit=0;ijack_fit<=njack;ijack_fit++)
    {
      //set pars
      minu.DefineParameter(0,"M",M[ijack_fit],M.err(),0,0);
      minu.DefineParameter(1,"ZA",ZA[ijack_fit],ZA.err(),0,0);
      minu.DefineParameter(2,"ZS",ZS[ijack_fit],ZS.err(),0,0);
      minu.DefineParameter(3,"ZL",ZL[ijack_fit],ZL.err(),0,0);
      
      //copy elements
      for(int iel=0;iel<=TH;iel++)
        {
          c_two_pts_ASL_fit[0][iel]=corrSA[iel][ijack_fit];
          c_two_pts_ASL_fit[1][iel]=corrSS[iel][ijack_fit];
          c_two_pts_ASL_fit[2][iel]=corrSL[iel][ijack_fit];
        }
      
      //minimize
      minu.Migrad();
      
      //get back parameters
      double dum;
      minu.GetParameter(0,M[ijack_fit],dum);
      minu.GetParameter(1,ZA[ijack_fit],dum);
      minu.GetParameter(2,ZS[ijack_fit],dum);
      minu.GetParameter(3,ZL[ijack_fit],dum);
    }
  
  //write plots
  write_constant_fit_plot(path,effective_mass(corrSA,TH,-1),M,tmin,tmax);
  append_constant_fit_plot(path,effective_mass(corrSS),M,tmin,tmax,3);
  append_constant_fit_plot(path,effective_mass(corrSL),M,tmin,tmax,6);
  
  //delete
  delete[] c_two_pts_ASL_fit[0];
  delete[] c_two_pts_ASL_fit[1];
  delete[] c_two_pts_ASL_fit[2];
  delete[] e_two_pts_ASL_fit[0];
  delete[] e_two_pts_ASL_fit[1];
  delete[] e_two_pts_ASL_fit[2];
}

void read_analysis_pars()
{
  debug_load=0;
  debug_fit=0;
  
  FILE *input_file=open_file("analysis_pars","r");
  read_formatted_from_file_expecting((char*)(&ib),input_file,"%d","IB");
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
  
  int ncorr=23;
  jvec Pth0_sm_lo=jvec_load(path,T,njacks,ncorr*im+0).simmetrized(1);
  jvec Pth0_sm_sm=jvec_load(path,T,njacks,ncorr*im+1).simmetrized(1);
  jvec Ath0_sm_lo=jvec_load(path,T,njacks,ncorr*im+2).simmetrized(-1);
  jvec Pth1_sm_lo=jvec_load(path,T,njacks,ncorr*im+3).simmetrized(1);
  jvec Pth1_sm_sm=jvec_load(path,T,njacks,ncorr*im+4).simmetrized(1);
  jvec Pth2_sm_lo=jvec_load(path,T,njacks,ncorr*im+5).simmetrized(1);
  jvec Pth2_sm_sm=jvec_load(path,T,njacks,ncorr*im+6).simmetrized(1);
  jvec Vth0_sm_lo=jvec_load(path,T,njacks,ncorr*im+7).simmetrized(1);
  jvec Vth0_sm_sm=jvec_load(path,T,njacks,ncorr*im+8).simmetrized(1);
  jvec Vth1_sm_lo=jvec_load(path,T,njacks,ncorr*im+9).simmetrized(1);
  jvec Vth1_sm_sm=jvec_load(path,T,njacks,ncorr*im+10).simmetrized(1);
  jvec Vth2_sm_lo=jvec_load(path,T,njacks,ncorr*im+11).simmetrized(1);
  jvec Vth2_sm_sm=jvec_load(path,T,njacks,ncorr*im+12).simmetrized(1);
  jvec Bth0_sm_lo=jvec_load(path,T,njacks,ncorr*im+13).simmetrized(1);
  jvec Bth0_sm_sm=jvec_load(path,T,njacks,ncorr*im+14).simmetrized(1);
  jvec Bth1_sm_lo=jvec_load(path,T,njacks,ncorr*im+15).simmetrized(1);
  jvec Bth1_sm_sm=jvec_load(path,T,njacks,ncorr*im+16).simmetrized(1);
  jvec Bth2_sm_lo=jvec_load(path,T,njacks,ncorr*im+17).simmetrized(1);
  jvec Bth2_sm_sm=jvec_load(path,T,njacks,ncorr*im+18).simmetrized(1);
  
  jvec P5VJVK=jvec_load(path,T,njacks,ncorr*im+19);
  jvec P5V0P5_A=jvec_load(path,T,njacks,ncorr*im+20);
  jvec P5V0P5_B=jvec_load(path,T,njacks,ncorr*im+21);
  jvec P5VJBK=jvec_load(path,T,njacks,ncorr*im+22);
  
  if(tsep!=T/2)
    {
      P5VJVK=P5VJVK.subset(0,tsep+1);
      P5V0P5_A=P5V0P5_A.subset(0,tsep+1);
      P5V0P5_B=P5V0P5_B.subset(0,tsep+1);
      P5VJBK=P5VJBK.subset(0,tsep+1);
    }
  else
    {
      P5VJVK=P5VJVK.simmetrized(-1);
      P5V0P5_A=P5V0P5_A.simmetrized(-1);
      P5V0P5_B=P5V0P5_B.simmetrized(-1);
      P5VJBK=P5VJBK.simmetrized(-1);
    }  
  
  //momentum as expected from theta
  double inputP1=th[im][0]*M_PI/L;
  double inputP2=th[im][1]*M_PI/L;
  
  jack dum(njacks);
  jack EP0(njacks),ZPS0(njacks),ZPL0(njacks);
  jack EA0(njacks),ZAS0(njacks),ZAL0(njacks);
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
  two_pts_ASL_fit(EA0,ZAL0,ZPS0,ZPL0,Ath0_sm_lo,Pth0_sm_sm,Pth0_sm_lo,tminPl,tmaxPl,
		  combine("plots/%d/ASL_th0.xmg",im).c_str());
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
  
  //check dispertion relation/speed of light
  cout<<"-speed of light extracted from pseudoscalar energy:"<<endl;
  cout<<" th1, Fitted energy: "<<smart_print(EP1)<<", Latt.disp.rel: "<<smart_print(latt_e(EP0,inputP1))
      <<", cont: "<<smart_print(cont_e(EP0,inputP1))<<", exp c: "<<smart_print(sqrt(EP1*EP1-EP0*EP0)/inputP1/sqrt(3))<<endl;
  cout<<" th2, Fitted energy: "<<smart_print(EP2)<<", Latt.disp.rel: "<<smart_print(latt_e(EP0,inputP2))
      <<", cont: "<<smart_print(cont_e(EP0,inputP2))<<", exp c: "<<smart_print(sqrt(EP2*EP2-EP0*EP0)/inputP2/sqrt(3))<<endl;
  
  //find target theta
  jack targ_p1V_2rf=find_p1_2rf(EP0,EV0);
  jack targ_p1B_2rf=find_p1_2rf(EP0,EB0);
  jack targ_p1V_2rf_cont=find_p1_2rf(EP0,EV0,cont_e);
  jack targ_p1B_2rf_cont=find_p1_2rf(EP0,EB0,cont_e);
  cout<<"-target theta:"<<endl;
  cout<<" in J/Psi rest frame: "<<smart_print(targ_p1V_2rf*L/M_PI)<<", assuming continuum disp rel: "<<smart_print(targ_p1V_2rf_cont*L/M_PI)<<endl;
  cout<<" in B rest frame:     "<<smart_print(targ_p1B_2rf*L/M_PI)<<", assuming continuum disp rel: "<<smart_print(targ_p1B_2rf_cont*L/M_PI)<<endl;
  
  //extract momentum from P and V
  jack exP1=sqrt((EP1*EP1-EP0*EP0)/3);
  jack exP2=sqrt((EP2*EP2-EP0*EP0)/3);
  jack exV1=sqrt((EV1*EV1-EV0*EV0)/3);
  jack exV2=sqrt((EV2*EV2-EV0*EV0)/3);
  jack exB1=sqrt((EB1*EB1-EB0*EB0)/3);
  jack exB2=sqrt((EB2*EB2-EB0*EB0)/3);
  
  cout<<"-extr_mom:"<<endl;
  cout<<" mom1: P="<<smart_print(exP1)<<"=ex*"<<smart_print(exP1/inputP1)<<", V="<<smart_print(exV1)<<", B="<<smart_print(exB1)<<", input: "<<inputP1<<endl;
  cout<<" mom2: P="<<smart_print(exP2)<<"=ex*"<<smart_print(exP2/inputP2)<<", V="<<smart_print(exV2)<<", B="<<smart_print(exB2)<<", input: "<<inputP2<<endl;
  cout<<"-targ_mom:"<<endl;
  cout<<" in J/Psi rest frame: "<<smart_print(targ_p1V_2rf)<<endl;
  cout<<" in B     rest frame: "<<smart_print(targ_p1B_2rf)<<endl;
  cout<<"-Q2:"<<endl;
  cout<<" in J/Psi rest frame: "<<smart_print((sqr(EV0-EP1)-3*sqr(inputP1)))<<endl;
  cout<<" in B rest frame:     "<<smart_print((sqr(EB0-EP2)-3*sqr(inputP2)))<<endl;
  
  //define dT for V and B
  jvec dTV_nu(tsep+1,njacks),dTV_sa(tsep+1,njacks);
  jvec dTB_nu(tsep+1,njacks),dTB_sa(tsep+1,njacks);
  for(int t=0;t<=tsep;t++)
    {
      dTV_nu[t]=Vth0_sm_lo[t]*Pth1_sm_lo[tsep-t]/ZPL1/ZVL0;
      dTV_sa[t]=ZPS1*ZVS0*exp(-EV0*t-EP1*(tsep-t))/(2*EP1*2*EV0);
      dTB_nu[t]=Bth0_sm_lo[t]*Pth2_sm_lo[tsep-t]/ZPL2/ZBL0;
      dTB_sa[t]=ZPS2*ZBS0*exp(-EB0*t-EP2*(tsep-t))/(2*EP2*2*EB0);
    }
  dTV_nu.print_to_file(combine("plots/%d/dTV_nu.xmg",im).c_str());
  dTV_sa.print_to_file(combine("plots/%d/dTV_sa.xmg",im).c_str());
  dTB_nu.print_to_file(combine("plots/%d/dTB_nu.xmg",im).c_str());
  dTB_sa.print_to_file(combine("plots/%d/dTB_sa.xmg",im).c_str());
  
  //Zv
  P5V0P5_A.print_to_file(combine("plots/%d/P5V0P5_A.xmg",im).c_str());
  P5V0P5_B.print_to_file(combine("plots/%d/P5V0P5_B.xmg",im).c_str());
  jvec Zv_corr_A=1/P5V0P5_A,Zv_corr_B=1/P5V0P5_B;
  for(int t=0;t<tsep;t++)
    {
      Zv_corr_A[t]*=EP1*ZPS0*ZPS1*exp(-EP0*t-EP1*(tsep-t))/(2*EP1*EP0);
      Zv_corr_B[t]*=EP2*ZPS0*ZPS2*exp(-EP0*t-EP2*(tsep-t))/(2*EP2*EP0);
    }
  Zv_corr_A.print_to_file(combine("plots/%d/Zv_A.xmg",im).c_str());
  Zv_corr_B.print_to_file(combine("plots/%d/Zv_B.xmg",im).c_str());
  jack Zv_A=-Zv_corr_A[tsep/2];
  jack Zv_B=-Zv_corr_B[tsep/2];
  jack Zv_slope=(Zv_B-Zv_A)/(inputP2-inputP1);
  jack Zv=Zv_A-Zv_slope*inputP1;
  cout<<"-Zv:"<<endl;
  cout<<" th0: "<<smart_print(Zv)<<" (extr)"<<endl;
  cout<<" th1: "<<smart_print(Zv_A)<<endl;
  cout<<" th2: "<<smart_print(Zv_B)<<endl;
  
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
  
  //compute F(J/Psi)
  jack afjpsi=Zv*ZVL0/EV0;
  cout<<"-f(J/Psi)"<<endl;
  cout<<" af: "<<smart_print(afjpsi)<<endl;
  cout<<" f:  "<<smart_print(afjpsi/a[ib])<<endl;
  
  //compute F(eta)
  jack afeta=-Zv*ZAL0/EP0;
  cout<<"-f(eta) (WARNING - remember to fix Zv->Za)"<<endl;
  cout<<" af: "<<smart_print(afeta)<<endl;
  cout<<" f:  "<<smart_print(afeta/a[ib])<<endl;
  
  //compute masses
  cout<<"-Masses:"<<endl;
  cout<<" M(eta):    "<<smart_print(EP0/a[ib])<<" GeV [experimental: 2.984 GeV]"<<endl;
  cout<<" M(J/Psi):  "<<smart_print(EV0/a[ib])<<" GeV [diff with eta: "<<smart_print((EV0-EP0)/a[ib]*1000)<<" MeV]"<<endl;
  cout<<" M(h_c):    "<<smart_print(EB0/a[ib])<<" GeV [diff with eta: "<<smart_print((EB0-EP0)/a[ib]*1000)<<" MeV]"<<endl;
}

int main(int narg,char **arg)
{
  read_analysis_pars();
  
  cout<<"Beta: "<<ib<<endl;
  
  for(int im=0;im<3;im++)
    {
      cout<<"====================================="<<endl;
      analysis(im);
    }
  
  return 0;
}
