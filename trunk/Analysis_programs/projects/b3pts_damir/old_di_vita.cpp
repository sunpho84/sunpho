#include "include.h"
#include "../nf2/common_pars.cpp"

const int nth=8;
int nm=12;
int i_sea=0,il[4];
const int RE=0,IM=1;
const int EVN=1,ODD=-1;
int T,TH,L;
double th_P[nth];
double xs[3],xs_phys;
double xc[3],xc_phys;
int ibeta;
int tmin,tmax; 
int tmin_D,tmin_P;

double momentum(double th)
{return 2*M_PI/TH*th;}

jack latt_en(jack M,double th)
{return 2*asinh(sqrt(3*sqr(sin(momentum(th)/2))+sqr(sinh(M/2))));}
jack cont_en(jack M,double th)
{return sqrt(M*M+3*sqr(momentum(th)));}

int icombo_3pts(int im_S1,int im_S0,int ith,int ri)
{return ri+2*(im_S0+nm*(im_S1+nm*ith));}

int icombo_2pts(int im1,int im2,int ith)
{return 0+2*(im2+nm*(im1+nm*ith));}

jvec load_3pts(const char *ME,int im_S1,int im_S0,int ith,int ri,int pa,int si)
{
  jvec c1(T,njack),c2(T,njack);
  c1.load(combine("CORRELATORS/3pts_%s.dat",ME).c_str(),icombo_3pts(im_S1,im_S0,ith,ri));
  if(im_S1!=im_S0) c2.load(combine("CORRELATORS/3pts_%s.dat",ME).c_str(),icombo_3pts(im_S0,im_S1,ith,ri));
  else c2=c1;
  
  return (c1.simmetrized(pa)+si*c2.simmetrized(pa).simmetric())/2;
}

jvec load_2pts(const char *ME,int im1,int im2,int ith)
{
  jvec c(T,njack);
  c.load("CORRELATORS/2pts_P5P5.dat",icombo_2pts(im1,im2,ith));
  return c.simmetrized(1);
}

void extr(jvec &ext_EP,jvec &ext_ED,jvec &ext_Q2,jvec &ext_fP,jvec &ext_fM,jvec &ext_f0,jvec &ext_fT,int il_sea,int il,int ic)
{
  ////////////////////////////////////////// R0 //////////////////////////////////////  

  jvec R0_corr;
  jack R0(njack);
  
  //load standing
  jvec ll0_st=load_3pts("V0",il,il,0,RE,ODD,1);
  jvec lc0_st=load_3pts("V0",ic,il,0,RE,ODD,1);
  jvec cc0_st=load_3pts("V0",ic,ic,0,RE,ODD,1);
  
  //build R0
  R0_corr=lc0_st*lc0_st.simmetric()/(cc0_st*ll0_st);
  
  //fit and plot
  R0=constant_fit(R0_corr,TH-tmax,tmax,combine("plots/R0_il_%d_ic_%d.xmg",il,ic).c_str());
  
  //////////////////////////////////////////// R2 ////////////////////////////////////
  
  jvec R2_corr[nth];
  jvec RT_corr[nth];
  jvec R2(nth,njack);
  jvec RT(nth,njack);
  
  ofstream out_R2(combine("plots/R2_il_%d_ic_%d.xmg",il,ic).c_str());
  ofstream out_RT(combine("plots/RT_il_%d_ic_%d.xmg",il,ic).c_str());
  jvec lcK_th[nth],lc0_th[nth],lcT_th[nth];
  for(int ith=0;ith<nth;ith++)
    {
      //load corrs
      lcK_th[ith]=load_3pts("VK",ic,il,ith,IM,EVN,-1)/(6*th_P[ith]);
      lc0_th[ith]=load_3pts("V0",ic,il,ith,RE,ODD,1);
      lcT_th[ith]=load_3pts("VTK",ic,il,ith,IM,ODD,1)/(6*th_P[ith]);
      
      //build ratios
      R2_corr[ith]=lcK_th[ith]/lc0_th[ith];
      RT_corr[ith]=lcT_th[ith]/lcK_th[ith];
      
      //fit
      R2[ith]=constant_fit(R2_corr[ith],tmin,tmax);
      RT[ith]=constant_fit(RT_corr[ith],tmin,tmax);
      
      //plot
      out_R2<<write_constant_fit_plot(R2_corr[ith],R2[ith],tmin,tmax);
      out_RT<<write_constant_fit_plot(RT_corr[ith],RT[ith],tmin,tmax);
    }
  
  ////////////////////////////////////////// R1 //////////////////////////////////////  
  
  jvec R1_corr[nth];
  jvec R1(nth,njack);

  ofstream out_P(combine("plots/out_P_il_%d_ic_%d.xmg",il,ic).c_str());
  out_P<<"@type xydy"<<endl;
  ofstream out_D(combine("plots/out_D_il_%d_ic_%d.xmg",il,ic).c_str());
  out_D<<"@type xydy"<<endl;
  ofstream out_R1(combine("plots/out_R1_il_%d_ic_%d.xmg",il,ic).c_str());
  out_R1<<"@type xydy"<<endl;
  
  //load Pi and D
  jvec P_corr[nth],D_corr[nth];
  jvec ED(nth,njack),EP(nth,njack);
  for(int ith=0;ith<nth;ith++)
    {
      //load moving pion
      P_corr[ith]=load_2pts("2pts_P5P5.dat",il_sea,il,ith);
      out_P<<"@type xydy"<<endl;
      EP[ith]=constant_fit(effective_mass(P_corr[ith]),tmin_P,TH,combine("plots/P_eff_mass_il_%d_ic_%d_ith_%d.xmg",
									 il,ic,ith).c_str());
      out_P<<write_constant_fit_plot(effective_mass(P_corr[ith]),EP[ith],tmin_P,TH);
      out_P<<"&"<<endl;
      
      //recompute EP and ED from standing one
      if(ith)
	{
	  ED[ith]=latt_en(ED[0],th_P[ith]);
	  EP[ith]=latt_en(EP[0],th_P[ith]);
	}

      //load moving D
      D_corr[ith]=load_2pts("2pts_P5P5.dat",il,ic,ith);
      out_D<<"@type xydy"<<endl;
      ED[ith]=constant_fit(effective_mass(D_corr[ith]),tmin_D,TH,combine("plots/D_eff_mass_il_%d_ic_%d_ith_%d.xmg",
									 il,ic,ith).c_str());
      out_D<<write_constant_fit_plot(effective_mass(D_corr[ith]),ED[ith],tmin_D,TH);
      out_D<<"&"<<endl;
      
      //build the ratio
      R1_corr[ith]=lc0_th[ith]/lc0_th[0];
      for(int t=0;t<TH;t++)
	{
	  int E_fit_reco_flag=1;

	  jack Dt(njack),Pt(njack);	  
	  if(E_fit_reco_flag==0)
	    {
	      Dt=D_corr[0][t]/D_corr[ith][t];
	      Pt=P_corr[0][TH-t]/P_corr[ith][TH-t];
	    }
	  else
	    {
	      jack ED_th=latt_en(ED[0],th_P[ith]),EP_th=latt_en(EP[0],th_P[ith]);
	      Dt=exp(-(ED[0]-ED_th)*t)*ED_th/ED[0];
	      Pt=exp(-(EP[0]-EP_th)*(TH-t))*EP_th/EP[0];
	    }
	  
	  R1_corr[ith][t]*=Dt*Pt;
	}
      
      //fit
      R1[ith]=constant_fit(R1_corr[ith],tmin,tmax);
      
      //plot
      out_R1<<write_constant_fit_plot(R1_corr[ith],R1[ith],tmin,tmax);
    }
  
  //////////////////////////////////////// solve the ratios //////////////////////////////
  
  //compute f0[q2max]
  jvec f0_r(nth,njack),fP_r(nth,njack),fT_r(nth,njack);
  f0_r[0]=sqrt(R0*4*ED[0]*EP[0])/(ED[0]+EP[0]);
  cout<<"f0_r[q2max]: "<<f0_r[0]<<endl;
  
  //compute QK and Q2
  double mom[nth];
  jvec PK(nth,njack),QK(nth,njack);
  jvec P0(nth,njack),Q0(nth,njack),Q2(nth,njack),P2(nth,njack);
  jvec P0_r(nth,njack),Q0_r(nth,njack),Q2_r(nth,njack),P2_r(nth,njack);
  for(int ith=0;ith<nth;ith++)
    {
      P0[ith]=ED[ith]+EP[ith]; //P=initial+final
      Q0[ith]=ED[ith]-EP[ith]; //Q=initial-final
      P0_r[ith]=latt_en(ED[0],th_P[ith])+latt_en(EP[0],th_P[ith]);
      Q0_r[ith]=latt_en(ED[0],th_P[ith])-latt_en(EP[0],th_P[ith]);

      //we are describing the process D->Pi
      mom[ith]=momentum(th_P[ith]);
      double P_D=-mom[ith];
      double P_Pi=mom[ith];
  
      PK[ith]=P_D+P_Pi;
      QK[ith]=P_D-P_Pi;
      
      P2[ith]=sqr(P0[ith])-3*sqr(PK[ith]);
      Q2[ith]=sqr(Q0[ith])-3*sqr(QK[ith]);
      
      //reconstruct Q2
      P2_r[ith]=sqr(P0_r[ith])-3*sqr(PK[ith]);
      Q2_r[ith]=sqr(Q0_r[ith])-3*sqr(QK[ith]);
    }

  //checking Pion dispertion relation
  ofstream out_disp_P(combine("plots/Pion_disp_rel_il_%d_ic_%d.xmg",il,ic).c_str());
  out_disp_P<<"@type xydy"<<endl;
  for(int ith=0;ith<nth;ith++) out_disp_P<<3*sqr(mom[ith])<<" "<<sqr(EP[ith])<<endl;
  out_disp_P<<"&"<<endl;
  for(int ith=0;ith<nth;ith++) out_disp_P<<3*sqr(mom[ith])<<" "<<sqr(cont_en(EP[0],th_P[ith]))<<endl;
  out_disp_P<<"&"<<endl;
  for(int ith=0;ith<nth;ith++) out_disp_P<<3*sqr(mom[ith])<<" "<<sqr(latt_en(EP[0],th_P[ith]))<<endl;
  out_disp_P<<"&"<<endl;
  
  //checking D dispertion relation
  ofstream out_disp_D(combine("plots/D_disp_rel_il_%d_ic_%d.xmg",il,ic).c_str());
  out_disp_D<<"@type xydy"<<endl;
  for(int ith=0;ith<nth;ith++) out_disp_D<<3*sqr(mom[ith])<<" "<<sqr(ED[ith])<<endl;
  out_disp_D<<"&"<<endl;
  for(int ith=0;ith<nth;ith++) out_disp_D<<3*sqr(mom[ith])<<" "<<sqr(cont_en(ED[0],th_P[ith]))<<endl;
  out_disp_D<<"&"<<endl;
  for(int ith=0;ith<nth;ith++) out_disp_D<<3*sqr(mom[ith])<<" "<<sqr(latt_en(ED[0],th_P[ith]))<<endl;
  out_disp_D<<"&"<<endl;
  
  //compute xi
  jvec xi(nth,njack);
  for(int ith=1;ith<nth;ith++)
    {
      int E_fit_reco_flag=0; //it makes no diff
      
      jack P0_th=E_fit_reco_flag?P0_r[ith]:P0[ith];
      jack Q0_th=E_fit_reco_flag?Q0_r[ith]:Q0[ith];
      
      xi[ith]=R2[ith]*P0_th;
      xi[ith]/=QK[ith]-R2[ith]*Q0_th;
    }
  
  //compute fP
  ofstream out_fP_r(combine("plots/fP_r_il_%d_ic_%d.xmg",il,ic).c_str());
  out_fP_r<<"@type xydy"<<endl;
  for(int ith=1;ith<nth;ith++)
    {
      int E_fit_reco_flag=1; //it makes no diff
      
      jack P0_th=E_fit_reco_flag?P0_r[ith]:P0[ith];
      jack Q0_th=E_fit_reco_flag?Q0_r[ith]:Q0[ith];
      
      jack c=P0_th/(ED[0]+EP[0])*(1+xi[ith]*Q0_th/P0_th);
      fP_r[ith]=R1[ith]/c*f0_r[0];

      out_fP_r<<Q2[ith].med()<<" "<<fP_r[ith]<<endl;
    }
  
  //compute f0 and fT
  ofstream out_f0_r(combine("plots/f0_r_il_%d_ic_%d.xmg",il,ic).c_str());
  ofstream out_fT_r(combine("plots/fT_r_il_%d_ic_%d.xmg",il,ic).c_str());;
  out_f0_r<<"@type xydy"<<endl;
  out_f0_r<<Q2[0].med()<<" "<<f0_r[0]<<endl;
  out_fT_r<<"@type xydy"<<endl;
  for(int ith=1;ith<nth;ith++)
    {
      //it seems better here to solve using reconstructed energies
      int E_fit_reco_flag=0;
  
      jack EP_th=E_fit_reco_flag?latt_en(EP[0],th_P[ith]):EP[ith];
      jack ED_th=E_fit_reco_flag?latt_en(ED[0],th_P[ith]):ED[ith];
      jack Q2_th=E_fit_reco_flag?Q2_r[ith]:Q2[ith];
      
      jack fM_r=xi[ith]*fP_r[ith]; //checked
      f0_r[ith]=fP_r[ith]+fM_r[ith]*Q2_th/(sqr(ED_th)-sqr(EP_th));
      
      out_f0_r<<Q2[ith].med()<<" "<<f0_r[ith]<<endl;
      
      fT_r[ith]=fM_r[ith]*RT[ith]*Zt_med[ibeta]/Zv_med[ibeta]*(EP[0]+ED[0])/(ED[ith]+EP[ith]); //ADD
      
      out_fT_r<<Q2[ith].med()<<" "<<fT_r[ith]<<endl;
    }
  
  
  //////////////////////////////////////// analytic method /////////////////////////////  
  
  jvec fP_a(nth,njack),fM_a(nth,njack),f0_a(nth,njack),fT_a(nth,njack);
  jvec fP_n(nth,njack),fM_n(nth,njack),f0_n(nth,njack),fT_n(nth,njack);
  
  //determine M and Z for pion and D
  jvec ZP(nth,njack),ZD(nth,njack);
  for(int ith=0;ith<nth;ith++)
    {
      jack E,Z2;
      two_pts_fit(E,Z2,P_corr[ith],tmin_P,TH);
      ZP[ith]=sqrt(Z2);
      two_pts_fit(E,Z2,D_corr[ith],tmin_D,TH);
      ZD[ith]=sqrt(Z2);
    }
  
  //compute V
  jvec VK_a(nth,njack),V0_a(nth,njack),TK_a(nth,njack);
  jvec VK_n(nth,njack),V0_n(nth,njack),TK_n(nth,njack);
  for(int ith=0;ith<nth;ith++)
    {
      ofstream out_V0(combine("plots/V0_il_%d_ic_%d_ith_%d_analytic_numeric.xmg",il,ic,ith).c_str());
      out_V0<<"@type xydy"<<endl;
      ofstream out_VK(combine("plots/VK_il_%d_ic_%d_ith_%d_analytic_numeric.xmg",il,ic,ith).c_str());
      out_VK<<"@type xydy"<<endl;
      ofstream out_TK(combine("plots/TK_il_%d_ic_%d_ith_%d_analytic_numeric.xmg",il,ic,ith).c_str());
      out_TK<<"@type xydy"<<endl;
      ofstream out_dt(combine("plots/dt_il_%d_ic_%d_ith_%d.xmg",il,ic,ith).c_str());
      out_dt<<"@type xydy"<<endl;
      
      //computing time dependance
      jvec dt_a(TH+1,njack),dt_n(TH+1,njack);
      {
	//it seems better here to use fitted energies
	int E_fit_reco_flag=1;
	jack EP_th=E_fit_reco_flag?latt_en(EP[0],th_P[ith]):EP[ith];
	jack ED_th=E_fit_reco_flag?latt_en(ED[0],th_P[ith]):ED[ith];
	
	for(int t=0;t<=TH;t++)
	  {
	    dt_a[t]=exp(-(ED_th*t+EP_th*(TH-t)))*ZP[0]*ZD[0]/(4*EP_th*ED_th);
	    dt_n[t]=D_corr[ith][t]*P_corr[ith][TH-t]/(ZD[0]*ZP[0]);
	  }
      }
      
      //remove time dependance using analytic or numeric expression
      jvec VK_corr_a=Zv_med[ibeta]*lcK_th[ith]/dt_a,V0_corr_a=Zv_med[ibeta]*lc0_th[ith]/dt_a;
      jvec VK_corr_n=Zv_med[ibeta]*lcK_th[ith]/dt_n,V0_corr_n=Zv_med[ibeta]*lc0_th[ith]/dt_n;
      jvec TK_corr_n=Zt_med[ibeta]*lcT_th[ith]/dt_n,TK_corr_a=Zt_med[ibeta]*lcT_th[ith]/dt_a;
      
      //fit V0
      V0_a[ith]=constant_fit(V0_corr_a,tmin,tmax);
      V0_n[ith]=constant_fit(V0_corr_n,tmin,tmax);
      out_V0<<write_constant_fit_plot(V0_corr_a,V0_a[ith],tmin,tmax)<<"&"<<endl;
      out_V0<<write_constant_fit_plot(V0_corr_n,V0_n[ith],tmin,tmax)<<"&"<<endl;
      
      //fit VK
      VK_a[ith]=constant_fit(VK_corr_a,tmin,tmax);
      VK_n[ith]=constant_fit(VK_corr_n,tmin,tmax);
      out_VK<<write_constant_fit_plot(VK_corr_a,VK_a[ith],tmin,tmax)<<"&"<<endl;
      out_VK<<write_constant_fit_plot(VK_corr_n,VK_n[ith],tmin,tmax)<<"&"<<endl;

      //fit TK
      TK_a[ith]=constant_fit(TK_corr_a,tmin,tmax);
      TK_n[ith]=constant_fit(TK_corr_n,tmin,tmax);
      out_TK<<write_constant_fit_plot(TK_corr_a,TK_a[ith],tmin,tmax)<<"&"<<endl;
      out_TK<<write_constant_fit_plot(TK_corr_n,TK_n[ith],tmin,tmax)<<"&"<<endl;
    }
  
  //compute f0(q2max)
  f0_a[0]=V0_a[0]/(ED[0]+EP[0]);
  f0_n[0]=V0_n[0]/(ED[0]+EP[0]);
  cout<<"f0_a["<<Q2[0].med()<<"]: "<<f0_a[0]<<endl;
  cout<<"f0_n["<<Q2[0].med()<<"]: "<<f0_n[0]<<endl;
  
  //solve for fP and f0
  for(int ith=1;ith<nth;ith++)
    {
      jack delta=P0[ith]*QK[ith]-Q0[ith]*PK[ith];

      //solve using analytic fit
      jack deltaP_a=V0_a[ith]*QK[ith]-Q0[ith]*VK_a[ith];
      jack deltaM_a=P0[ith]*VK_a[ith]-V0_a[ith]*PK[ith];  
      fP_a[ith]=deltaP_a/delta;
      fM_a[ith]=deltaM_a/delta;
      
      //solve using numeric fit
      jack deltaP_n=V0_n[ith]*QK[ith]-Q0[ith]*VK_n[ith];
      jack deltaM_n=P0[ith]*VK_n[ith]-V0_n[ith]*PK[ith];  
      fP_n[ith]=deltaP_n/delta;
      fM_n[ith]=deltaM_n/delta;

      //compute f0
      f0_a[ith]=fP_a[ith]+fM_a[ith]*Q2[ith]/(ED[0]*ED[0]-EP[0]*EP[0]);
      f0_n[ith]=fP_n[ith]+fM_n[ith]*Q2[ith]/(ED[0]*ED[0]-EP[0]*EP[0]);

      //solve fT
      fT_a[ith]=-TK_a[ith]*(EP[0]+ED[0])/(2*(ED[ith]+EP[ith]))/mom[ith];
      fT_n[ith]=-TK_n[ith]*(EP[0]+ED[0])/(2*(ED[ith]+EP[ith]))/mom[ith];
    }
  
  //write analytic and umeric plot of fP and f0
  ofstream out_fP_a("plots/fP_a.xmg"),out_fP_n("plots/fP_n.xmg");
  ofstream out_fM_a("plots/fM_a.xmg"),out_fM_n("plots/fM_n.xmg");
  ofstream out_f0_a("plots/f0_a.xmg"),out_f0_n("plots/f0_n.xmg");
  ofstream out_fT_a("plots/fT_a.xmg"),out_fT_n("plots/fT_n.xmg");
  out_fP_a<<"@type xydy"<<endl;
  out_fP_n<<"@type xydy"<<endl;
  out_f0_a<<"@type xydy"<<endl;
  out_f0_n<<"@type xydy"<<endl;
  out_fM_a<<"@type xydy"<<endl;
  out_fM_n<<"@type xydy"<<endl;
  out_fT_a<<"@type xydy"<<endl;
  out_fT_n<<"@type xydy"<<endl;
  out_f0_a<<Q2[0].med()<<" "<<f0_a[0]<<endl;
  out_f0_n<<Q2[0].med()<<" "<<f0_n[0]<<endl;
  for(int ith=1;ith<nth;ith++)
    {
      out_fP_a<<Q2[ith].med()<<" "<<fP_a[ith]<<endl;
      out_fP_n<<Q2[ith].med()<<" "<<fP_n[ith]<<endl;
      out_fM_a<<Q2[ith].med()<<" "<<fM_a[ith]<<endl;
      out_fM_n<<Q2[ith].med()<<" "<<fM_n[ith]<<endl;
      out_f0_a<<Q2[ith].med()<<" "<<f0_a[ith]<<endl;
      out_f0_n<<Q2[ith].med()<<" "<<f0_n[ith]<<endl;
      out_fT_a<<Q2[ith].med()<<" "<<fT_a[ith]<<endl;
      out_fT_n<<Q2[ith].med()<<" "<<fT_n[ith]<<endl;
    }
  
  ext_EP=EP;
  ext_ED=ED;
  ext_Q2=Q2;
  ext_fP=fP_a;
  ext_fM=fM_a;
  ext_f0=f0_a;
  ext_fT=fT_a;
}

jvec parab_spline(jvec *A,double *x,double xint,const char *path=NULL)
{
  jvec out(nth,njack);
  ofstream plot;
  if(path!=NULL) plot.open(path);
  
  jvec a(nth,njack),b(nth,njack),c(nth,njack);
  for(int ith=0;ith<nth;ith++)
    {
      for(int ijack=0;ijack<=njack;ijack++)
	{
	  double y[3];
	  for(int i=0;i<3;i++) y[i]=A[i][ith][ijack];
	  parabolic_spline(a[ith].data[ijack],b[ith].data[ijack],c[ith].data[ijack],x,y);
	  out[ith].data[ijack]=a[ith][ijack]*xint*xint+b[ith][ijack]*xint+c[ith][ijack];
	}
    }
  
  if(path)
    {
      for(int ith=0;ith<nth;ith++)
	{
	  jvec inte(100,njack);
	  double xin=min(xint,x[0]),xfin=max(xint,x[2]),dx=(xfin-xin)/99;
	  double xn[100];
	  for(int ix=0;ix<100;ix++)
	    {
	      xn[ix]=xin+dx*ix;
	      inte[ix]=a[ith]*xn[ix]*xn[ix]+b[ith]*xn[ix]+c[ith];
	    }
	  plot<<write_polygon(xn,inte);
	  plot<<"&"<<endl;
	  plot<<"@type xydy"<<endl;
	  for(int ic=0;ic<3;ic++) plot<<x[ic]<<" "<<A[ic][ith]<<endl;
	  plot<<"&"<<endl;
	  plot<<xint<<" "<<out[ith]<<endl;
	  plot<<"&"<<endl;
	}
    }
  return out;
}

void read_input()
{
  FILE *fin=open_file("analysis_pars","r");
  read_formatted_from_file_expecting((char*)(&T),fin,"%d","T");
  read_formatted_from_file_expecting((char*)(&ibeta),fin,"%d","ibeta");
  TH=L=T/2;
  read_formatted_from_file_expecting((char*)(th_P),fin,"%lg","th");
  for(int ith=1;ith<nth;ith++) read_formatted_from_file((char*)(th_P+ith),fin,"%lg",combine("th%d",ith).c_str());
  read_formatted_from_file_expecting((char*)(xc),fin,"%lg","xc");
  for(int ic=1;ic<3;ic++) read_formatted_from_file((char*)(xc+ic),fin,"%lg",combine("xc%d",ic).c_str());
  read_formatted_from_file_expecting((char*)(&xc_phys),fin,"%lg","xc_phys");
  read_formatted_from_file_expecting((char*)(il+0),fin,"%d","il_list");
  read_formatted_from_file((char*)(il+1),fin,"%d","il[1]");
  read_formatted_from_file((char*)(il+2),fin,"%d","il[2]");
  read_formatted_from_file((char*)(il+3),fin,"%d","il[3]");
  read_formatted_from_file_expecting((char*)(xs),fin,"%lg","xs");
  for(int is=1;is<3;is++) read_formatted_from_file((char*)(xs+is),fin,"%lg",combine("xs%d",is).c_str());
  read_formatted_from_file_expecting((char*)(&xs_phys),fin,"%lg","xs_phys");
  read_formatted_from_file_expecting((char*)(&tmin),fin,"%d","tint");
  read_formatted_from_file((char*)(&tmax),fin,"%d","tint");
  
  read_formatted_from_file_expecting((char*)(&tmin_P),fin,"%d","tmin_P");  
  read_formatted_from_file_expecting((char*)(&tmin_D),fin,"%d","tmin_D");
  
  fclose(fin);
}

void study_il(jvec &EP,jvec &MD,jvec &Q2,jvec &fP,jvec &fM,jvec &f0,jvec &fT,int il_sea,int il)
{
  int ic_base=8;
  jvec EP_ic[3],MD_ic[3],Q2_ic[3],fP_ic[3],fM_ic[3],f0_ic[3],fT_ic[3];
  for(int dic=2;dic>=0;dic--)
    {
      int ic=ic_base+dic;
      Q2_ic[dic]=fP_ic[dic]=f0_ic[dic]=fT_ic[dic]=jvec(nth,njack);
      extr(EP_ic[dic],MD_ic[dic],Q2_ic[dic],fP_ic[dic],fM_ic[dic],f0_ic[dic],fT_ic[dic],il_sea,il,ic);
    }
  
  EP=parab_spline(EP_ic,xc,xc_phys,combine("plots/EP_spline_%d.xmg",il).c_str());
  MD=parab_spline(MD_ic,xc,xc_phys,combine("plots/MD_spline_%d.xmg",il).c_str());
  Q2=parab_spline(Q2_ic,xc,xc_phys,combine("plots/Q2_spline_%d.xmg",il).c_str());
  fP=parab_spline(fP_ic,xc,xc_phys,combine("plots/fP_spline_%d.xmg",il).c_str());
  fM=parab_spline(fM_ic,xc,xc_phys,combine("plots/fM_spline_%d.xmg",il).c_str());
  f0=parab_spline(f0_ic,xc,xc_phys,combine("plots/f0_spline_%d.xmg",il).c_str());
  fT=parab_spline(fT_ic,xc,xc_phys,combine("plots/fT_spline_%d.xmg",il).c_str());
  fP[0]*=0;
  fT[0]*=0;
  
  ofstream out_fP(combine("plots/fP_%d.xmg",il).c_str());
  ofstream out_fM(combine("plots/fM_%d.xmg",il).c_str());
  ofstream out_f0(combine("plots/f0_%d.xmg",il).c_str());
  ofstream out_fT(combine("plots/fT_%d.xmg",il).c_str());
  out_fP<<"@type xydy"<<endl;
  out_fM<<"@type xydy"<<endl;
  out_f0<<"@type xydy"<<endl;
  out_fT<<"@type xydy"<<endl;
  for(int ith=0;ith<nth;ith++)
    {
      if(ith) out_fP<<Q2[ith].med()<<" "<<fP[ith]<<endl;
      if(ith) out_fM<<Q2[ith].med()<<" "<<fM[ith]<<endl;
      out_f0<<Q2[ith].med()<<" "<<f0[ith]<<endl;
      if(ith) out_fT<<Q2[ith].med()<<" "<<fT[ith]<<endl;
    }
}

void write(const char *dir,jvec &EP,jvec &MD,jvec &Q2,jvec fP,jvec &fM,jvec &f0,jvec &fT)
{
  char results_path[100];
  sprintf(results_path,"%s/EP_Q2_fP_f0_fT_f0s",dir);
  EP.write_to_binfile(results_path);
  Q2.append_to_binfile(results_path);
  fP.append_to_binfile(results_path);
  f0.append_to_binfile(results_path);
  fT.append_to_binfile(results_path);
  f0.append_to_binfile(results_path); //fake
  
  MD[0].write_to_binfile(combine("%s/M_H_P5",dir).c_str());
  EP[0].write_to_binfile(combine("%s/M_K_P5",dir).c_str());
}

int main()
{
  read_input();
  
  jvec EP_l[4],MD_l[4],Q2_l[4],fP_l[4],fM_l[4],f0_l[4],fT_l[4];
  for(int i=0;i<4;i++) study_il(EP_l[i],MD_l[i],Q2_l[i],fP_l[i],fM_l[i],f0_l[i],fT_l[i],il[0],il[i]);
  
  write("D_Pi",EP_l[0],MD_l[0],Q2_l[0],fP_l[0],fM_l[0],f0_l[0],fT_l[0]);

  {
  ofstream out_fP("plots/fP.xmg");
  ofstream out_fM("plots/fM.xmg");
  ofstream out_f0("plots/f0.xmg");
  ofstream out_fT("plots/fT.xmg");
  out_fP<<"@type xydy"<<endl;
  out_fM<<"@type xydy"<<endl;
  out_f0<<"@type xydy"<<endl;
  out_fT<<"@type xydy"<<endl;
  for(int ith=0;ith<nth;ith++)
    {
      if(ith) out_fP<<Q2_l[0][ith].med()<<" "<<fP_l[0][ith]<<endl;
      if(ith) out_fM<<Q2_l[0][ith].med()<<" "<<fM_l[0][ith]<<endl;
      out_f0<<Q2_l[0][ith].med()<<" "<<f0_l[0][ith]<<endl;
      if(ith) out_fT<<Q2_l[0][ith].med()<<" "<<fT_l[0][ith]<<endl;
    }
  }
  /////
  
  jvec EP,MD,Q2,fP,fM,f0,fT;
  EP=parab_spline(EP_l+1,xs,xs_phys,"plots/EP_spline_s.xmg");
  MD=parab_spline(MD_l+1,xs,xs_phys,"plots/MD_spline_s.xmg");
  Q2=parab_spline(Q2_l+1,xs,xs_phys,"plots/Q2_spline_s.xmg");
  fP=parab_spline(fP_l+1,xs,xs_phys,"plots/fP_spline_s.xmg");
  fM=parab_spline(fM_l+1,xs,xs_phys,"plots/fM_spline_s.xmg");
  f0=parab_spline(f0_l+1,xs,xs_phys,"plots/f0_spline_s.xmg");
  fT=parab_spline(fT_l+1,xs,xs_phys,"plots/fT_spline_s.xmg");
  fP[0]*=0;
  fT[0]*=0;
  
  ofstream out_fP("plots/fP_s.xmg");
  ofstream out_fM("plots/fM_s.xmg");
  ofstream out_f0("plots/f0_s.xmg");
  ofstream out_fT("plots/fT_s.xmg");
  out_fP<<"@type xydy"<<endl;
  out_fM<<"@type xydy"<<endl;
  out_f0<<"@type xydy"<<endl;
  out_fT<<"@type xydy"<<endl;
  for(int ith=0;ith<nth;ith++)
    {
      if(ith) out_fP<<Q2[ith].med()<<" "<<fP[ith]<<endl;
      if(ith) out_fM<<Q2[ith].med()<<" "<<fM[ith]<<endl;
      out_f0<<Q2[ith].med()<<" "<<f0[ith]<<endl;
      if(ith) out_fT<<Q2[ith].med()<<" "<<fT[ith]<<endl;
    }
  
  write("D_K",EP,MD,Q2,fP,fM,f0,fT);

  return 0;
}
