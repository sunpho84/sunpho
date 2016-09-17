#include "../../src/include.h"

int T=48,TH=T/2,L=24;
int spat_vol=L*L*L;
int tmin=11;
int tmax=TH-tmin;
#define NPROJ 1
double Zv=0.6103;
double aMu=0.0464;
int njacks=22;
int phi=0,eta=1,alt=2,nphi_eta_alt=3;
int nw=8;
int norie=2;
int nr=2;

const double bc=+0.4416405016774463;
const double en_ph=2*M_PI/L;
const double pi_ext_mu=M_PI*bc/L;
const double pi_int_mu=pi_ext_mu-en_ph;
const double en_nu=sqrt(3)*pi_ext_mu;
const double m_mu=0.0464;
const double en_ext_mu=sqrt(sqr(m_mu)+3*sqr(pi_ext_mu));
const double en_int_mu=sqrt(sqr(m_mu)+2*sqr(pi_ext_mu)+sqr(pi_int_mu));
const double m_pi=0.210;
const double e_m=sqrt(sqr(m_pi)+sqr(en_ph));

////////////////////////////// ib on quarks only ///////////////////////

const int RE=0,IM=1;
const int LI=0;
int nmass=1;

//insertion is done on the second passed mass
int icombo_qonly(int imrev,int imins,int r,int REIM)
{return REIM+2*(r+nr*(imins+nmass*imrev));}

jvec load(const char *name,int im1,int im2,int REIM,int rpar,int par)
{
  jvec a=jvec_load(combine("corr_%s",name).c_str(),T,njacks,icombo_qonly(im1,im2,0,REIM)).simmetrized(par);
  if(nr==1 || rpar==0) return a;
  else return (a+rpar*jvec_load(combine("corr_%s",name).c_str(),T,njacks,icombo_qonly(im1,im2,1,REIM)).simmetrized(par))/2;
}

///////////////////////////////////////////////////////////////////////////

int icombo(int ri,int iproj,int iw,int rl,int orie,int r2,int phi_eta_alt,int irev,int qins,int il=0)
{return ri+2*(iproj+NPROJ*(iw+nw*(rl+nr*(orie+norie*(r2+nr*(phi_eta_alt+nphi_eta_alt*(irev+2*(qins+2*il))))))));}

//int icombo_chris(int ri,int tmu,int iproj,int iw,int orie,int phi_eta,int irev,int qins,int il=0)
//{return ri+2*(tmu+T*(iproj+NPROJ*(iw+nw*(orie+2*(phi_eta+2*(irev+2*(qins+2*il)))))));}

jvec load_alone(int i)
{return jvec_load("corr_hl",T,njacks,i);}

//jvec load_alone_chris(int i)
//{return jvec_load("corr_hl_chris",T,njacks,i);}

jvec load(int ibase)
{
  jvec tot(T,njacks);
  
  tot=0;
  int n=8;
  for(int i=0;i<n;i++)
    {
      tot+=load_alone(ibase+(i+0)*16);
      tot-=load_alone(ibase+(i+8)*16);
    }
  return tot.simmetrized(+1)/2/n;
}

const int r_unk[2]={+2,+0};
const int r_evn[2]={+1,+1};
const int r_odd[2]={+1,-1};

jvec load_improved(int iw,int iproj,int ri,const int *sign_r,int phi_eta_alt,int il=0,const char *path=NULL)
{
  jvec out(T,njacks);
  out=0;
  
  int sign_q[2]={+1,-1};
  int n=0;
  int lim=2;
  for(int qins=0;qins<1;qins++)
    for(int irev=0;irev<1;irev++) //tbf
      for(int r2=0;r2<nr;r2++)
	for(int orie=0;orie<lim;orie++)
	  for(int rl=0;rl<nr;rl++)
	    //if(r2==rl)
	    {
	      int ic=icombo(ri,iproj,iw,rl,orie,r2,phi_eta_alt,irev,qins,il);
	      //if(debug_load) cout<<"Loading combo: "<<ic<<endl;
	      jvec corr=load_alone(ic);
	      out+=sign_q[qins]*sign_r[rl]*sign_r[r2]*corr;
	      if(path)
		{
		  ostringstream app;
		  app<<"qins_"<<qins<<"_irev_"<<irev<<"_phietaalt_"<<phi_eta_alt<<"_r2_"<<r2<<"_orie_"<<orie<<"_rl_"<<rl;
		  ofstream fout(combine(path,app.str().c_str()).c_str());
		  fout<<"@type xydy"<<endl;
		  fout<<corr<<endl;
		  fout<<"@title \""<<ic<<"\""<<endl;
		}
	      n++;
	    }
  
  //return out.simmetrized(-1)/n;
  return out.subset(0,TH+1)/n;
  //return out/n;
}

// jvec load_improved_chris(int tmu,int iw,int iproj,int ri,const int *sign_r,int il=0,const char *path=NULL)
// {
//   jvec out(T,njacks);
//   out=0;
  
//   int sign_q[2]={+1,-1};
//   int sign_rev[2]={+1,-1};
//   int n=0;
//   int lim=2;
//   for(int qins=0;qins<1;qins++)
//     for(int irev=0;irev<1;irev++) //tbf
//       for(int phi_eta=0;phi_eta<1;phi_eta++) //eta is more noisy
// 	for(int orie=1;orie<2;orie++)
// 	  {
// 	    int ic=icombo_chris(ri,tmu,iproj,iw,orie,phi_eta,irev,qins,il);
// 	    //if(debug_load) cout<<"Loading combo: "<<ic<<endl;
// 	    jvec corr=load_alone_chris(ic);
// 	    out+=sign_rev[irev]*sign_q[qins]*corr;
// 	    if(path)
// 	      {
// 		ostringstream app;
// 		app<<"tmu_"<<tmu<<"_qins_"<<qins<<"_irev_"<<irev<<"_phieta_"<<phi_eta<<"_orie_"<<orie;
// 		ofstream fout(combine(path,app.str().c_str()).c_str());
// 		fout<<"@type xydy"<<endl;
// 		fout<<corr<<endl;
// 		fout<<"@title \""<<ic<<"\""<<endl;
// 	      }
// 	    n++;
// 	  }
  
//   //return out.simmetrized(-1)/n;
//   //return out.subset(0,TH+1)/n;
//   return out/n;
// }

// jvec load_improved_chris(int iw,int iproj,int ri,const int *sign_r,int icase,int int12,const char *path=NULL,int il=0,const char *path_int=NULL)
// {
//   jvec *temp=new jvec [T];
//   for(int t2=0;t2<T;t2++) temp[t2]=load_improved_chris(t2,iw,iproj,ri,sign_r,il,path_int);
  
//   jvec out(T,njacks);
//   out=0;
//   const int tw=14;
//   int tfix;
  
//   int diff=2;
//   if(int12==4) diff=-diff;
//   if(icase<=3) tfix=tw-diff;
//   else         tfix=tw+diff;
  
//   for(int t1=0;t1<T;t1++)
//     for(int t2=0;t2<T;t2++)
//       {
// 	bool cons;
// 	switch(icase)
// 	  {
// 	  case 1:cons=(tw>=t1 && t2>=t1 && t2>=tw);break; // t1 tw t2
// 	  case 2:cons=(tw>=t1 && t2>=t1 && tw>=t2);break; // t1 t2 tw
// 	  case 3:cons=(tw>=t1 && t1>=t2 && tw>=t2);break; // t2 t1 tw
// 	  case 4:cons=(t1>=tw && t2>=t1 && t2>=tw);break; // tw t1 t2
// 	  case 5:cons=(t1>=tw && t1>=t2 && t2>=tw);break; // tw t2 t1
// 	  case 6:cons=(t1>=tw && t1>=t2 && tw>=t2);break; // t2 tw t1
// 	  default:cons=false;crash("unknown case %d",icase);
// 	  }
	
// 	int t3;
// 	switch(int12)
// 	  {
// 	  case 1:t3=t2;break;
// 	  case 2:t3=t1;break;
// 	  case 3:t3=t2;break;
// 	  case 4:t3=t1;break;
// 	  default:t3=0;crash("unknwon integration time %d",int12);
// 	  }
	
// 	if(int12==3) cons&=(t1==tfix);
// 	if(int12==4) cons&=(t2==tfix);
	
// 	if(cons) out[t3]+=temp[t2][t1];
//       }
  
//   //out=0;
//   //for(int t2=0;t2<T;t2++)
//   //for(int t1=0;t1<14;t1++)
//   //out[t2]+=temp[t2][t1];
  
//   //for(int t2=0;t2<T;t2++)
//   //for(int t1=0;t1<T;t1++)
//   //out[t1]+=temp[t2][t1];
  
//   //for(int t1=0;t1<T;t1++)
//   //for(int t2=t1;t2<T;t2++)
//   //out[t1]+=temp[t2][t1];
  
//   //for(int t1=0;t1<T;t1++)
//   //for(int t2=14;t2<t1;t2++)
//   //out[t1]+=temp[t2][t1];
  
//   //int t1=14+6;
//   //for(int t2=0;t2<T;t2++)
//   //for(int t3=t1;t3<t2;t3++)
//   //out[t2]+=temp[t3][t1];
  
//   //int t1=16;
//   //for(int t2=0;t2<T;t2++)
//   //out[t2]+=temp[t2][t1];
  
//   if(path)
//     {
//       //tags
//       const char tags_order[6][40]={
// 	"t\\s1\\N<=t\\sw\\N<=t\\s2\\N",
// 	"t\\s1\\N<=t\\s2\\N<=t\\sw\\N",
// 	"t\\s2\\N<=t\\s1\\N<=t\\sw\\N",
// 	"t\\sw\\N<=t\\s1\\N<=t\\s2\\N",
// 	"t\\sw\\N<=t\\s2\\N<=t\\s1\\N",
// 	"t\\s2\\N<=t\\sw\\N<=t\\s1\\N"};
      
//       ofstream fout(combine(path,icase,int12).c_str());
//       fout<<"@type xydy"<<endl;
//       fout<<out<<endl;
//       fout<<"@title \""<<tags_order[icase-1]<<"\""<<endl;
//       if(int12<3)
// 	{
// 	  fout<<"@subtitle \"integrating t\\s"<<int12<<"\""<<endl;
// 	  fout<<"@xaxis label \"t\\s"<<3-int12<<"\""<<endl;
// 	}
//       else
// 	if(int12==3)
// 	  {
// 	    fout<<"@subtitle \"fixing t1=tw+"<<tfix-tw<<"\""<<endl;
// 	    fout<<"@xaxis label \"t2\""<<endl;
// 	  }
// 	else
// 	  {
// 	    fout<<"@subtitle \"fixing t2=tw+"<<tfix-tw<<"\""<<endl;
// 	    fout<<"@xaxis label \"t1\""<<endl;
// 	  }
      
//       fout<<"@xaxis label char size 1.5"<<endl;
//     }
  
//   return out;
// }

double *c_fit,*e_fit;
template <class T> T fun_case_5_int_t2(T ca,T cb,double t,T a,T b)
{
  double dt=t-14;
  return ca*exp(-dt*a)-cb*exp(-dt*b);
}

void ch2_case_5_int_t2_fit(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  double ca=p[0];
  double cb=p[1];
  double a=p[2];
  double b=p[3];
  
  for(int t=14;t<T;t++)
    {
      double num=c_fit[t];
      double teo=fun_case_5_int_t2(ca,cb,t,a,b);
      double diff=num-teo;
      double err=e_fit[t];
      double cont=sqr(diff/err);
      ch+=cont;
    }
}

void case_5_int_t2_fit(jack &CA,jack &CB,jack &A,jack &B,jvec corr)
{
  TMinuit minu;
  //minu.SetPrintLevel(-1);
  minu.SetFCN(ch2_case_5_int_t2_fit);

  minu.DefineParameter(0,"CA",7e-06,0.001,0,0);
  minu.DefineParameter(1,"CB",7e-06,0.001,0,0);
  minu.DefineParameter(2,"A",0.2,0.001,0,0);
  minu.DefineParameter(3,"B",0.2,0.001,0,0);
  
  c_fit=new double[T];
  e_fit=new double[T];
  
  for(int iel=0;iel<T;iel++)
    e_fit[iel]=corr[iel].err();
  
  for(int ijack_fit=0;ijack_fit<=njacks;ijack_fit++)
    {
      //minu.FixParameter(0);
      for(int iel=0;iel<=TH;iel++)
          c_fit[iel]=corr[iel][ijack_fit];
      minu.Migrad();
      double dum;
      minu.GetParameter(0,CA.data[ijack_fit],dum);
      minu.GetParameter(1,CB.data[ijack_fit],dum);
      minu.GetParameter(2,A.data[ijack_fit],dum);
      minu.GetParameter(3,B.data[ijack_fit],dum);
    }
  
  double ch2,grad[4],par[4]={CA[njacks],CB[njacks],A[njacks],B[njacks]};
  minu.Eval(4,grad,ch2,par,4);
  cout<<"A: "<<A<<" B: "<<B<<", ch2: "<<ch2<<endl;
  
  jvec fun_fit(T,njacks);
  fun_fit=0;
  for(int t=14;t<T;t++) fun_fit[t]=fun_case_5_int_t2(CA,CB,t,A,B);
  ofstream fout("fit_case5_int2.xmg");
  fout<<"@type xydy"<<endl;
  fout<<corr<<"&"<<endl<<"@type xy"<<endl;
  for(int t=14;t<T;t++) fout<<t<<" "<<fun_fit[t].med()-fun_fit[t].err()<<endl;
  for(int t=T-1;t>=14;t--) fout<<t<<" "<<fun_fit[t].med()+fun_fit[t].err()<<endl;
  
  delete[] c_fit;
  delete[] e_fit;
}

int main()
{
  //debug_load=false;
  debug_fit=false;
  
  //pure pion
  jvec P_corr_0=load("P5P5_00",LI,LI,RE,1,1);
  jack M_P=constant_fit(effective_mass(P_corr_0),12,23,"/tmp/pion_effmass.xmg");
  
  //pion exchange
  jvec P_corr_AB=load("P5P5_AB",LI,LI,RE,1,1);
  jvec P_ratio_AB=P_corr_AB/P_corr_0;
  P_ratio_AB.print_to_file("pion_exchange.xmg");
  
  jvec P_corr_LL=load("P5P5_LL",LI,LI,RE,1,1);
  jvec P_ratio_LL=P_corr_LL/P_corr_0;
  P_ratio_LL.print_to_file("pion_exchange_sqrt.xmg");
  
  //pion self
  jvec P_corr_0X=load("P5P5_0X",LI,LI,RE,1,1);
  jvec P_ratio_0X=P_corr_0X/P_corr_0;
  P_ratio_0X.print_to_file("pion_self.xmg");
  
  jvec P_corr_0M=load("P5P5_0M",LI,LI,RE,1,1);
  jvec P_ratio_0M=P_corr_0M/P_corr_0;
  P_ratio_0M.print_to_file("pion_self_sqrt.xmg");
  
  ifstream fin("../pure_L");
  double L[T];
  for(int t=0;t<T;t++)
    fin>>L[t];
  fin.close();
  
  //compute normalization
  jvec cP5P5=jvec_load("corr_P5P5_00",T,njacks,0).simmetrized(1);
  jvec cA0P5=-jvec_load("corr_A0P5_00",T,njacks,0).simmetrized(-1);
  jack aM(njacks),ZP5(njacks),ZA0(njacks);
  cP5P5.print_to_file("corr_P5P5.xmg");
  two_pts_SL_fit(aM,ZA0,ZP5,cA0P5,cP5P5,tmin,TH-1,tmin,TH,"eff_mass_P5P5.xmg",NULL,"/tmp/list.txt",-1,+1);
  jack aMbis(njacks);
  jack Z2(njacks);
  two_pts_fit(aMbis,Z2,cP5P5,tmin,TH);
  double ml=0.0100,a=0.087/0.197;
  jack afpi_P5=2*ml*ZP5/(aM*aM);
  jack norm=2*aM/ZP5*afpi_P5;
  jack afpi_A0=ZA0/aM*Zv;
  cout<<" fpi from P5: "<<afpi_P5/a<<endl;
  cout<<" fpi from A0: "<<afpi_A0/a<<endl;
  cout<<" aMpi: "<<smart_print(aM)<<" or "<<smart_print(aMbis)<<" from p5p5 only "<<endl;
  
  //try leading order
  jack A02=2*sqr(afpi_P5/a*aM/a*aMu/a)*(1-sqr(aMu/aM));
  ofstream lead("leading_order.xmg");
  lead<<" @s0  symbol 1"<<endl;
  lead<<" @s0  symbol color 14"<<endl;
  lead<<" @s0  symbol fill color 14"<<endl;
  lead<<" @s0  symbol linewidth 2.0"<<endl;
  lead<<" @s0  line type 0"<<endl;
  lead<<" @s0  line color 14"<<endl;
  lead<<" @s0  errorbar color 14"<<endl;
  lead<<" @s0  errorbar linewidth 2.0"<<endl;
  lead<<" @s0  errorbar riser linewidth 2.0"<<endl;
  lead<<" @s1  symbol color 10"<<endl;
  lead<<" @s1  symbol fill color 10"<<endl;
  lead<<" @s1  line color 10"<<endl;
  lead<<" @s1  fill type 1"<<endl;
  lead<<" @s1  fill color 10"<<endl;
  lead<<" @s1  errorbar color 10"<<endl;
  lead<<" @s2  symbol color 2"<<endl;
  lead<<" @s2  symbol fill color 2"<<endl;
  lead<<" @s2  line linewidth 2.0"<<endl;
  lead<<" @s2  line color 2"<<endl;
  lead<<" @s2  errorbar color 2"<<endl;
  lead<<"@type xydy"<<endl;
  jvec s=cA0P5/ZP5*Zv/spat_vol;
  //for(int t=0;t<=TH;t++) s[t]*=L[t];;
  //for(int t=0;t<=TH;t++) s[t]*=exp(aM*t);
  for(int t=0;t<TH+1;t++)
    {
      s[t]=L[t]*exp(-aM*t)/spat_vol*sqr(afpi_P5/a);
      if(t<TH) lead<<t<<" "<<s[t]<<"&"<<endl;
    }
  lead<<write_constant_with_error(A02,0,TH-1);
  
  cout<<"A02: "<<A02<<" mu mass: "<<aMu/a<<endl;
  
#if NPROJ == 4
  const int ipA0=0,ipV0=1,ipP5=2,ipS0=3;
#endif

#if NPROJ == 2
  const int ipA0=0,ipV0=1;
#endif

#if NPROJ == 1
  const int ipV0=0;
#endif
  
  const int iqVi_lVi=0,iqV0_lV0=1,iqAi_lAi=2,iqA0_lA0=3,iqVi_lAi=4,iqV0_lA0=5,iqAi_lVi=6,iqA0_lV0=7;
  const int RE=0,IM=1;

#if NPROJ >= 2
  jvec tot_qVi_lVi_pA0=load_improved(iqVi_lVi,ipA0,RE,r_odd,0,"/tmp/tmp_ViVi_pA0%s.xmg");
  jvec tot_qV0_lV0_pA0=load_improved(iqV0_lV0,ipA0,RE,r_odd); //
  jvec tot_qAi_lAi_pA0=load_improved(iqAi_lAi,ipA0,RE,r_evn);
  jvec tot_qA0_lA0_pA0=load_improved(iqA0_lA0,ipA0,RE,r_evn);
  jvec tot_qAi_lVi_pA0=load_improved(iqAi_lVi,ipA0,RE,r_unk);
  jvec tot_qA0_lV0_pA0=load_improved(iqA0_lV0,ipA0,RE,r_unk);
  jvec tot_qVi_lAi_pA0=load_improved(iqVi_lAi,ipA0,RE,r_unk); ////
  jvec tot_qV0_lA0_pA0=load_improved(iqV0_lA0,ipA0,RE,r_unk);
  jack qVi_lVi_pA0=constant_fit(tot_qVi_lVi_pA0,tmin,tmax,"qVi_lVi_pA0.xmg");
  jack qV0_lV0_pA0=constant_fit(tot_qV0_lV0_pA0,tmin,tmax,"qV0_lV0_pA0.xmg");
  jack qAi_lAi_pA0=constant_fit(tot_qAi_lAi_pA0,tmin,tmax,"qAi_lAi_pA0.xmg");
  jack qA0_lA0_pA0=constant_fit(tot_qA0_lA0_pA0,tmin,tmax,"qA0_lA0_pA0.xmg");
  jack qAi_lVi_pA0=constant_fit(tot_qAi_lVi_pA0,tmin,tmax,"qAi_lVi_pA0.xmg");
  jack qA0_lV0_pA0=constant_fit(tot_qA0_lV0_pA0,tmin,tmax,"qA0_lV0_pA0.xmg");
  jack qVi_lAi_pA0=constant_fit(tot_qVi_lAi_pA0,tmin,tmax,"qVi_lAi_pA0.xmg");
  jack qV0_lA0_pA0=constant_fit(tot_qV0_lA0_pA0,tmin,tmax,"qV0_lA0_pA0.xmg");
  
#if NPROJ == 4
  jvec tot_qVi_lVi_pP5=load_improved(iqVi_lVi,ipP5,RE,r_odd);
  jack qVi_lVi_pP5=constant_fit(tot_qVi_lVi_pP5,tmin,tmax,"qVi_lVi_pP5.xmg");
  // {
  //   ofstream fout("chris_test.xmg");
  //   fout<<"@type xydy"<<endl;
  //   fout<<tot_qVi_lVi_pP5*sinh(aMu)<<endl;
  //   //fout<<asinh(((tot_qVi_lVi_pP5/tot_qVi_lVi_pA0).subset(0,TH))*sinh(aMu))<<endl;
  //   //fout<<"&"<<endl;
  //   //fout<<effective_mass(cP5P5)/aMu<<endl;
  //   fout<<"&"<<endl;
  //   fout<<tot_qVi_lVi_pA0*sinh(aM)<<endl;
  //   //fout<<effective_mass(cP5P5)<<endl;
  //   //fout<<"&"<<endl;
  //   //fout<<write_constant_with_error(aM,tmin,TH)<<endl;
  //   //fout<<"&"<<endl;
  //   //fout<<write_constant_fit_plot(sinh(effective_mass(cP5P5))/sinh(aMu),constant_fit(sinh(effective_mass(cP5P5))/sinh(aMu),tmin,TH),tmin,TH)<<endl;
  // }
  cout<<smart_print(asinh(qVi_lVi_pP5/qVi_lVi_pA0*sinh(aMu)))<<endl;
  cout<<smart_print(aM)<<endl;
  cout<<smart_print(constant_fit(effective_mass(cP5P5),tmin,TH))<<endl;
#endif
#endif
    
  jvec tot_qVi_lVi_pV0=load_improved(iqVi_lVi,ipV0,IM,r_evn,phi,0,"/tmp/tmp_ViVi%s.xmg"); //null
  jvec tot_qVi_lVi_pV0_alt=load_improved(iqVi_lVi,ipV0,IM,r_evn,alt,0,"/tmp/tmp_ViVi%s.xmg"); //null
  jvec tot_qV0_lV0_pV0=load_improved(iqV0_lV0,ipV0,IM,r_evn,phi,0,"/tmp/tmp_V0V0%s.xmg"); //cut-off effects?
  jvec tot_qAi_lAi_pV0=load_improved(iqAi_lAi,ipV0,IM,r_evn,phi,0,"/tmp/tmp_AiAi%s.xmg"); //doubly odd?! c.e.?
  // {
  //   for(int icase=1;icase<=6;icase++)
  //     for(int int12=1;int12<=2;int12++)
  // 	{
  // 	  jvec tot_qAi_lAi_pV0_chris=load_improved_chris(iqAi_lAi,ipV0,IM,r_evn,icase,int12,"AiAi_chris_test_fixed_case_%d_int_t%d.xmg");
  // 	  if(icase==5 && int12==2)
  // 	    {
  // 	      jack A(njacks),B(njacks),CA(njacks),CB(njacks);
  // 	      case_5_int_t2_fit(CA,CB,A,B,tot_qAi_lAi_pV0_chris);
  // 	      cout<<A<<" "<<B<<" "<<A/B<<endl;
  // 	    }
  // 	}
  //   for(int i=1;i<=6;i++) load_improved_chris(iqAi_lAi,ipV0,IM,r_evn,i,3,"AiAi_chris_test_fixed_case_%d_int_t%d.xmg");
  //   for(int i=1;i<=6;i++) load_improved_chris(iqAi_lAi,ipV0,IM,r_evn,i,4,"AiAi_chris_test_fixed_case_%d_int_t%d.xmg");
  // }
  jvec tot_qA0_lA0_pV0=load_improved(iqA0_lA0,ipV0,IM,r_evn,phi,0,"/tmp/tmp_A0A0%s.xmg"); //singly odd... c.e?
  jvec tot_qAi_lVi_pV0=load_improved(iqAi_lVi,ipV0,RE,r_evn,phi,0,"/tmp/tmp_AiVi%s.xmg"); //ok
  // for(int icase=1;icase<=6;icase++)
  //   for(int int12=1;int12<=2;int12++)
  //     {
  // 	jvec tot_qAi_lVi_pV0_chris=load_improved_chris(iqAi_lVi,ipV0,RE,r_evn,icase,int12); //ok
  // 	//(tot_qAi_lVi_pV0_chris).print_to_file("AiVi_chris_test_fixed_case_%d_int_t%d.xmg",icase,int12);
  // 	//rectangle_integrate(tot_qAi_lVi_pV0_chris).print_to_file("AiVi_chris_test_fixed_int.xmg");
  //     }
  // load_improved_chris(iqAi_lVi,ipV0,RE,r_evn,5,3,"AiVi_chris_test_fixed_case_%d_int_t%d.xmg");
  // load_improved_chris(iqAi_lVi,ipV0,RE,r_evn,4,3,"AiVi_chris_test_fixed_case_%d_int_t%d.xmg");
  jvec tot_qA0_lV0_pV0=load_improved(iqA0_lV0,ipV0,RE,r_evn,phi,0,"/tmp/tmp_A0V0%s.xmg"); //ok
  jvec tot_qVi_lAi_pV0=load_improved(iqVi_lAi,ipV0,RE,r_odd,phi,0,"/tmp/tmp_ViAi%s.xmg"); //p.o
  //jvec tot_qVi_lAi_pV0_chris=load_improved_chris(iqVi_lAi,ipV0,RE,r_odd,0,"/tmp/tmp_ViAi%s_chris.xmg"); //ok
  //(tot_qVi_lAi_pV0_chris).print_to_file("ViAi_chris_test_fixed.xmg");
  //rectangle_integrate(tot_qVi_lAi_pV0_chris).print_to_file("ViAi_chris_test_fixed_int.xmg");
  jvec tot_qV0_lA0_pV0=load_improved(iqV0_lA0,ipV0,RE,r_odd,phi,0,"/tmp/tmp_V0A0%s.xmg"); //p.o
  jack qVi_lVi_pV0=constant_fit(tot_qVi_lVi_pV0,tmin,tmax,"qVi_lVi_pV0.xmg");
  jack qVi_lVi_pV0_alt=constant_fit(tot_qVi_lVi_pV0_alt,tmin,tmax,"qVi_lVi_pV0_alt.xmg");
  //jack qV0_lV0_pV0=constant_fit(tot_qV0_lV0_pV0,tmin,tmax,"qV0_lV0_pV0.xmg");
  jack qAi_lAi_pV0=constant_fit(tot_qAi_lAi_pV0,tmin,tmax,"qAi_lAi_pV0.xmg");
  //jack qA0_lA0_pV0=constant_fit(tot_qA0_lA0_pV0,tmin,tmax,"qA0_lA0_pV0.xmg");
  jack qAi_lVi_pV0=constant_fit(tot_qAi_lVi_pV0,tmin,tmax,"qAi_lVi_pV0.xmg");
  jack qA0_lV0_pV0=constant_fit(tot_qA0_lV0_pV0,tmin,tmax,"qA0_lV0_pV0.xmg");
  jack qA_lV_pV0=constant_fit(tot_qAi_lVi_pV0+tot_qA0_lV0_pV0,tmin,tmax,"qA_lV_pV0.xmg");
  jack qVi_lAi_pV0=constant_fit(tot_qVi_lAi_pV0,tmin,tmax,"qVi_lAi_pV0.xmg");
  jack qV0_lA0_pV0=constant_fit(tot_qV0_lA0_pV0,tmin,tmax,"qV0_lA0_pV0.xmg");
  jack qV_lA_pV0=constant_fit(tot_qVi_lAi_pV0+tot_qV0_lA0_pV0,tmin,tmax,"qV_lA_pV0.xmg");
  
  cout<<norm<<endl;
  
  debug_load=true;
  
  int ri=0;
  int iw=3;
  int rl=0;
  int orie=0;
  int r2=0;
  int phi_eta=0;
  int irev=0;
  int qins=0;
  int il=0;
  ofstream fout("test.xmg");
  fout<<"@type xydy"<<endl;
  int is=0;
  fout<<load_alone(icombo(ri,0,iw,rl,orie,r2,phi_eta,irev,qins,il))<<"@s"<<is++<<" legend \"bare\""<<endl;
  //fout<<load_alone(icombo(ri,iw,rl,!orie,r2,phi_eta,irev,qins,il))<<"@s"<<is++<<" legend \"orie\""<<endl;
  //fout<<load_alone(icombo(ri,iw,rl,orie,r2,!phi_eta,irev,qins,il))<<"@s"<<is++<<" legend \"phi_eta\""<<endl;
  //fout<<load_alone(icombo(ri,iw,rl,orie,r2,phi_eta,!irev,qins,il))<<"@s"<<is++<<" legend \"irev\""<<endl;
  //fout<<load_alone(icombo(ri,iw,rl,orie,r2,phi_eta,irev,!qins,il))<<"@s"<<is++<<" legend \"qins\""<<endl;
  //fout<<load_alone(icombo(ri,iw,rl,orie,!r2,phi_eta,irev,qins,il))<<"@s"<<is++<<" legend \"qr\""<<endl;
  //fout<<load_alone(icombo(ri,iw,!rl,orie,r2,phi_eta,irev,qins,il))<<"@s"<<is++<<" legend \"lr\""<<endl;
  fout<<load_alone(icombo(ri,0,iw,!rl,orie,!r2,phi_eta,irev,qins,il))<<"@s"<<is++<<" legend \"lr\""<<endl;
    
  cout<<e_m+en_int_mu-en_ext_mu<<" "<<e_m+en_ph<<endl;
  
  //Em + EM - Emu;
  //Em + K;
  
  double c1=-en_int_mu+en_ext_mu+en_ph;
  double c2=-e_m-en_int_mu+en_ext_mu;
  cout<<"pi int mu: "<<pi_int_mu<<endl;
  cout<<"pi ext mu: "<<pi_ext_mu<<endl;
  cout<<"en int mu: "<<en_int_mu<<endl;
  cout<<"en ext mu: "<<en_ext_mu<<endl;
  cout<<"en ph: "<<en_ph<<endl;
  cout<<c1<<" "<<c2<<endl;
  for(int t=14;t<T;t++)
    {
      
      //cout<<t<<" "<<-exp(c1*t)+exp(c2*t)<<endl; //E5 int
      cerr<<t<<" "<<-exp((en_ext_mu-en_int_mu+en_ph)*t)<<endl;
    }
  
  cout<<"E_nu+E_ext_mu: "<<en_nu+en_ext_mu<<", m_pi: "<<m_pi<<endl;
  
  return 0;
}
