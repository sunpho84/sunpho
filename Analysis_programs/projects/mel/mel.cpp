#include "include.h"

int T=32,TH=T/2,njacks=15,nmass=3,nr=2;
int tmin_P5P5=12;
const int RE=0,IM=1;
const int LI=0,ST=1;

//fit stuff
int _fit_tmin,_fit_tmax,_fit_nrati;
double *_fit_corr,**_fit_rati;
double *_fit_dcorr,**_fit_drati;

const double e2=4*M_PI/137;

typedef char path_type[1024];

//fit function
double fun_corr_P5P5(double Z2_P5P5,double M_P5P5,double t)
{return Z2_P5P5*exp(-M_P5P5*TH)*cosh(M_P5P5*(TH-t))/sinh(M_P5P5);}
double fun_rati_P5P5(double A,double SL,double M,double t)
{return A+SL*(t-TH)*tanh(M*(t-TH));}

//calculate the chi square
double chi2_mass_ratio(double *p)
{
  double ch2=0;

  for(int t=_fit_tmin;t<=_fit_tmax;t++)
    {
      double ch2_corr=pow((_fit_corr[t]-fun_corr_P5P5(p[0],p[1],t))/_fit_dcorr[t],2);
      ch2+=ch2_corr;
      for(int irati=0;irati<_fit_nrati;irati++)
        {
          double ch2_rati=pow((_fit_rati[irati][t]-fun_rati_P5P5(p[2+2*irati],p[3+2*irati],p[1],t))/_fit_drati[irati][t],2);
          ch2+=ch2_rati;
        }
    }
  
  return ch2;
}

//wrapper for the calculation of the chi2
void ch2_mass_ratio_wr(int &npar,double *fuf,double &ch,double *p,int flag)
{ch=chi2_mass_ratio(p);}

void fit_mass_and_ratio(jack &C,jack &M,jvec &A,jvec &SL,jvec &corr,jvec *rati,int nrati,double tmin,double tmax,const char *path_plot_effcorr=NULL,path_type *path_plot_rati=NULL)
{
  tmin=min(TH,max(tmin,0));
  tmax=min(TH,max(tmax,0));
  
  _fit_tmin=tmin;
  _fit_tmax=tmax;
  _fit_nrati=nrati;

  _fit_corr=new double[TH+1];
  _fit_dcorr=new double[TH+1];
  _fit_rati=new double*[nrati];
  _fit_drati=new double*[nrati];
  for(int irati=0;irati<nrati;irati++)
    {
      _fit_rati[irati]=new double[TH+1];
      _fit_drati[irati]=new double[TH+1];
    }
  
  //first extimate
  two_pts_fit(M,C,corr,tmin,tmax);
  for(int irati=0;irati<nrati;irati++) linear_fit(A[irati],SL[irati],rati[irati],tmin,tmax);
  
  //calculate errors
  for(int t=tmin;t<=tmax;t++)
    {
      _fit_dcorr[t]=corr[t].err();
      for(int irati=0;irati<nrati;irati++) _fit_drati[irati][t]=rati[irati][t].err(); 
    }
  
  TMinuit minu;
  minu.SetPrintLevel(-1);
  minu.DefineParameter(0,"C",C[0],C.err(),C[0]-6*C.err(),C[0]+6*C.err());
  minu.DefineParameter(1,"M",M[0],M.err(),M[0]-6*M.err(),M[0]+6*M.err());
  for(int irati=0;irati<nrati;irati++)
    {
      minu.DefineParameter(2+2*irati,combine("A%02d",irati ).c_str(),A[irati][0],A[irati].err(),0,0);
      minu.DefineParameter(3+2*irati,combine("SL%02d",irati).c_str(),SL[irati][0],A[irati].err(),0,0);
    }
  
  minu.SetFCN(ch2_mass_ratio_wr);
  
  for(int ijack=0;ijack<njacks+1;ijack++)
    {
      //copy data so that glob function may access it
      for(int t=tmin;t<=tmax;t++)
        {
          _fit_corr[t]=corr[t][ijack];
          for(int irati=0;irati<nrati;irati++) _fit_rati[irati][t]=rati[irati][t][ijack];
        }
      
      minu.Migrad();
      
      //get back parameters
      double dum;
      minu.GetParameter(0,C.data[ijack],dum);
      minu.GetParameter(1,M.data[ijack],dum);
      for(int irati=0;irati<nrati;irati++)
        {
          minu.GetParameter(2+2*irati,A[irati].data[ijack],dum);
          minu.GetParameter(3+2*irati,SL[irati].data[ijack],dum);
        }
    }
  
  if(path_plot_effcorr!=NULL)
    {
      ofstream plot(path_plot_effcorr);
      
      plot<<"@type xydy"<<endl;
      plot<<"@s0 line type 0"<<endl;
      plot<<"@s0 symbol 1"<<endl;
      plot<<effective_mass(corr)<<endl;
      plot<<"&\n@type xy"<<endl;
      plot<<write_constant_with_error(M,tmin,tmax);
    }
  
  if(path_plot_rati!=NULL)
    for(int irati=0;irati<nrati;irati++)
      {
        int npoints=100;
        double x[npoints];
        jvec y(npoints,njacks);
        for(int ip=0;ip<npoints;ip++)
          {
            x[ip]=tmin+(tmax-tmin)/(npoints-1)*ip;
            for(int ijack=0;ijack<=njacks;ijack++) y[ip].data[ijack]=fun_rati_P5P5(A[irati][ijack],SL[irati][ijack],M[ijack],x[ip]);
          }
        
        ofstream plot(path_plot_rati[irati]);
        plot<<"@type xydy"<<endl;
        plot<<"@s0 line type 0"<<endl;
        plot<<"@s0 symbol 1"<<endl;
        plot<<rati[irati]<<endl;
        plot<<"&\n@type xy"<<endl;
        plot<<write_polygon(x,y);
      }
}

//insertion is done on the second passed mass
int icombo(int imrev,int imins,int r,int REIM)
{return REIM+2*(r+nr*(imins+nmass*imrev));}

jvec load(const char *name,int im1,int im2,int REIM,int rpar,int par)
{
  jvec a=jvec_load(combine("data/%s",name).c_str(),T,njacks,icombo(im1,im2,0,REIM)).simmetrized(par);
  if(rpar==0) return a;
  else return (a+rpar*jvec_load(combine("data/%s",name).c_str(),T,njacks,icombo(im1,im2,1,REIM)).simmetrized(par))/2;
}

int main()
{
  debug_load=false;
  debug_fit=false;
  
  //pure pion
  jvec P_corr_0=load("P5P5_00",LI,LI,RE,1,1);
  jack M_P=constant_fit(effective_mass(P_corr_0),12,23,"plots/pion_effmass.xmg");
  
  //pure kaon
  jvec K_corr_0=load("P5P5_00",LI,ST,RE,1,1);
  jack M_K=constant_fit(effective_mass(K_corr_0),12,23,"plots/kaon_effmass.xmg");
  
  //insertion of scalar density on the light quark
  jvec K_corr_l_S=load("P5P5_0S",ST,LI,RE,1,1);
  jvec K_ratio_l_S=(K_corr_l_S/K_corr_0);
  
  //pion exchange
  jvec P_corr_AB=load("P5P5_AB",LI,LI,RE,1,1);
  jvec P_ratio_AB=P_corr_AB/P_corr_0;
  
  //insertion of the pseudoscalar on the light and strange quark
  jvec K_corr_l_P=load("P5P5_0P",ST,LI,IM,-1,1);
  jvec K_corr_s_P=load("P5P5_0P",LI,ST,IM,-1,1);
  jvec K_ratio_l_P=K_corr_l_P/K_corr_0;
  jvec K_ratio_s_P=K_corr_s_P/K_corr_0;
  
  //insertion of the tadpole
  jvec K_corr_l_T=load("P5P5_0T",ST,LI,RE,1,1);
  jvec K_corr_s_T=load("P5P5_0T",LI,ST,RE,1,1);
  jvec K_ratio_l_T=-K_corr_l_T/K_corr_0;
  jvec K_ratio_s_T=-K_corr_s_T/K_corr_0;
  
  //self energy
  jvec K_corr_l_X=load("P5P5_0X",ST,LI,RE,1,1);
  jvec K_corr_s_X=load("P5P5_0X",LI,ST,RE,1,1);
  jvec K_ratio_l_X=-K_corr_l_X/K_corr_0;
  jvec K_ratio_s_X=-K_corr_s_X/K_corr_0;
  
  //kaon exchange
  jvec K_corr_AB=(load("P5P5_AB",ST,LI,RE,1,1)+load("P5P5_AB",LI,ST,RE,1,1))/2;
  jvec K_ratio_AB=-K_corr_AB/K_corr_0;
  
  ///////////////////////// fit all slopes for kaon ///////////////////////////
  
  jack K_Z2(njacks),K_M(njacks);
  jvec K_A(8,njacks),K_SL(8,njacks);
  jvec K_ratios[8]={K_ratio_l_S, K_ratio_l_X, K_ratio_s_X, K_ratio_AB, K_ratio_l_T, K_ratio_s_T, K_ratio_l_P, K_ratio_s_P};
  enum i_kaon_SL{   iKlS,        iKlX,        iKsX,        iKAB,       iKlT,        iKsT,        iKlP,        iKsP};
  
  path_type K_path_ratios[8]={"plots/K_l_mib_slope.xmg","plots/K_l_self_slope.xmg","plots/K_s_self_slope.xmg","plots/K_exchange_slope.xmg","plots/K_l_tad_slope.xmg","plots/K_s_tad_slope.xmg","plots/K_l_kib_slope.xmg","plots/K_s_kib_slope.xmg"};
  fit_mass_and_ratio(K_Z2,K_M,K_A,K_SL,K_corr_0,K_ratios,8,tmin_P5P5,TH,"plots/K_effmass.xmg",K_path_ratios);
  
  cout<<"Kaon mass: "<<smart_print(K_M)<<endl;
  cout<<"Kaon slope ib: "<<smart_print(K_SL[0])<<endl;
  cout<<"Kaon slope l self energy: "<<smart_print(K_SL[1])<<endl;
  cout<<"Kaon slope s self energy: "<<smart_print(K_SL[2])<<endl;
  cout<<"Kaon slope exchange: "<<smart_print(K_SL[3])<<endl;
  cout<<"Kaon slope l tad: "<<smart_print(K_SL[4])<<endl;
  cout<<"Kaon slope s tad: "<<smart_print(K_SL[5])<<endl;
  cout<<"Kaon slope l kib: "<<smart_print(K_SL[6])<<endl;
  cout<<"Kaon slope s kib: "<<smart_print(K_SL[7])<<endl;
  cout<<endl;
  
  ///////////////////////// fit all slopes for pion ///////////////////////////
  
  //pion exchange
  jvec P_corr_AB_PEC=(load("P5P5_AB_PEC",ST,LI,RE,1,1)+load("P5P5_AB_PEC",LI,ST,RE,1,1))/2;
  jvec P_ratio_AB_PEC=P_corr_AB_PEC/P_corr_0;
  
  jack P_Z2(njacks),P_M(njacks);
  jvec P_A(2,njacks),P_SL(2,njacks);
  jvec P_ratios[2]={P_ratio_AB,P_ratio_AB_PEC};
  path_type P_path_ratios[2]={"plots/P_exchange_slope","plots/P_exchange_slope_PEC"};
  fit_mass_and_ratio(P_Z2,P_M,P_A,P_SL,P_corr_0,P_ratios,2,tmin_P5P5,TH,"plots/P_effmass.xmg",P_path_ratios);
  cout<<"Pion mass: "<<smart_print(P_M)<<endl;
  cout<<"Pion slope exchange: "<<smart_print(P_SL[0])<<endl;
  cout<<"Pion slope exchange PECIONE: "<<smart_print(P_SL[1])<<endl;
  cout<<"PEC_UNN_DIFF: "<<smart_print((P_SL[1]-P_SL[0])/P_M/2)<<" "<<smart_print(1/sqr(P_M*24)/4)<<endl;
  
  //double a=0.086/0.197;
  //jack dM2P=-P_SL[0]*P_M*e2/2/a/a*1000000;
  //cout<<"Relative pion mass difference: "<<smart_print(dM2P)<<endl;
  
  return 0;
}
