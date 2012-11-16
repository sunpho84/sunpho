#include "include.h"

#include "../nf2/common_pars.cpp"

typedef char path_type[1024];

//0:A, 1:B, 2:AB 

int njacks;
int nmass;
int L,T,TH;
int ibeta;
int spat_vol;

int tmin,tmax;

const int PION=0,KAON_SINS=1,KAON_LINS=2;
const int LI=0,ST=1;

jvec load_loop_corr(int ik1,int ik2,int ik3,int ri)
{
  int icorr=ri+2*(ik1+nmass*(ik2+nmass*ik3));
  
  jvec a(T,njacks);
  
  a.load("data/oPPoLoop-sel_conf.1.dat",icorr);
  a/=spat_vol;
  
  return a;
}

jvec load_corr(const char *inter,const char *name,int ik1,int ik2,int tw,int ri)
{
  ik1=2*ik1+tw;
  ik2=2*ik2+tw;
  
  int icorr=ri+2*(ik1+nmass*ik2);
  
  jvec a(T,njacks);
      
  a.load(combine("data/o%so-%s_conf.1.dat",inter,name).c_str(),icorr);
  a/=-spat_vol;
  
  return a;
}

jvec load_corr_tma(const char *inter,const char *name,int ik1,int ik2,int ri)
{return 0.5*(load_corr(inter,name,ik1,ik2,0,ri)+load_corr(inter,name,ik1,ik2,0,ri));}

jvec load(const char *inter,const char *name,int what,int ri)
{
  jvec a;
  
  switch(what)
    {
    case PION:      a=load_corr_tma(inter,name,LI,LI,ri);break;
    case KAON_LINS: a=load_corr_tma(inter,name,ST,LI,ri);break;
    case KAON_SINS: a=load_corr_tma(inter,name,LI,ST,ri);break;
    }
  
  return a;
}

//fit stuff
int _fit_tmin,_fit_tmax,_fit_nrati;
double *_fit_corr,**_fit_rati;
double *_fit_dcorr,**_fit_drati;

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

int main(int narg,char **arg)
{
  init_latpars();
  
  FILE *fin=fopen("input","r");
  read_formatted_from_file_expecting((char*)&L,fin,"%d","L");
  T=2*L;
  TH=L;
  spat_vol=L*L*L;

  read_formatted_from_file_expecting((char*)&nmass,fin,"%d","nmass");
  read_formatted_from_file_expecting((char*)&njacks,fin,"%d","njacks");
  read_formatted_from_file_expecting((char*)&tmin,fin,"%d","tmin");
  read_formatted_from_file_expecting((char*)&tmax,fin,"%d","tmax");
  read_formatted_from_file_expecting((char*)&ibeta,fin,"%d","ibeta");
  
  ////////////// load mass original pion and kaon //////////////////
  
  //pion correlator
  jvec c_pion=load("PP","ss",PION,0).simmetrized(1);
  jvec c_pion_effmass=effective_mass(c_pion);
  c_pion.print_to_file("plots/pion.xmg");
  
  //kaon correlator
  jvec c_kaon=load("PP","ss",KAON_LINS,0).simmetrized(1);
  jvec c_kaon_effmass=effective_mass(c_kaon);
  
  
  ////////////// load mass ib correction to the kaon //////////////////
  
  //isospin breaking in the kaon
  jvec c_kaon_mib=load("PP","sd",KAON_LINS,0).simmetrized(1);
  jvec c_kaon_mib_slope=c_kaon_mib/c_kaon;
  
  
  //////////////// load self energy of the pion and compute slopes /////////////
  
  //self energy of the pion
  jvec c_pion_self=load("PP","es",PION,0).simmetrized(1);
  jvec c_pion_self_slope=c_pion_self/c_pion;
  c_pion_self_slope.print_to_file("plots/pion_self_slope.xmg");
  jack pion_self_slope=constant_fit(numerical_derivative(c_pion_self_slope),tmin,tmax,"plots/pion_self_slope_deriv.xmg");
  
  
  //////////////// load self energy of the kaon and compute slopes /////////////
  
  //self energy of the s for the kaon
  jvec c_kaon_s_self=load("PP","es",KAON_SINS,0).simmetrized(1);
  jvec c_kaon_s_self_slope=c_kaon_s_self/c_kaon;
  
  //self energy of the l for the kaon
  jvec c_kaon_l_self=load("PP","es",KAON_LINS,0).simmetrized(1);
  jvec c_kaon_l_self_slope=c_kaon_l_self/c_kaon;
  
  
  //////////////// load tadpole of the kaon and compute slopes /////////////
  
  //tadpole of the s for the kaon
  jvec c_kaon_s_tad=load("PP","ts",KAON_SINS,0).simmetrized(1);
  jvec c_kaon_s_tad_slope=c_kaon_s_tad/c_kaon;
  
  //self energy of the l for the kaon
  jvec c_kaon_l_tad=load("PP","ts",KAON_LINS,0).simmetrized(1);
  jvec c_kaon_l_tad_slope=c_kaon_l_tad/c_kaon;
  
  
  //////////////// load photon exchange for the pion and compute slopes /////////////
  
  //photon exchange with A on the strange for the pion
  jvec c_pion_exchange=load("PP","ee",PION,0).simmetrized(1);
  jvec c_pion_exchange_slope=c_pion_exchange/c_pion;
  c_pion_exchange_slope.print_to_file("plots/pion_exchange_slope.xmg");
  
  
  //////////////// load photon exchange for the kaon and compute slopes /////////////
  
  //photon exchange with A on the strange for the kaon
  jvec c_kaon_sA_exchange=load("PP","ee",KAON_LINS,0).simmetrized(1);
  jvec c_kaon_sA_exchange_slope=c_kaon_sA_exchange/c_kaon;
  c_kaon_sA_exchange_slope.print_to_file("plots/kaon_sA_exchange_slope.xmg");

  //photon exchange with A on the strange for the kaon
  jvec c_kaon_lA_exchange=load("PP","ee",KAON_SINS,0).simmetrized(1);
  jvec c_kaon_lA_exchange_slope=c_kaon_lA_exchange/c_kaon;
  c_kaon_lA_exchange_slope.print_to_file("plots/kaon_lA_exchange_slope.xmg");  
  
  //average kaon photon exchange with A inserted on l and s
  jvec c_kaon_exchange=0.5*(c_kaon_sA_exchange+c_kaon_lA_exchange);
  jvec c_kaon_exchange_slope=0.5*(c_kaon_sA_exchange_slope+c_kaon_lA_exchange_slope);
  cout<<endl;
  
  
  ///////////////////////// fit all slopes for kaon ///////////////////////////  
  
  jack kaon_Z2(njacks),kaon_M(njacks);
  jvec kaon_A(6,njacks),kaon_SL(6,njacks);
  jvec c_kaon_ratios[6]={c_kaon_mib_slope,c_kaon_l_self_slope,c_kaon_s_self_slope,c_kaon_exchange_slope,c_kaon_l_tad_slope,c_kaon_s_tad_slope};
  path_type kaon_path_ratios[6];
  sprintf(kaon_path_ratios[0],"plots/kaon_mib_slope.xmg");
  sprintf(kaon_path_ratios[1],"plots/kaon_l_self_slope.xmg");
  sprintf(kaon_path_ratios[2],"plots/kaon_s_self_slope.xmg");
  sprintf(kaon_path_ratios[3],"plots/kaon_exchange_slope.xmg");
  sprintf(kaon_path_ratios[4],"plots/kaon_l_tad_slope.xmg");
  sprintf(kaon_path_ratios[5],"plots/kaon_s_tad_slope.xmg");
  fit_mass_and_ratio(kaon_Z2,kaon_M,kaon_A,kaon_SL,c_kaon,c_kaon_ratios,6,tmin,tmax,"plots/kaon_effmass.xmg",kaon_path_ratios);
  cout<<"Kaon mass: "<<kaon_M<<endl;
  cout<<"Kaon slope ib: "<<kaon_SL[0]<<endl;
  cout<<"Kaon slope l self energy: "<<kaon_SL[1]<<endl;
  cout<<"Kaon slope s self energy: "<<kaon_SL[2]<<endl;
  cout<<"Kaon slope exchange: "<<kaon_SL[3]<<endl;
  cout<<"Kaon slope l tad: "<<kaon_SL[4]<<endl;
  cout<<"Kaon slope s tad: "<<kaon_SL[5]<<endl;
  cout<<endl;
  
  
  ///////////////////////// fit all slopes for pion ///////////////////////////  
  
  double e2=4*M_PI/137;
  
  jack pion_Z2(njacks),pion_M(njacks);
  jvec pion_A(2,njacks),pion_SL(2,njacks);
  jvec c_pion_ratios[2]={c_pion_self_slope,c_pion_exchange_slope};
  path_type pion_path_ratios[2];
  sprintf(pion_path_ratios[0],"plots/pion_self_slope.xmg");
  sprintf(pion_path_ratios[1],"plots/pion_exchange_slope.xmg");
  fit_mass_and_ratio(pion_Z2,pion_M,pion_A,pion_SL,c_pion,c_pion_ratios,2,tmin,tmax,"plots/pion_effmass.xmg",pion_path_ratios);
  cout<<"Pion mass: "<<pion_M<<endl;
  cout<<"Pion slope self energy: "<<pion_SL[0]<<endl;
  cout<<"Pion slope exchange: "<<pion_SL[1]<<endl;
  jack delta=-pion_SL[1]*2*pion_M*e2/2;
  cout<<"Squared pion mass difference: "<<delta<<endl;
  cout<<endl;
  
  delta.write_to_binfile("DeltaM2pi");
  (pion_M*pion_M).write_to_binfile("M2pi");
  
  ///////////////////////// load the disconnected //////////////////////////
  
  return 0;
}
