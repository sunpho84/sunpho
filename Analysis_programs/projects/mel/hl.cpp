#include "../../src/include.h"

int T,TH,L;
int spat_vol;
double a,Zv,Za;

int nquarks;
double *amq;

int njacks;
int nr;

//parameters used for HVP
int nsubpoints=1;
double m_mu_phys=105.658e-3;

const int NPROJ=1;
const int nw=9;
const int norie=2; //number of orientation of the meson
const int nrev=2; //number of possible reversion
const int nins=3; //number of quarks insertion: 0(none), 1 or 2
const int RE=0,IM=1;

double eu=2.0/3,ed=-1.0/3;
double e2=4*M_PI/137.04;

int version=0;

//insertion is done on the second passed mass
int icombo_qonly(int imrev,int imins,int r,int REIM)
{return REIM+2*(r+nr*(imins+nquarks*imrev));}

jvec load_qonly(const char *name,int im1,int im2,int REIM,int rpar,int par,const char *outpath=NULL)
{
  jvec a=jvec_load(combine("data/corr_%s",name).c_str(),T,njacks,icombo_qonly(im1,im2,0,REIM));
  if(nr==1 || rpar==0) a=a.simmetrized(par);
  else
    {
      jvec b=jvec_load(combine("data/corr_%s",name).c_str(),T,njacks,icombo_qonly(im1,im2,1,REIM));
      a=(a+rpar*b).simmetrized(par)/2;
    }
  
  if(outpath) a.print_to_file(outpath);
  
  if(version==0 && REIM==1) a*=-1;
  
  return a;
}

struct dirstr_t
{
  char name[20];
  int rpar;
  int par;
  int ri;
  int glb_si;
  dirstr_t(const char *_name,int rpar,int par,int ri,int glb_si) : rpar(rpar),par(par),ri(ri),glb_si(glb_si) {memcpy(name,_name,20);}
private:
  dirstr_t(){}
};
dirstr_t P5P5("P5P5",+1,+1,RE,1);
dirstr_t V0P5("V0P5",-1,-1,IM,(version==1)?-1:1);
dirstr_t A0P5("A0P5",+1,-1,RE,1);

struct hl_t
{
  int imes;
  double aM_lep;
};

struct mes_t
{
  int iq1,iq2;
  double aM_mes;
  char name[20];
  int tmin,tmax;
};

struct mes_corrs_t
{
  //pure
  jvec c00;
  //self
  jvec c0M_ins1;
  jvec c0M_ins2;
  //tad
  jvec c0T_ins1;
  jvec c0T_ins2;
  //exchange
  jvec cLL;
  //scalar insertion
  jvec c0S_ins1;
  jvec c0S_ins2;
  //pseudo insertion
  jvec c0P_ins1;
  jvec c0P_ins2;
  mes_corrs_t(dirstr_t &d,mes_t &mes)
  {
    int s=d.glb_si;
    c00=s*load_qonly(combine("%s_00",d.name).c_str(),mes.iq1,mes.iq2,d.ri,d.rpar,d.par,combine("plots/%s/%s_00_corr.xmg",mes.name,d.name).c_str());
    c0M_ins2=s*load_qonly(combine("%s_0M",d.name).c_str(),mes.iq1,mes.iq2,d.ri,d.rpar,d.par,combine("plots/%s/%s_0M_ins2_corr.xmg",mes.name,d.name).c_str());
    c0M_ins1=s*load_qonly(combine("%s_0M",d.name).c_str(),mes.iq2,mes.iq1,d.ri,d.rpar,d.par,combine("plots/%s/%s_0M_ins1_corr.xmg",mes.name,d.name).c_str());
    c0T_ins2=s*load_qonly(combine("%s_0T",d.name).c_str(),mes.iq1,mes.iq2,d.ri,d.rpar,d.par,combine("plots/%s/%s_0T_ins2_corr.xmg",mes.name,d.name).c_str());
    c0T_ins1=s*load_qonly(combine("%s_0T",d.name).c_str(),mes.iq2,mes.iq1,d.ri,d.rpar,d.par,combine("plots/%s/%s_0T_ins1_corr.xmg",mes.name,d.name).c_str());
    cLL=s*load_qonly(combine("%s_LL",d.name).c_str(),mes.iq1,mes.iq2,d.ri,d.rpar,d.par,combine("plots/%s/%s_LL_corr.xmg",mes.name,d.name).c_str());
    c0S_ins2=s*load_qonly(combine("%s_0S",d.name).c_str(),mes.iq1,mes.iq2,d.ri,d.rpar,d.par,combine("plots/%s/%s_0S_ins2_corr.xmg",mes.name,d.name).c_str());
    c0S_ins1=s*load_qonly(combine("%s_0S",d.name).c_str(),mes.iq2,mes.iq1,d.ri,d.rpar,d.par,combine("plots/%s/%s_0S_ins1_corr.xmg",mes.name,d.name).c_str());
    //fix sign due to having forgotten "i"
    int pri=!d.ri;
    c0P_ins2=s*load_qonly(combine("%s_0P",d.name).c_str(),mes.iq1,mes.iq2,pri,-d.rpar,d.par,combine("plots/%s/%s_0P_ins2_corr.xmg",mes.name,d.name).c_str());
    c0P_ins1=s*load_qonly(combine("%s_0P",d.name).c_str(),mes.iq2,mes.iq1,pri,-d.rpar,d.par,combine("plots/%s/%s_0P_ins1_corr.xmg",mes.name,d.name).c_str());
    if(pri==IM) //change the sign when reading the imaginary part to fix "i"
      {
	c0P_ins2*=-1;
	c0P_ins1*=-1;
      }
  }
private:
  mes_corrs_t(){};
};

int nmes;
mes_t *mes,*deg_mes;
jack *dMcrit;
int nhl;
hl_t *hl;

////////////////////////////// ib on quarks only ///////////////////////

jvec load_spatav(const char *name,int im1,int im2,int REIM,int rpar,int par,const char *outpath=NULL)
{
  jvec a=(load_qonly(combine(name,0,0).c_str(),im1,im2,REIM,rpar,par)+
	  load_qonly(combine(name,1,1).c_str(),im1,im2,REIM,rpar,par)+
	  load_qonly(combine(name,2,2).c_str(),im1,im2,REIM,rpar,par))/3;
  if(outpath) a.print_to_file(outpath);
  return a;
}

///////////////////////////////////////////////////////////////////////////

int icombo(int ri,int iproj,int iw,int rl,int orie,int r2,int irev,int qins,int ilepton=0)
{return ri+2*(iproj+NPROJ*(iw+nw*(rl+nr*(orie+norie*(r2+nr*(irev+nrev*(qins+nins*ilepton)))))));}

jvec load_alone(int i)
{return jvec_load("data/corr_hl",T,njacks,i);}

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

const int unk[2]={+2,+0};
const int unk2[2]={0,+2};
const int evn[2]={+1,+1};
const int odd[2]={+1,-1};

jvec load_improved(int iw,int iproj,const int *orie_par,const int *rev_par,int qins,const char *path=NULL)
{
  jvec out(T,njacks);
  out=0;
  
  int n=0;
  for(int irev=0;irev<nrev;irev++)
    for(int r2=0;r2<nr;r2++)
      for(int orie=0;orie<norie;orie++)
	for(int rl=0;rl<nr;rl++)
	  for(int ri=0;ri<2;ri++)
	    {
	      int ic=icombo(ri,iproj,iw,rl,orie,r2,irev,qins,0);
	      jvec corr=load_alone(ic);
	      
	      //insertion on
	      if(ri==RE && r2==0) //nb keeping r2 and rl identical
		{
		  out+=orie_par[orie]*rev_par[irev]*corr;
		  n++;
		}
	      
	      if(path)
		{
		  ostringstream app;
		  app<<"qins_"<<qins<<"_irev_"<<irev<<"_r2_"<<r2<<"_orie_"<<orie<<"_rl_"<<rl<<"_ri_"<<ri;
		  ofstream fout(combine(path,app.str().c_str()).c_str());
		  fout<<"@type xydy"<<endl;
		  fout<<corr<<endl;
		  fout<<"@title \""<<ic<<"\""<<endl;
		}
	    }
  
  return out.simmetrized(1)/n;
}

void read_input(const char *path)
{
  FILE *fin=open_file(path,"r");
  read_formatted_from_file_expecting((char*)&T,fin,"%d","T");
  TH=L=T/2;
  spat_vol=L*L*L;
  read_formatted_from_file_expecting((char*)&njacks,fin,"%d","NJacks");
  read_formatted_from_file_expecting((char*)&a,fin,"%lg","a");
  a=a/0.19731;
  cout<<"Lattice spacing inverse: "<<1/a<<" GeV"<<endl;
  read_formatted_from_file_expecting((char*)&Zv,fin,"%lg","Zv");
  read_formatted_from_file_expecting((char*)&Za,fin,"%lg","Za");
  read_formatted_from_file_expecting((char*)&nquarks,fin,"%d","nquarks");
  amq=new double[nquarks];
  deg_mes=new mes_t[nquarks];
  for(int i=0;i<nquarks;i++)
    {
      read_formatted_from_file((char*)&(amq[i]),fin,"%lg","mq");
      deg_mes[i].iq1=deg_mes[i].iq2=i;
      char name[10];
      read_formatted_from_file((char*)&name,fin,"%s","name");
      sprintf(deg_mes[i].name,"%s-%s",name,name);
      read_formatted_from_file((char*)&(deg_mes[i].tmin),fin,"%d","tmin");
      read_formatted_from_file((char*)&(deg_mes[i].tmax),fin,"%d","tmax");
    }
  //
  read_formatted_from_file_expecting((char*)&nmes,fin,"%d","nmes");
  mes=new mes_t[nmes];
  for(int i=0;i<nmes;i++)
    {
      read_formatted_from_file((char*)&(mes[i].iq1),fin,"%d","iq1");
      read_formatted_from_file((char*)&(mes[i].iq2),fin,"%d","iq2");
      read_formatted_from_file((char*)&(mes[i].aM_mes),fin,"%lg","aM_mes");
      read_formatted_from_file((char*)&(mes[i].name),fin,"%s","name");
      read_formatted_from_file((char*)&(mes[i].tmin),fin,"%d","tmin");
      read_formatted_from_file((char*)&(mes[i].tmax),fin,"%d","tmax");
    }
  //
  read_formatted_from_file_expecting((char*)&nhl,fin,"%d","nhl");
  hl=new hl_t[nhl];
  for(int i=0;i<nhl;i++)
    {
      read_formatted_from_file((char*)&(hl[i].imes),fin,"%d","imes");
      read_formatted_from_file((char*)&(hl[i].aM_lep),fin,"%lg","aM_lep");
    }
  
  read_formatted_from_file_expecting((char*)&nr,fin,"%d","NR");
}

namespace ratfit
{
  int tmin,tmax,nrati;
  double *corr,**rati;
  double *dcorr,**drati;
  
  //fit function
  double fun_corr_P5P5(double Z2,double M,double t)
  {return Z2*exp(-M*TH)*cosh(M*(TH-t))/sinh(M);}
  double fun_rati_P5P5(double A,double SL,double M,double t)
  {return A+SL*(t-TH)*tanh(M*(t-TH));}
  
  //calculate the chi square
  double chi2_mass_ratio(double *p)
  {
    double ch2=0;
    
    for(int t=tmin;t<=tmax;t++)
      {
	double ch2_corr=pow((corr[t]-fun_corr_P5P5(p[0],p[1],t))/dcorr[t],2);
	ch2+=ch2_corr;
	for(int irati=0;irati<nrati;irati++)
	  {
	    double ch2_rati=pow((rati[irati][t]-fun_rati_P5P5(p[2+2*irati],p[3+2*irati],p[1],t))/drati[irati][t],2);
	    ch2+=ch2_rati;
	  }
      }
    
    return ch2;
  }
  
  //wrapper for the calculation of the chi2
  void ch2_mass_ratio_wr(int &npar,double *fuf,double &ch,double *p,int flag)
  {ch=chi2_mass_ratio(p);}
}

void fit_mass_and_ratio(jack &C,jack &M,jvec &A,jvec &SL,jvec &corr,vector<jvec> &rati,double tmin,double tmax,const char *path_plot_effcorr,vector<string> &path_plot_rati)
{
  int nrati=rati.size();
  
  tmin=min(TH,max(tmin,0));
  tmax=min(TH,max(tmax,0));
  
  ratfit::tmin=tmin;
  ratfit::tmax=tmax;
  ratfit::nrati=nrati;
  
  ratfit::corr=new double[TH+1];
  ratfit::dcorr=new double[TH+1];
  ratfit::rati=new double*[nrati];
  ratfit::drati=new double*[nrati];
  for(int irati=0;irati<nrati;irati++)
    {
      ratfit::rati[irati]=new double[TH+1];
      ratfit::drati[irati]=new double[TH+1];
    }
  
  //first extimate
  two_pts_fit(M,C,corr,tmin,tmax);
  for(int irati=0;irati<nrati;irati++) linear_fit(A[irati],SL[irati],rati[irati],tmin,tmax);
  
  //calculate errors
  for(int t=tmin;t<=tmax;t++)
    {
      ratfit::dcorr[t]=corr[t].err();
      for(int irati=0;irati<nrati;irati++) ratfit::drati[irati][t]=rati[irati][t].err();
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
  
  minu.SetFCN(ratfit::ch2_mass_ratio_wr);
  
  for(int ijack=0;ijack<njacks+1;ijack++)
    {
      //copy data so that glob function may access it
      for(int t=tmin;t<=tmax;t++)
        {
	  ratfit::corr[t]=corr[t][ijack];
          for(int irati=0;irati<nrati;irati++) ratfit::rati[irati][t]=rati[irati][t][ijack];
        }
      
      minu.Migrad();
      
      //get back parameters
      double dum;
      minu.GetParameter(0,C.data[ijack],dum);
      minu.GetParameter(1,M.data[ijack],dum);
      for(int irati=0;irati<nrati;irati++)
	{
	  minu.GetParameter(2+2*irati,A[irati][ijack],dum);
	  minu.GetParameter(3+2*irati,SL[irati][ijack],dum);
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
  
  for(int irati=0;irati<nrati;irati++)
    {
      int npoints=100;
      double x[npoints];
      jvec y(npoints,njacks);
      for(int ip=0;ip<npoints;ip++)
	{
	  x[ip]=tmin+(tmax-tmin)/(npoints-1)*ip;
	  for(int ijack=0;ijack<=njacks;ijack++) y[ip].data[ijack]=ratfit::fun_rati_P5P5(A[irati][ijack],SL[irati][ijack],M[ijack],x[ip]);
	}
      
      ofstream plot(path_plot_rati[irati]);
      plot<<"@type xydy"<<endl;
      plot<<"@s0 line type 0"<<endl;
      plot<<"@s0 symbol 1"<<endl;
      plot<<rati[irati]<<endl;
      plot<<"&\n@type xy"<<endl;
      plot<<write_polygon(x,y);
    }
  
  SL*=-1;
}

void fit_mass_and_ratio(jack &C,jack &M,jack &A,jack &SL,jvec &corr,jvec &rati,double tmin,double tmax,const char *path_plot_effcorr,const char *path_plot_rati)
{
  vector<string> path_plot_rati_vect(1,path_plot_rati);
  jvec Av(1,njacks),SLv(1,njacks);
  vector<jvec> rativ(1,rati);
  fit_mass_and_ratio(C,M,Av,SLv,corr,rativ,tmin,tmax,path_plot_effcorr,path_plot_rati_vect);
  A=Av[0];
  SL=SLv[0];
}

double HVP_kernel(double Q2)
{
  if(Q2==0) return 0;
  
  double Z=-(1-sqrt(1+4*sqr(m_mu_phys)/Q2))/(2*sqr(m_mu_phys));
  return sqr(m_mu_phys)*Q2*pow(Z,3)*(1-Q2*Z)/(1+sqr(m_mu_phys*Z)*Q2);
}

jvec calc_HVP(jvec corr)
{
  jvec HVP(T*nsubpoints,njacks);
  
  //point at q=0
  HVP[0]=0;
  
  //all other q
  for(int i=1;i<T*nsubpoints;i++)
    {
      double q0=2*M_PI*i/(T*nsubpoints);
      HVP[i]=0;
      for(int t=0;t<=T/2;t++) HVP[i]+=2*((cos(q0*t)-1)/sqr(q0)+sqr(t)/2)*corr[t];
    }
  return HVP*sqr(Za)/sqr(137.0*M_PI);
}

void write_HVP(const char *path,jvec corr,double val_0=0)
{
  corr[0]=val_0;
  ofstream out(path);
  out<<"@type xydy"<<endl;
  for(int i=0;i<T*nsubpoints;i++) out<<sqr(2*M_PI*i/(T*nsubpoints)/a)<<" "<<corr[i]<<endl;
}

int main(int narg,char **arg)
{
  debug_load=false;
  debug_fit=false;
  
  //read input
  string str("input");
  if(narg>2) str=arg[1];
  read_input(str.c_str());
  
  const int use_tad=1;
  cout<<"Use tad: "<<use_tad<<endl;
  
  /////////////////////////////////////// determine critical mass corrections ///////////////////////////////
  
  //eq.76 of 1303.4896
  dMcrit=new jack[nquarks];
  for(int iquark=0;iquark<nquarks;iquark++)
    {
      //load the correlation function
      mes_corrs_t cV0P5(V0P5,deg_mes[iquark]);
      
      //take the numerical derivative in time, so to check kappa is the critical one
      numerical_derivative(cV0P5.c00).print_to_file(combine("plots/%s/V0P5_00_tder.xmg",deg_mes[iquark].name).c_str());
      
      //now compute the correction
      jvec num=numerical_derivative(/*cV0P5.cLL+*/2*(cV0P5.c0M_ins2+use_tad*cV0P5.c0T_ins2)); //LL is zero due to charge symmetry
      jvec den=numerical_derivative(cV0P5.c0P_ins2);
      constant_fit(num,deg_mes[iquark].tmin,deg_mes[iquark].tmax,combine("plots/%s/num_dMcrit_corr.xmg",deg_mes[iquark].name).c_str());
      constant_fit(den,deg_mes[iquark].tmin,deg_mes[iquark].tmax,combine("plots/%s/den_dMcrit_corr.xmg",deg_mes[iquark].name).c_str());
      dMcrit[iquark]=constant_fit(num/den,deg_mes[iquark].tmin,deg_mes[iquark].tmax,combine("plots/%s/dMcrit_corr.xmg",deg_mes[iquark].name).c_str());
      dMcrit[iquark]/=-2; //including factor -1/2
      
      cout<<"dMcrit["<<iquark<<"]: "<<smart_print(dMcrit[iquark])<<endl;
    }
  
  cout<<"-------------------------------------------------------------"<<endl;
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  double mu_ms=2;//GeV
  double Zud=32*sqr(M_PI)/(22.596-6*log(a*mu_ms)); //sign flipped to collect e^2(eu^2-ed^2)
  cout<<"Zud: "<<Zud<<endl;
  
  //fit masses (then, also decay constants)
  jvec aM(nmes,njacks);
  for(int imes=0;imes<nmes;imes++)
    {
      mes_corrs_t cP5P5(P5P5,mes[imes]);
      mes_corrs_t cA0P5(A0P5,mes[imes]);
      
      aM[imes]=constant_fit(effective_mass(cP5P5.c00),mes[imes].tmin,mes[imes].tmax,combine("plots/%s/P5P5_00_meff.xmg",mes[imes].name).c_str());
      
      jvec P5P5_LL_ratio=cP5P5.cLL/cP5P5.c00;
      P5P5_LL_ratio.print_to_file(combine("plots/%s/P5P5_LL_ratio.xmg",mes[imes].name));
      
      jvec P5P5_0M_ratio=cP5P5.c0M_ins1/cP5P5.c00;
      P5P5_0M_ratio.print_to_file(combine("plots/%s/P5P5_0M_ratio.xmg",mes[imes].name));
      
      cout<<"Mes "<<imes<<" mismatch: "<<smart_print(aM[imes])<<"/"<<mes[imes].aM_mes<<"-1 = "<<smart_print(aM[imes]/mes[imes].aM_mes-1)<<endl;
    }
  
  //coefficients for fit
  jack C(njacks),M(njacks),A(njacks);
  
  //define mesons
  const int iPi=0;
  const int iK=1;
  const int iD=2;
  
  //corrections for pion
  mes_corrs_t Pi_P5P5(P5P5,mes[iPi]);
  jvec Pi_P5P5_LL_ratio=Pi_P5P5.cLL/Pi_P5P5.c00;
  jack SL_pion;
  fit_mass_and_ratio(C,M,A,SL_pion,Pi_P5P5.c00,Pi_P5P5_LL_ratio,mes[iPi].tmin,mes[iPi].tmax,"/dev/null","plots/Pi/P5P5_LL_slope.xmg");
  jack dM2Pi_de2=sqr(eu-ed)*SL_pion*2*aM[iPi]/2;
  cout<<"Pion exchange slope: "<<smart_print(SL_pion)<<endl;
  cout<<"Pion mass (latt units): "<<smart_print(aM[iPi])<<endl;
  cout<<"daM2Pi/de^2: "<<smart_print(dM2Pi_de2)<<endl;
  cout<<"dM2Pi (phys=0.001261 GeV^2): "<<smart_print(dM2Pi_de2*e2/sqr(a))<<" GeV^2"<<endl;
  
  cout<<"-----------"<<endl;
  
  //corrections for K
  mes_corrs_t K_P5P5(P5P5,mes[iK]);
  jvec K_P5P5_QED_ratio=(-2*amq[0]/Zud*K_P5P5.c0S_ins2+(K_P5P5.cLL-(K_P5P5.c0M_ins2+use_tad*K_P5P5.c0T_ins2)-dMcrit[0]*K_P5P5.c0P_ins2))/K_P5P5.c00;
  jvec K_P5P5_QCD_ratio=K_P5P5.c0S_ins2/K_P5P5.c00;
  vector<jvec> K_P5P5_ratios(2);
  K_P5P5_ratios[0]=K_P5P5_QED_ratio;
  K_P5P5_ratios[1]=K_P5P5_QCD_ratio;
  vector<string> K_P5P5_ratios_tag(2);
  K_P5P5_ratios_tag[0]="plots/K/P5P5_QED_slope.xmg";
  K_P5P5_ratios_tag[1]="plots/K/P5P5_QCD_slope.xmg";
  jvec A_K(2,njacks),SL_K(2,njacks);
  fit_mass_and_ratio(C,M,A_K,SL_K,K_P5P5.c00,K_P5P5_ratios,mes[iK].tmin,mes[iK].tmax,"/dev/null",K_P5P5_ratios_tag);
  jack SL_QED_K=SL_K[0];
  jack SL_QCD_K=SL_K[1];
  jack dM2K_QED_de2=SL_QED_K*2*aM[iK]*(sqr(eu)-sqr(ed));
  cout<<"daM2K_QED/de^2: "<<smart_print(dM2K_QED_de2)<<endl;
  cout<<"dM2K_QED: "<<smart_print(dM2K_QED_de2*e2/sqr(a))<<" GeV^2"<<endl;
  
  cout<<"Epsilon: "<<smart_print(dM2K_QED_de2/dM2Pi_de2-1)<<endl;
  
  cout<<"-----------"<<endl;
  
  //corrections for D
  mes_corrs_t D_P5P5(P5P5,mes[iD]);
  jvec D_P5P5_QED_ratio=(-2*amq[0]/Zud*D_P5P5.c0S_ins2+(D_P5P5.cLL-(D_P5P5.c0M_ins2+use_tad*D_P5P5.c0T_ins2)-dMcrit[0]*D_P5P5.c0P_ins2))/D_P5P5.c00;
  vector<jvec> D_P5P5_ratios(1);
  D_P5P5_ratios[0]=D_P5P5_QED_ratio;
  vector<string> D_P5P5_ratios_tag(1);
  D_P5P5_ratios_tag[0]="plots/D/P5P5_QED_slope.xmg";
  jvec A_D(1,njacks),SL_D(1,njacks);
  fit_mass_and_ratio(C,M,A_D,SL_D,D_P5P5.c00,D_P5P5_ratios,mes[iD].tmin,mes[iD].tmax,"/dev/null",D_P5P5_ratios_tag);
  jack SL_QED_D=SL_D[0];
  cout<<"QED slope: "<<smart_print(SL_QED_D)<<endl;
  jack dMD_QED_de2=SL_QED_D*(sqr(ed)-sqr(eu));
  cout<<"MD: "<<smart_print(M/a)<<" GeV or "<<smart_print(aM[iD]/a)<<" GeV"<<endl;
  cout<<"daMD_QED/de^2: "<<smart_print(dMD_QED_de2)<<endl;
  cout<<"dMD_QED: "<<smart_print(dMD_QED_de2*e2/a*1000)<<" MeV"<<endl;

  
  //(numerical_derivative(temp)*aM[iK]/(numerical_derivative(Pi_P5P5_LL_ratio)*aM[iPi])-1).print_to_file("/tmp/test.xmg");
  //jack K_CORR=(eu*eu-ed*ed)*e2*(K_EM_CORR-K_KA_CORR)*(2*kaon_M);
  //cout<<"Corr: "<<K_CORR/a/a<<" MeV^2"<<endl;
  //cout<<"EM: "<<K_EM_CORR<<" KA: "<<K_KA_CORR<<endl;
  
  /*
  //////////////////////////
  
  //pure kaon
  jvec Ka_P5P5_00_corr=load("P5P5_00",LI,ST,RE,1,1,"plots/Ka_P5P5_00_corr.xmg");
  
  //kaon self
  jvec Ka_P5P5_0M_sins_corr=load("P5P5_0M",LI,ST,RE,1,1,"plots/Ka_P5P5_0M_sins_corr.xmg");
  jvec Ka_P5P5_0M_lins_corr=load("P5P5_0M",ST,LI,RE,1,1,"plots/Ka_P5P5_0M_lins_corr.xmg");
  jvec Ka_P5P5_0M_sins_ratio=Ka_P5P5_0M_sins_corr/Ka_P5P5_00_corr;
  jvec Ka_P5P5_0M_lins_ratio=Ka_P5P5_0M_lins_corr/Ka_P5P5_00_corr;
  Ka_P5P5_0M_sins_ratio.print_to_file("plots/Ka_P5P5_0M_sins_ratio.xmg");
  Ka_P5P5_0M_lins_ratio.print_to_file("plots/Ka_P5P5_0M_lins_ratio.xmg");
  //kaon tad
  jvec Ka_P5P5_0T_sins_corr=load("P5P5_0T",LI,ST,RE,1,1,"plots/Ka_P5P5_0T_sins_corr.xmg");
  jvec Ka_P5P5_0T_lins_corr=load("P5P5_0T",ST,LI,RE,1,1,"plots/Ka_P5P5_0T_lins_corr.xmg");
  jvec Ka_P5P5_0T_sins_ratio=Ka_P5P5_0T_sins_corr/Ka_P5P5_00_corr;
  jvec Ka_P5P5_0T_lins_ratio=Ka_P5P5_0T_lins_corr/Ka_P5P5_00_corr;
  Ka_P5P5_0T_sins_ratio.print_to_file("plots/Ka_P5P5_0T_sins_ratio.xmg");
  Ka_P5P5_0T_lins_ratio.print_to_file("plots/Ka_P5P5_0T_lins_ratio.xmg");
  //kaon exchange
  jvec Ka_P5P5_LL_corr=load("P5P5_LL",LI,ST,RE,1,1,"plots/Ka_P5P5_LL_corr.xmg");
  jvec Ka_P5P5_LL_ratio=Ka_P5P5_LL_corr/Ka_P5P5_00_corr;
  Ka_P5P5_LL_ratio.print_to_file("plots/Ka_P5P5_LL_ratio.xmg");
  //kaon scalar insertion
  jvec Ka_P5P5_0S_lins_corr=load("P5P5_0S",ST,LI,RE,1,1,"plots/Ka_P5P5_0S_lins_corr.xmg");
  jvec Ka_P5P5_0S_sins_corr=load("P5P5_0S",LI,ST,RE,1,1,"plots/Ka_P5P5_0S_sins_corr.xmg");
  jvec Ka_P5P5_0S_lins_ratio=Ka_P5P5_0S_lins_corr/Ka_P5P5_00_corr;
  jvec Ka_P5P5_0S_sins_ratio=Ka_P5P5_0S_sins_corr/Ka_P5P5_00_corr;
  Ka_P5P5_0S_lins_ratio.print_to_file("plots/Ka_P5P5_0S_lins_ratio.xmg");
  Ka_P5P5_0S_sins_ratio.print_to_file("plots/Ka_P5P5_0S_sins_ratio.xmg");
  //kaon pseudo insertion
  jvec Ka_P5P5_0P_lins_corr=load("P5P5_0P",ST,LI,IM,-1,1,"plots/Ka_P5P5_0P_lins_corr.xmg");
  jvec Ka_P5P5_0P_sins_corr=load("P5P5_0P",LI,ST,IM,-1,1,"plots/Ka_P5P5_0P_sins_corr.xmg");
  jvec Ka_P5P5_0P_lins_ratio=Ka_P5P5_0P_lins_corr/Ka_P5P5_00_corr;
  jvec Ka_P5P5_0P_sins_ratio=Ka_P5P5_0P_sins_corr/Ka_P5P5_00_corr;
  Ka_P5P5_0P_lins_ratio.print_to_file("plots/Ka_P5P5_0P_lins_ratio.xmg");
  Ka_P5P5_0P_sins_ratio.print_to_file("plots/Ka_P5P5_0P_sins_ratio.xmg");
  
  jack aM_Pi_A0(njacks),ZA0_Pi(njacks),ZP5_Pi(njacks);
  two_pts_SL_fit(aM_Pi_A0,ZA0_Pi,ZP5_Pi,Pi_A0P5_00_corr,Pi_P5P5_00_corr,tmin,TH-1,tmin,TH,"plots/Pi_P5P5_00_effmass.xmg",NULL,"/tmp/list.txt",-1,+1);
  jack aM_Pi_P5(njacks);
  jack Z2P5_Pi(njacks);
  two_pts_fit(aM_Pi_P5,Z2P5_Pi,Pi_P5P5_00_corr,tmin,TH);
  jack af_Pi_P5=2*amq[LI]*ZP5_Pi/sqr(aM_Pi_P5);
  jack af_Pi_A0=-ZA0_Pi/aM_Pi_A0*Zv;
  cout<<" ZP5"<<endl;
  cout<<"  from A0P5: "<<smart_print(sqrt(Z2P5_Pi))<<endl;
  cout<<"  from P5P5: "<<smart_print(ZP5_Pi)<<endl;
  cout<<" aM_Pi"<<endl;
  cout<<"  from P5P5: "<<smart_print(aM_Pi_P5)<<endl;
  cout<<"  from A0P5: "<<smart_print(aM_Pi_A0)<<endl;
  cout<<" fpi"<<endl;
  cout<<"  from P5: "<<smart_print(af_Pi_P5/a)<<endl;
  cout<<"  from A0: "<<smart_print(af_Pi_A0/a)<<endl;
  
  //compute normalization
  jack hl_norm=2*aM_Pi_P5/ZP5_Pi*(af_Pi_P5);
  
  const int ipV0=0;
  
  const int iqVi_lVi=0,iqV0_lV0=1; //all null
  const int iqAi_lVi=6,iqA0_lV0=7;
  //const int iqAi_lAi=2,iqA0_lA0=3,iqVi_lAi=4,iqV0_lA0=5; //redundant because we put 1-g5 in the lepton side
  //---------                                        ori  q  rev
  //jvec tot_qAi_lAi_pV0=load_improved(iqAi_lAi,ipV0,evn,odd,1,"plots/tmp_AiAi%s.xmg"); //null
  //jvec tot_qA0_lA0_pV0=load_improved(iqA0_lA0,ipV0,evn,odd,1,"plots/tmp_A0A0%s.xmg"); //null
  jvec tot_qVi_lVi_pV0=load_improved(iqVi_lVi,ipV0,unk,unk,2,"plots/tmp_ViVi%s.xmg"); //tiny but not null!!!!
  jvec tot_qV0_lV0_pV0=load_improved(iqV0_lV0,ipV0,unk,unk,2,"plots/tmp_V0V0%s.xmg"); //tiny/compatible with 0
  jvec tot_qAi_lVi_pV0=load_improved(iqAi_lVi,ipV0,evn,odd,2,"plots/tmp_AiVi%s.xmg"); //ok
  jvec tot_qA0_lV0_pV0=load_improved(iqA0_lV0,ipV0,evn,odd,2,"plots/tmp_A0V0%s.xmg"); //ok
  //jvec tot_qVi_lAi_pV0=load_improved(iqVi_lAi,ipV0,evn,evn,1,"plots/tmp_ViAi%s.xmg"); //ok
  //jvec tot_qV0_lA0_pV0=load_improved(iqV0_lA0,ipV0,unk,unk,1,"plots/tmp_V0A0%s.xmg"); //null?!
  
  //jack qAi_lAi_pV0=constant_fit(tot_qAi_lAi_pV0,tmin,tmax,"plots/qAi_lAi_pV0.xmg");
  //jack qA0_lA0_pV0=constant_fit(tot_qA0_lA0_pV0,tmin,tmax,"plots/qA0_lA0_pV0.xmg");
  //jack qVi_lVi_pV0=constant_fit(tot_qVi_lVi_pV0,tmin,tmax,"plots/qVi_lVi_pV0.xmg");
  //jack qV0_lV0_pV0=constant_fit(tot_qV0_lV0_pV0,tmin,tmax,"plots/qV0_lV0_pV0.xmg");
  jack qAi_lVi_pV0=constant_fit(tot_qAi_lVi_pV0,tmin,tmax,"plots/qAi_lVi_pV0.xmg");
  jack qA0_lV0_pV0=constant_fit(tot_qA0_lV0_pV0,tmin,tmax,"plots/qA0_lV0_pV0.xmg");
  //jack qVi_lAi_pV0=constant_fit(tot_qVi_lAi_pV0,tmin,tmax,"plots/qVi_lAi_pV0.xmg");
  //jack qV0_lA0_pV0=constant_fit(tot_qV0_lA0_pV0,tmin,tmax,"plots/qV0_lA0_pV0.xmg");
  //jack qA_lV_pV0=constant_fit(tot_qAi_lVi_pV0+tot_qA0_lV0_pV0,tmin,tmax,"plots/qA_lV_pV0.xmg");
  //jack qV_lA_pV0=constant_fit(tot_qVi_lAi_pV0+tot_qV0_lA0_pV0,tmin,tmax,"plots/qV_lA_pV0.xmg");
  
  jvec TL_tot_qA0_lV0_pV0=load_improved(iqA0_lV0,ipV0,evn,evn,0,"plots/tmp_TL_A0V0%s.xmg"); //ok
  jack TL_qA0_lV0_pV0=constant_fit(TL_tot_qA0_lV0_pV0,tmin,tmax,"plots/TL_qA0_lV0_pV0.xmg");
  
  jvec test=Zv*TL_tot_qA0_lV0_pV0*spat_vol/(2*aM_Pi_P5*sqr(hl[0].aM_lep)*(1-sqr(hl[0].aM_lep/aM_Pi_P5))*ZP5_Pi/(2*aM_Pi_P5));
  test.print_to_file("/tmp/test.xmg");
  cout<<"afpi: "<<af_Pi_A0<<"       "<<test[10]<<endl;
  
  // Gamma_0_alpha/Gamma_0_tree [e+f, gamma5 gamma0xgamma0(1-gamma5)]
  auto printer=[](jvec c,const char *path,const char *ghadr,const char *glept)->void
    {
      ofstream fout(path);
      fout<<"@type xydy"<<endl;
      c[0]=c[1];
      fout<<c.subset(0,TH)<<endl;
      fout<<"@title \"\\xG\\0\\s0\\N\\S\\xa\\0\\N/\\xG\\0\\s0\\N\\Stree\\N [e+f, ("<<ghadr<<")\\shadr\\N x \\xg\\0\\s"<<glept<<"\\N(1-\\xg\\0\\s5\\N)\\slept\\N]\""<<endl;
      fout<<
      "@world 0.498788927, -0.0437564264, 23.6620675, -0.0235958891\n"
      "@xaxis  label \"t\"\n"
      "@yaxis  label \"ratio\"\n"
      "@xaxis  label char size 1.500000\n"
      "@yaxis  label char size 1.500000\n"
      "@yaxis  tick major 0.005\n"
      "@s0 symbol 1\n"
      "@s0 symbol size 0.510000\n"
      "@s0 symbol color 2\n"
      "@s0 symbol fill color 2\n"
      "@s0 symbol fill pattern 1\n"
      "@s0 symbol linewidth 2.0\n"
      "@s0 line type 0\n"
      "@s0 line color 2\n"
      "@    s0 errorbar color 2\n"
      "@    s0 errorbar linewidth 2.0\n"
      "@    s0 errorbar riser linewidth 2.0"<<endl;
    };
  printer(tot_qA0_lV0_pV0/TL_tot_qA0_lV0_pV0,"plots/e_lus_f_A0_over_TL.xmg","\\xg\\0\\s5\\N \\xg\\0\\s0\\N","0");
  printer(tot_qAi_lVi_pV0/TL_tot_qA0_lV0_pV0,"plots/e_lus_f_Ai_over_TL.xmg","\\xg\\0\\s5\\N \\xg\\0\\si\\N","i");
  printer(tot_qV0_lV0_pV0/TL_tot_qA0_lV0_pV0,"plots/e_lus_f_V0_over_TL.xmg","\\xg\\0\\s0\\N","0");
  printer(tot_qVi_lVi_pV0/TL_tot_qA0_lV0_pV0,"plots/e_lus_f_Vi_over_TL.xmg","\\xg\\0\\si\\N","i");
  //constant_fit((tot_qA0_lV0_pV0+tot_qAi_lVi_pV0)/TL_tot_qA0_lV0_pV0,tmin,tmax,"e_lus_f_A_TL.xmg");
  
  int ri=0;
  int iw=3;
  int rl=0;
  int orie=0;
  int r2=0;
  int irev=0;
  int qins=0;
  int il=0;
  ofstream fout("plots/test.xmg");
  fout<<"@type xydy"<<endl;
  int is=0;
  fout<<load_alone(icombo(ri,0,iw,rl,orie,r2,irev,qins,il))<<"@s"<<is++<<" legend \"bare\""<<endl;
  fout<<load_alone(icombo(ri,0,iw,!rl,orie,!r2,irev,qins,il))<<"@s"<<is++<<" legend \"lr\""<<endl;
  
  //VV
  */
  jvec V_corr_0=(load_qonly("V1V1_00",1,1,RE,+1,1)+load_qonly("V2V2_00",1,1,RE,+1,1)+load_qonly("V3V3_00",1,1,RE,+1,1))/3;
  const int LI=0;
  jvec V_corr_L=(load_qonly("V1V1_LL",LI,LI,RE,+1,1)+load_qonly("V2V2_LL",LI,LI,RE,+1,1)+load_qonly("V3V3_LL",LI,LI,RE,+1,1))/3;
  jvec V_corr_M=(load_qonly("V1V1_0M",LI,LI,RE,+1,1)+load_qonly("V2V2_0M",LI,LI,RE,+1,1)+load_qonly("V3V3_0M",LI,LI,RE,+1,1))/3;
  jvec V_corr_P=(load_qonly("V1V1_0P",LI,LI,IM,-1,1)+load_qonly("V2V2_0P",LI,LI,IM,-1,1)+load_qonly("V3V3_0P",LI,LI,IM,-1,1))/3;
  jvec V_corr_S=(load_qonly("V1V1_0S",LI,LI,RE,+1,1)+load_qonly("V2V2_0S",LI,LI,RE,+1,1)+load_qonly("V3V3_0S",LI,LI,RE,+1,1))/3;
  jvec V_corr_T=(load_qonly("V1V1_0T",LI,LI,RE,+1,1)+load_qonly("V2V2_0T",LI,LI,RE,+1,1)+load_qonly("V3V3_0T",LI,LI,RE,+1,1))/3;
  V_corr_0.print_to_file("plots/V.xmg");
  effective_mass(V_corr_0).print_to_file("plots/V_effmass.xmg");
  V_corr_L.print_to_file("plots/V_L.xmg");
  V_corr_M.print_to_file("plots/V_M.xmg");
  V_corr_P.print_to_file("plots/V_P.xmg");
  V_corr_S.print_to_file("plots/V_S.xmg");
  V_corr_T.print_to_file("plots/V_T.xmg");
  (V_corr_L/V_corr_0).print_to_file("plots/V_L_ratio.xmg");
  (V_corr_M/V_corr_0).print_to_file("plots/V_M_ratio.xmg");
  (V_corr_P/V_corr_0).print_to_file("plots/V_P_ratio.xmg");
  (V_corr_S/V_corr_0).print_to_file("plots/V_S_ratio.xmg");
  (V_corr_T/V_corr_0).print_to_file("plots/V_T_ratio.xmg");
  
  //compute HVP
  jvec HVP_0=calc_HVP(V_corr_0);
  
  write_HVP("plots/HVP_0.xmg",HVP_0);
  
  write_HVP("plots/HVP_L_ratio.xmg",calc_HVP(V_corr_L)/HVP_0,1);
  write_HVP("plots/HVP_M_ratio.xmg",calc_HVP(V_corr_M)/HVP_0,1);
  write_HVP("plots/HVP_P_ratio.xmg",calc_HVP(V_corr_P)/HVP_0,1);
  write_HVP("plots/HVP_S_ratio.xmg",calc_HVP(V_corr_S)/HVP_0,1);
  write_HVP("plots/HVP_T_ratio.xmg",calc_HVP(V_corr_T)/HVP_0,1);
  
  write_HVP("plots/HVP_L.xmg",calc_HVP(V_corr_L),1);
  write_HVP("plots/HVP_M.xmg",calc_HVP(V_corr_M),1);
  write_HVP("plots/HVP_P.xmg",calc_HVP(V_corr_P),1);
  write_HVP("plots/HVP_S.xmg",calc_HVP(V_corr_S),1);
  write_HVP("plots/HVP_T.xmg",calc_HVP(V_corr_T),1);
  
  //print weight function
  ofstream HVP_weight_out("plots/HVP_weight.xmg");
  ofstream HVP_0_weighted_out("plots/HVP_0_weighted.xmg");
  HVP_0_weighted_out<<"@type xydy"<<endl;
  for(int i=0;i<T*nsubpoints;i++)
    {
      double Q2=sqr(2*M_PI*i/(T*a*nsubpoints));
      double t=1/(1+log(5/Q2));
      //HVP_weight_out<<Q2<<" "<<" "<<HVP_kernel(Q2)<<endl;
      //HVP_0_weighted_out<<Q2<<" "<<" "<<HVP_kernel(Q2)*HVP_0[i]<<endl;
      if(Q2) HVP_0_weighted_out<<t<<" "<<Q2/sqr(t)*HVP_kernel(Q2)*HVP_0[i]<<endl;
    }

  /*
  //for debug
  // jvec_load("data/corr_V1P5_00",T,njacks,icombo_qonly(0,0,0,RE)).print_to_file("plots/V1P5_00.xmg");
  // jvec_load("data/corr_V1P5_0M",T,njacks,icombo_qonly(0,0,0,RE)).print_to_file("plots/V1P5_0M.xmg");
  // jvec_load("data/corr_V1P5_0S",T,njacks,icombo_qonly(0,0,0,RE)).print_to_file("plots/V1P5_0S.xmg");
  // jvec_load("data/corr_V1P5_0P",T,njacks,icombo_qonly(0,0,0,IM)).print_to_file("plots/V1P5_0P.xmg");
  // jvec_load("data/corr_V1P5_LL",T,njacks,icombo_qonly(0,0,0,RE)).print_to_file("plots/V1P5_LL.xmg");
  // jvec_load("data/corr_V1P5_0T",T,njacks,icombo_qonly(0,0,0,RE)).print_to_file("plots/V1P5_0T.xmg");
  */  
  return 0;
}
