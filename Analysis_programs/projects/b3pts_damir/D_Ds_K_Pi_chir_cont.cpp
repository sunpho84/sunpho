#include "include.h"
#include "../nf2/common_pars.cpp"

int nth,nens,nth_to_use;
const int nbeta=4,nth_max=8,nens_max=13,njacks=16,nboots=100,ncoeff_z=2;
const double MD_ph=1.868,MP_ph=0.135,EP_ph=(MD_ph*MD_ph+MP_ph*MP_ph)/(2*MD_ph);
int ibeta[nens_max],use[nens_max];
double lmass[nens_max];
const int start_th[6]={1,0,1,1,0,1};

bvec Q2,FP,F0,FT,F0S,Z,MD,MP,EP,FNU,FPE,FP_WP,F0_WP,F0S_WP;

int ref_ml_beta[4]={-1,-1,-1,-1};

const char set_color[nbeta][1024]={"black","blue","red","green4"};
const char set_fill_color[nbeta][1024]={"grey","turquoise","yellow","green"};
const char set_symbol[nbeta][1024]={"square","circle","triangle up","triangle left"};
const char set_legend[nbeta][1024]={"\\xb\\0=3.80","\\xb\\0=3.90","\\xb\\0=4.05","\\xb\\0=4.20"};
const char set_legend_fm[nbeta][1024]={"a = 0.098 fm","a = 0.085 fm","a = 0.067 fm","a = 0.054 fm"};

double MDSTAR=2.01,MD0STAR=2.32; //add D to make diff in the case D->Pi
double MDSSTAR=2.112,MD0SSTAR=2.32;
int plot_iboot;
double DeltaVP=0.137; //V-P 
int include_log=0,include_f0s=1,nff=2+include_f0s;
int Z_expansion=0;
int parametrization=7,npars_poss_non_Z[]={8,9,10,11,13,12,13,8},npars_poss_Z[]={8};

int icombo(int iens,int ith)
{return iens*nth+ith;}

double theta(double x)
{return x>=0;}

boot theta(boot in)
{
  boot out=in;
  for(int i=0;i<=nboot;i++) out.data[i]=theta(in[i]);
  return out;
}

template <class T> T F18(T x)
{
  T x2=x*x,c=1-x2,t=theta(c);
  c*=2*t-1;
  
  return t*
    (sqrt(c)*(-0.5*log(x2)+log(1+sqrt(c))))
    +(1-t)*
    (-sqrt(c))*(M_PI/2-atan(1/sqrt(c)))
    ;
}

const double f0_int=0.122,la=4*M_PI*f0_int,gc=0.54,g2=gc*gc;

template <class T> T I1(T m)
{return 2*m*m*log(m/la);}
template <class T> T I2(T m,T E)
{return -4*E*E*log(m/la)-4*E*E*F18(m/E)+2*E*E;}
template <class T> T J1(T m,T E)
{return 2*(-m*m+2.0/3*E*E)*log(m/la)+4.0/3*(E*E-m*m)*F18(m/E)-10.0/9*E*E+4.0/3*m*m;}

template <class T> T DFP(T m,T E)
{return (-3.0/4*(1+3*g2)*I1(m)+4*g2*J1(m,E)-8*M_PI*g2/3*m*m*m/E)/la/la;}
template <class T> T DFV(T m,T E)
{return ((5-9*g2)/4*I1(m)+2*I2(m,E))/la/la;}

template <class T1,class T2,class T3> T1 b_z(T1 *p,T2 a1,T3 MP)
{return p[0]*(1+p[1]*a1*a1+p[2]*MP*MP+p[3]*MP*MP*log(MP*MP));}

//fit for FPE (iff==0) or FNU (iff==1)
template <class T1,class T2,class T3> T1 fun_fit_F(int iff,T1 *p,T2 *x,T3 a)
{
  T2 MP=x[0];
  T3 a1=a/lat[1][nboots];
  
  if(Z_expansion)
    {
      T2 Z=x[1];
      T1 out=p[0]*0;
      int ibase;

      if(iff==1) ibase=ncoeff_z+1;
      
      for(int i=0;i<=ncoeff_z;ibase++) out+=pow(Z,2*i)*b_z(p+4*(ibase+i),a1,MP);
      
      return out;
    }
  else
    {
      T2 EP=x[1];
      T2 logP=DFP(MP,EP)*include_log;
      T2 logV=DFV(MP,EP)*include_log;
      
      switch(parametrization)
	{
	case 0: //fv misses a visible E^2 dependence
	  if(iff==0) return p[0]/(EP+DeltaVP)*(1+p[1]*MP*MP+p[2]*EP+p[3]*a1*a1);
	  else       return p[4]*(1+p[5]*MP*MP+p[6]*EP+p[7]*a1*a1);
	  break;
	  //////
	case 1: //minimal dependence for fv satisfied
	  if(iff==0) return p[0]/(EP+DeltaVP)*(1+p[1]*MP*MP+p[2]*EP+p[3]*a1*a1);
	  else       return p[4]*(1+p[5]*MP*MP+p[6]*EP+p[7]*a1*a1+p[8]*EP*EP);
	  break;
	  ///// 
	case 2: //MP*MP dependence added
	  if(iff==0) return p[0]/(EP+DeltaVP)*(1+p[1]*MP*MP+p[2]*EP+p[3]*a1*a1);
	  else       return p[4]*(1+p[5]*MP*MP+p[6]*EP+p[7]*a1*a1+p[8]*EP*EP+p[9]*EP*MP*MP);
	  break;
	  /////
	case 3: //wout ch logs work
	  if(iff==0) return p[0]/(EP+DeltaVP)*(1+logP+p[1]*MP*MP+p[2]*EP+p[3]*a1*a1+p[4]*EP*EP);
	  else       return p[5]*(1+logV+p[6]*MP*MP+p[7]*EP+p[8]*a1*a1+p[9]*EP*EP+p[10]*EP*MP*MP);
	  break;
	  /////
	case 4: //this makes error huge and not justified
	  if(iff==0) return p[0]/(EP+DeltaVP)*(1+p[1]*MP*MP+p[2]*EP+p[3]*a1*a1+p[4]*EP*EP+p[11]*EP*a1*a1);
	  else       return p[5]*(1+p[6]*MP*MP+p[7]*EP+p[8]*a1*a1+p[9]*EP*EP+p[10]*EP*MP*MP+p[12]*EP*a1*a1);
	  break;
	case 5: // adding E^4 in fv to cancel logs
	  if(iff==0) return p[0]/(EP+DeltaVP)*(1+logP+p[1]*MP*MP+p[2]*EP+p[3]*a1*a1+p[4]*EP*EP);
	  else       return p[5]*(1+logV+p[6]*MP*MP+p[7]*EP+p[8]*a1*a1+p[9]*EP*EP+p[10]*EP*MP*MP+p[11]*EP*EP*EP*EP);
	  break;
	case 6: // adding E^3 and E^4 in fv to cancel logs
	  if(iff==0) return p[0]/(EP+DeltaVP)*(1+logP+p[1]*MP*MP+p[2]*EP+p[3]*a1*a1+p[4]*EP*EP);
	  else       return p[5]*(1+logV+p[6]*MP*MP+p[7]*EP+p[8]*a1*a1+p[9]*EP*EP+p[10]*EP*MP*MP+(p[11]*MP*MP+p[12]*EP)*EP*EP);
	  break;
	case 7: //testing
	  if(iff==0) return p[0]*(1+p[1]*MP*MP+p[2]*EP+p[3]*a1*a1);
	  else       return p[4]*(1+p[5]*MP*MP+p[6]*EP+p[7]*a1*a1);
	  break;
	  /////
	  /////
	  ///// default case
	default:
	  crash("should not reach here");
	  return p[0]*0;
	  break;
	}
    }
  
  return p[0]*0;
}

template <class T1,class T2,class T3> void FP0_to_PENU(T1 &FPE,T1 &FNU,T1 &FP,T1 &F0,T2 &MD,T2 &MP,T3 &EP,T3 &Q2)
{
  T1 t=(MD*MD-MP*MP)/Q2*(FP-F0);
  FNU=(MD+EP)*FP-(MD-EP)*t;
  FNU/=sqrt(2*MD);
  
  FPE=FP+t;
  FPE/=sqrt(2*MD);
}

template <class T1,class T2,class T3> void FPENU_to_FP0(T1 &FP,T1 &F0,T1 &FPE,T1 &FNU,T2 &MD,T2 &MP,T3 &EP)
{
  T2 PPE2=-(EP*EP-MP*MP);
  FP=FNU+(MD-EP)*FPE;
  FP/=sqrt(2*MD);
  
  F0=(MD-EP)*FNU-PPE2*FPE;
  F0*=sqrt(2*MD)/(MD*MD-MP*MP);
}

template <class T1,class T2> T1 fun_Z(T1 MP,T1 MD,T2 Q2)
{
  double T0=(MD_ph+MP_ph)*sqr(sqrt(MD_ph)-sqrt(MP_ph)); //0
  T1 TP=sqr(MD+MP);
  T1 Z=sqrt(TP-Q2)-sqrt(TP-T0);
  Z/=sqrt(TP-Q2)+sqrt(TP-T0);
  
  return Z;
}

double X_fit[nens_max*nth_max][2];
double Y_fit[nens_max*nth_max][3],err_Y_fit[nens_max*nth_max][3];
double ext_ch2;

//calculate the chi square
int contr_flag=0;
double chi2(double *p,int npars)
{
  double ch2=0;
  
  double Y_teo[nens*nth][3];
  
  for(int iff=0;iff<nff;iff++)
    for(int iens=0;iens<nens;iens++)
      if(use[iens])
	for(int ith=start_th[iff];ith<nth_to_use;ith++)
	  {
	    double a1=lat[ibeta[iens]][plot_iboot]/lat[1][nboots];
	    int ic=icombo(iens,ith);
	    if(iff!=2) Y_teo[ic][iff]=fun_fit_F(iff,p,X_fit[ic],lat[ibeta[iens]][plot_iboot]);
	    else
	      {
		if(Z_expansion)
		  {
		    //tbi
		  }
		else
		  {
		    double dum;
		    FPENU_to_FP0(dum,Y_teo[ic][iff],Y_teo[ic][0],Y_teo[ic][1],MD[iens].data[plot_iboot],MP[iens].data[plot_iboot],EP[icombo(iens,ith)].data[plot_iboot]);
		    Y_teo[ic][iff]+=p[npars-1]*a1*a1;
		  }
	      }
	    double contr=pow((Y_fit[ic][iff]-Y_teo[ic][iff])/err_Y_fit[ic][iff],2);
	    ch2+=contr;
	    if(contr_flag==1)
	      cout<<"contr"<<iff<<" ("<<iens<<","<<ith<<"): "<<contr<<" = ("<<Y_fit[ic][iff]<<" - "<<
		Y_teo[ic][iff]<<") / "<<err_Y_fit[ic][iff]<<endl;
	  }
  
  //export ch2
  ext_ch2=ch2;
  
  return ch2;
}

//wrapper for the calculation of the chi2
void chi2_wr(int &npar,double *fuf,double &ch,double *p,int flag)
{ch=chi2(p,npar);}

void fit(bvec &pars,bvec &F0,bvec &F1,bvec &F2,bvec &X0,bvec &X1)
{
  //copy Y errs
  for(int i=0;i<nth*nens;i++)
    {
      err_Y_fit[i][0]=F0[i].err();
      err_Y_fit[i][1]=F1[i].err();
      if(include_f0s) err_Y_fit[i][2]=F2[i].err();
    }
  
  //set minuit function
  TMinuit minu;
  minu.SetFCN(chi2_wr);
  //minu.SetPrintLevel(-1);
  
  //set pars
  int npars=pars.nel;
  for(int ipar=0;ipar<npars;ipar++) minu.DefineParameter(ipar,combine("P%d",ipar).c_str(),0,0.0001,0,0);
  
  for(int iboot=0;iboot<=nboots;iboot++)
    {
      plot_iboot=iboot;
      //quiet if not iboot 0
      if(iboot>0) minu.SetPrintLevel(-1);
      
      for(int i=0;i<nens*nth;i++)
	{
	  X_fit[i][0]=X0[i/nth][iboot]; //MP only defined for iens
	  X_fit[i][1]=X1[i][iboot];
	  Y_fit[i][0]=F0[i][iboot];
	  Y_fit[i][1]=F1[i][iboot];
	  Y_fit[i][2]=F2[i][iboot];
	}
      
      //minimize
      if(iboot==0)
	{
	  minu.Migrad();
	  minu.mnimpr();
	}
      minu.Migrad();
      
      //get back parameters
      double dum;
      for(int ipar=0;ipar<npars;ipar++)
	minu.GetParameter(ipar,pars.data[ipar].data[iboot],dum);
      
      //make it print separate contributions
      if(iboot==0)
	{
	  double t[npars];
	  for(int ipar=0;ipar<npars;ipar++) t[ipar]=pars.data[ipar].data[iboot];
	  contr_flag=1;
	  chi2(t,npars);
	  contr_flag=0;
	}
    }
  
  //print pars
  for(int ipar=0;ipar<npars;ipar++) cout<<"P"<<ipar<<"=("<<pars[ipar]<<")"<<endl;
  cout<<endl;
  
  //count data and print ch2
  int npoints=((2+include_f0s)*nth-1)*nens;
  for(int iens=0;iens<nens;iens++) if(!use[iens]) npoints-=(2+include_f0s)*nth-1;
  cout<<"Chi2 = "<<ext_ch2<<" / "<<npoints-npars<<" = "<<ext_ch2/(npoints-npars)<<endl;
}

void boot_from_jack(boot &out,jack in,int *iboot_jack)
{
  int nboot=out.nboot;
  int njack=in.njack;
  for(int iboot=0;iboot<nboot;iboot++) out.data[iboot]=in.data[iboot_jack[iboot]];
  out.data[nboot]=in.data[njack];
}

void load_iboot(int *iboot_jack,char *ens_name)
{
  char path[1024];
  sprintf(path,"/Users/francesco/QCD/LAVORI/RADIATIVE_CHARMONIUM/DATA1/%s/iboot",ens_name);
  
  FILE *fiboot=fopen(path,"r");
  if(fiboot==NULL)
    {
      perror(combine("Error opening file iboot for ensamble %s",path).c_str());
      exit(1);
    }
  int nr=fread(iboot_jack,sizeof(int),100,fiboot);
  if(nr!=100)
    {
      perror(combine("Error loading iboot data for ensamble %s",path).c_str());
      exit(1);
    }
  fclose(fiboot);
}

int main(int narg,char **arg)
{
  init_latpars();
  
  //read ensemble list, meson masses and meson name
  FILE *an_input_file=open_file("analysis_pars","r");
  char base_path_ff[200],base_path_MD[200],base_path_MP[200];
  
  read_formatted_from_file_expecting(base_path_ff,an_input_file,"%s","base_path_ff");
  read_formatted_from_file_expecting(base_path_MD,an_input_file,"%s","base_path_MD");
  read_formatted_from_file_expecting(base_path_MP,an_input_file,"%s","base_path_MP");
  
  read_formatted_from_file_expecting((char*)(&nth),an_input_file,"%d","nth");
  read_formatted_from_file_expecting((char*)(&nth_to_use),an_input_file,"%d","nth_to_use");
  read_formatted_from_file_expecting((char*)(&nens),an_input_file,"%d","nens");
  if(nth>nth_max) crash("nth>nth_max");
  if(nth_to_use>nth) crash("nth_to_use>nth");
  if(nens>nens_max) crash("nens>nens_max");
  
  FP=F0=FT=F0S=FP_WP=F0_WP=F0S_WP=Z=EP=FNU=FPE=Q2=bvec(nens*nth,nboots,njacks);
  MD=MP=bvec(nens,nboots,njacks);
  
  //read data
  ofstream bare_data_table("bare_data_table");
  for(int iens=0;iens<nens;iens++)
    {
      char ens_path[1024];
      read_formatted_from_file((char*)&(use[iens]),an_input_file,"%d","use");
      read_formatted_from_file((char*)&(ibeta[iens]),an_input_file,"%d","ibeta");
      read_formatted_from_file((char*)&(lmass[iens]),an_input_file,"%lg","lmass");
      read_formatted_from_file(ens_path,an_input_file,"%s","ens_path");
      
      //set the heavier mass
      int ib=ibeta[iens];
      if(use[iens]) if(ref_ml_beta[ib]==-1||lmass[iens]>lmass[ref_ml_beta[ib]]) ref_ml_beta[ib]=iens;

      //load iboot
      int iboot_jack[100];
      load_iboot(iboot_jack,ens_path);
      
      //read data
      jvec EP_j(nth,njacks),Q2_j(nth,njacks),FP_j(nth,njacks),F0_j(nth,njacks),FT_j(nth,njacks),F0S_j(nth,njacks);
      EP_j.load(combine(base_path_ff,ens_path).c_str(),0);
      Q2_j.load(combine(base_path_ff,ens_path).c_str(),1);
      FP_j.load(combine(base_path_ff,ens_path).c_str(),2);
      F0_j.load(combine(base_path_ff,ens_path).c_str(),3);
      FT_j.load(combine(base_path_ff,ens_path).c_str(),4);
      F0S_j.load(combine(base_path_ff,ens_path).c_str(),5);
      
      //convert jack to boot and pass to dimensionful quantity
      jack MD_j(njacks),MP_j(njacks);
      MD_j.load(combine(base_path_MD,ens_path).c_str(),0);
      MP_j.load(combine(base_path_MP,ens_path).c_str(),0);
      boot_from_jack(MD[iens],MD_j,iboot_jack);
      boot_from_jack(MP[iens],MP_j,iboot_jack);
      MD[iens]/=lat[ib];
      MP[iens]/=lat[ib];
      
      //write the bare data table while converting jack to boot
      for(int ith=0;ith<nth;ith++)
	{
	  int ic=icombo(iens,ith),ib=ibeta[iens];
	  boot_from_jack(EP.data[ic],EP_j[ith],iboot_jack);
	  boot_from_jack(Q2.data[ic],Q2_j[ith],iboot_jack);
	  boot_from_jack(FP.data[ic],FP_j[ith],iboot_jack);
	  boot_from_jack(F0.data[ic],F0_j[ith],iboot_jack);
	  boot_from_jack(FT.data[ic],FT_j[ith],iboot_jack);
	  boot_from_jack(F0S.data[ic],F0S_j[ith],iboot_jack);
	  
	  //pass to dimensionful quantity
	  Q2[ic]/=sqr(lat[ib]);
	  
	  EP[ic]/=lat[ib];
	  
	  Z[ic]=fun_Z(MP[iens],MD[iens],Q2[ic]);
	  
	  //compute FNU and FPE
	  FP0_to_PENU(FPE[ic],FNU[ic],FP[ic],F0[ic],MD[iens],MP[iens],EP[ic],Q2[ic]);
	  FP_WP[ic]=FP[ic]*(1-Q2[ic]/sqr(MDSTAR));
	  F0_WP[ic]=F0[ic]*(1-Q2[ic]/sqr(MD0STAR));
	  F0S_WP[ic]=F0S[ic]*(1-Q2[ic]/sqr(MD0STAR));
	  
	  bare_data_table<<ens_path<<" "<<EP[ic]<<" "<<Q2[ith]<<" "<<FP[ic]<<" "<<F0[ic]<<" "<<FT[ic]<<" "<<endl;
	}
    }
  fclose(an_input_file);
  
  //print ref ml mass
  cout<<"---"<<endl;
  for(int ib=0;ib<nbeta;ib++) cout<<"Ref "<<ib<<" = "<<ref_ml_beta[ib]<<", "<<lmass[ref_ml_beta[ib]]<<endl;
  cout<<"---"<<endl;
  
  //////////////////////////////////////// perform the fit ////////////////////////////////
  
  //perform the fit
  int npars=(Z_expansion?npars_poss_Z:npars_poss_non_Z)[parametrization]+include_f0s;
  bvec pars(npars,nboots,njacks);
  if(Z_expansion) fit(pars,FP_WP,F0_WP,F0S_WP,MP,EP);
  else            fit(pars,FPE,FNU,F0S,MP,EP);
  
  //extrapolate to the continuum, physical Pi
  double p2_ref[nth],EP_ref[nth],Q2_ref[nth],Z_ref[nth],x_ref[3][nth];
  bvec FPE_ph_ref(nth,nboot,njack),FNU_ph_ref(nth,nboot,njack);
  bvec FP_ph_ref(nth,nboot,njack),F0_ph_ref(nth,nboot,njack);
  bvec FP_WP_ph_ref(nth,nboot,njack),F0_WP_ph_ref(nth,nboot,njack);
  bvec F_ph_ref[6];
  for(int ith=0;ith<nth;ith++)
    {
      //compute ref points
      int iens_base=4; //compute back p2 for a given ens
      p2_ref[ith]=(sqr(EP[iens_base*nth+ith])-sqr(MP[iens_base])).med();
      x_ref[0][ith]=EP_ref[ith]=sqrt(sqr(MP_ph)+p2_ref[ith]);
      x_ref[1][ith]=Q2_ref[ith]=sqr(MD_ph)+sqr(MP_ph)-2*MD_ph*EP_ref[ith];
      x_ref[2][ith]=Z_ref[ith]=fun_Z(MP_ph,MD_ph,Q2_ref[ith]);
      
      //compute y
      if(!Z_expansion)
	{
	  double x_ref[2]={MP_ph,EP_ref[ith]};
	  FPE_ph_ref[ith]=fun_fit_F(0,pars.data,x_ref,0);
	  FNU_ph_ref[ith]=fun_fit_F(1,pars.data,x_ref,0);
	  FPENU_to_FP0(FP_ph_ref[ith],F0_ph_ref[ith],FPE_ph_ref[ith],FNU_ph_ref[ith],MD_ph,MP_ph,EP_ref[ith]);
	  FP_WP_ph_ref[ith]=FP_ph_ref[ith]*(1-Q2_ref[ith]/sqr(MDSTAR));
	  F0_WP_ph_ref[ith]=F0_ph_ref[ith]*(1-Q2_ref[ith]/sqr(MDSTAR));
	}
      else
	{
	  /*
	    double x_ph[]={MP_ph,Z_ref[ith]};
	    FIX
	    FP_WP_ph_ref[ith]=fun_fit_F(0,pars.data,x_ref,0);
	    F0_WP_ph_ref[ith]=fun_fit_F(1,pars.data,x_ref,0);
	    FP_ph_ref[ith]=FP_WP_ph_ref[ith]/(1-Q2_ref[ith]/sqr(MDSTAR));
	    F0_ph_ref[ith]=F0_WP_ph_ref[ith]/(1-Q2_ref[ith]/sqr(MDSTAR));
	    FP0_to_FPENU(FPE_ph_ref[ith],FNU_ph_ref[ith],FP_ph_ref[ith],F0_ph_ref[ith],MD_ph,MP_ph,EP_ref[ith]);
	  */
	}
      
      cout<<"ref ith "<<ith<<": "<<EP_ref[ith]<<": "<<FP_ph_ref[ith]<<endl;
      F_ph_ref[0]=FPE_ph_ref;
      F_ph_ref[1]=FNU_ph_ref;
      F_ph_ref[2]=FP_ph_ref;
      F_ph_ref[3]=F0_ph_ref;
      F_ph_ref[4]=FP_WP_ph_ref;
      F_ph_ref[5]=F0_WP_ph_ref;
    }
  
  //extrapolate to the continuum, physical PI and Q2=0
  boot FP_ph_cont,F0_ph_cont;
  boot FP_WP_ph_cont,F0_WP_ph_cont;
  boot FPE_ph_cont,FNU_ph_cont;
  if(!Z_expansion)
    {
      double x_ph[]={MP_ph,EP_ph};
      FPE_ph_cont=fun_fit_F(0,pars.data,x_ph,0);
      FNU_ph_cont=fun_fit_F(1,pars.data,x_ph,0);
      FPENU_to_FP0(FP_ph_cont,F0_ph_cont,FPE_ph_cont,FNU_ph_cont,MD_ph,MP_ph,EP_ph);
      FP_WP_ph_cont=FP_ph_cont;
      F0_WP_ph_cont=F0_ph_cont;
    }
  else
    {
      /* fix
      double x_ph[]={MP_ph,Z_ph};
      FP_WP_ph_cont=fun_fit_F(0,pars.data,x_ph,0);
      F0_WP_ph_cont=fun_fit_F(1,pars.data,x_ph,0);
      FP_ph_cont=FP_WP_ph_cont;
      F0_ph_cont=F0_WP_ph_cont;
      FP0_to_FPENU(FPE_ph_cont,FNU_ph_cont,FP_ph_cont,F0_ph_cont,MD_ph,MP_ph,EP_ph);
      */
    }
  cout<<"phys point in the continuum: "<<EP_ph<<": "<<FP_ph_cont<<endl;
  boot F_ph_cont[6]={FPE_ph_cont,FNU_ph_cont,FP_ph_cont,F0_ph_cont,FP_WP_ph_cont,F0_WP_ph_cont};
  
  //////////////////////////////// plot fpe and fnu as a function of energy ///////////////
    
  char tag_ff[6][4]={"pe","nu","p","0","pwp","0wp"};
  char tag_EQ2Z[3][4]={"E","Q2","Z"};
  for(int flag_EQ2Z=0;flag_EQ2Z<3;flag_EQ2Z++)
    for(int iff=0;iff<6;iff++)
      for(int iens=0;iens<=nens;iens++)
	{
	  //open appropriate file
	  ofstream out(combine("plots/f%s_ens_%d_fun_%s.xmg",tag_ff[iff],iens,tag_EQ2Z[flag_EQ2Z]).c_str());
	  
	  //non-continuum point
	  if(iens!=nens)
	    {
	      //write data
	      out<<"@type xydy"<<endl;
	      for(int ith=start_th[iff];ith<nth_to_use;ith++)
		{
		  int ic=icombo(iens,ith);
		  switch(flag_EQ2Z)
		    {
		    case 0: out<<EP[ic].med()<<" ";break;
		    case 1: out<<Q2[ic].med()<<" ";break;
		    case 2: out<<Z[ic].med()<<" ";break;
		    }
		  
		  switch(iff)
		    {
		    case 0: out<<FPE[ic];break;
		    case 1: out<<FNU[ic];break;
		    case 2: out<<FP[ic];break;
		    case 3: out<<F0[ic];break;
		    case 4: out<<FP_WP[ic];break;
		    case 5: out<<F0_WP[ic];break;
		    }
		  out<<endl;
		}
	    }

	  //compute fit line by points
	  int cont_latt_flag=!(iens==nens);
	  int chir_latt_flag=!(iens==nens);
	  double MP_plot_poss[2]={MP_ph,MP[min(iens,nens-1)][nboots]},
	    MD_plot_poss[2]={MD_ph,MD[min(iens,nens-1)][nboots]};
	  double MP_plot=MP_plot_poss[chir_latt_flag],MD_plot=MD_plot_poss[chir_latt_flag];
	  int np=150;
	  double x[2]={MP_plot,0},xs=MP_plot,xe=1.0,dx=(xe-xs)/(np-1);
	  double outx[np];
	  bvec Y(np,nboots,njacks);
	  
	  //compute the curve
	  for(int ip=0;ip<np;ip++)
	    {
	      //set energy
	      x[1]=xs+dx*ip;

	      //compute all possible x
	      double q2x=(MD_plot*MD_plot+MP_plot*MP_plot-2*MD_plot*x[1]);
	      switch(flag_EQ2Z)
		{
		case 0: outx[ip]=x[1];break;
		case 1: outx[ip]=q2x;break;
		case 2: outx[ip]=fun_Z(MP_plot,MD_plot,q2x);break;
		}
	      
	      //set a
	      boot a(nboot,njacks);
	      if(cont_latt_flag==0) a=0.0;else a=lat[ibeta[iens]];
	      
	      //compute y
	      if(iff<2||iff>=4) Y[ip]=fun_fit_F(iff,pars.data,x,a);
	      else
		{
		  boot t[2],o[2];
		  for(int i=0;i<2;i++) t[i]=fun_fit_F(i,pars.data,x,a);
		  FPENU_to_FP0(o[0],o[1],t[0],t[1],MD_plot,MP_plot,x[1]);
		  Y[ip]=o[iff-2];
		}
	    }
	  
	  //write polygon
	  out<<"&\n@type xy\n";
	  for(int ip=0;ip<np;ip++) out<<outx[ip]<<" "<<Y[ip].med()-Y[ip].err()<<endl;
	  for(int ip=np-1;ip>=0;ip--) out<<outx[ip]<<" "<<Y[ip].med()+Y[ip].err()<<endl;
	  out<<"&\n@type xy\n";
	  for(int ip=0;ip<np;ip++) out<<outx[ip]<<" "<<Y[ip].med()<<endl;
	  
	  //add the point at ref p2 for continuum curve
	  if(iens==nens)
	    {
	      out<<"&\n@type xydy\n";
	      for(int ith=0;ith<nth;ith++)
		out<<x_ref[flag_EQ2Z][ith]<<" "<<F_ph_ref[iff][ith]<<endl;
	    }
	  
	  //add the point at q2=0
	  {
	    double q2_0=0,E_0=(MD_plot*MD_plot+MP_plot*MP_plot)/(2*MD_plot);
	    double outx_0;
	    switch(flag_EQ2Z)
	      {
	      case 0: outx_0=E_0;break;
	      case 1: outx_0=q2_0;break;
	      case 2: outx_0=fun_Z(MP_plot,MD_plot,q2_0);break;
	      }
	    
	    out<<"&\n@type xydy\n";
	    out<<outx_0<<" "<<F_ph_cont[iff]<<endl;
	  }
	}
  
  return 0;
  
}
