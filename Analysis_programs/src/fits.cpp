#pragma once

#include "effmass.cpp"

int debug_fit=1;

//jack-vec version
jvec effective_mass(jvec a,int TH=-1,int par=1)
{
  if(TH==-1) TH=a.nel-1;

  int njack=a.njack;
  
  jvec b(a.nel-1,a.njack);
  
  for(int t=0;t<a.nel-1;t++)
    {
      jack temp=-log(a[t+1]/a[t]);
      double miniz=temp.med();
      double einiz=temp.err();
      if(einiz==0) einiz=1.e-10;
      
      for(int ijack=0;ijack<=njack;ijack++)
	b.data[t].data[ijack]=effective_mass(a[t][ijack],a[t+1][ijack],t,TH,miniz,einiz,par);
    }
  return b;
}

jvec antonin_effective_mass(jvec in,int par=1)
{
  jvec out(in.nel-2,in.njack);
  if(par==1) for(int t=1;t<in.nel-1;t++) out[t-1]=acosh((in[t-1]+in[t+1])/(2*in[t]));
  else       for(int t=1;t<in.nel-1;t++) out[t-1]=asinh((in[t-1]+in[t+1])/(2*in[t]));

  return out;
}

jvec numerical_derivative(jvec a)
{
  jvec b(a.nel-1,a.njack);
  
  for(int iel=0;iel<a.nel-1;iel++)
    b[iel]=a[iel+1]-a[iel];
  
  return b;
}

jvec simmetric_derivative(jvec a)
{
  jvec b(a.nel,a.njack);
  
  for(int iel=1;iel<a.nel-1;iel++) b[iel]=(a[iel+1]-a[iel-1])/2;
  b[0]=b[a.nel-1]=0;
  
  return b;
}

//jack-vec
jvec aperiodic_effective_mass(const jvec a)
{
  int TH=a.nel-1;
  int njack=a.njack;
  
  jvec b(TH,njack);
  
  for(int t=0;t<TH;t++) b.data[t]=log(a.data[t]/a.data[t+1]);
  
  return b;
}

//fit the mass
jack mass_fit(jvec corr,int tmin,int tmax,const char *path=NULL,int TH=-1,int parity=1)
{
  jvec effe=effective_mass(corr,TH,parity);
  jack mass=constant_fit(effe,tmin,tmax-1,path);
  
  return mass;
}

//fit the mass and the matrix element
void two_pts_fit(jack &E,jack &Z2,jvec corr,int tmin,int tmax,const char *path1=NULL,const char *path2=NULL,int TH=-1,int parity=1)
{
  E=mass_fit(corr,tmin,tmax,path1,TH,parity);
  jvec temp(corr.nel,corr.njack);
  if(TH==-1) TH=temp.nel-1;
  if(parity==1) for(int t=0;t<=TH;t++) temp[t]=corr[t]/exp(-E*TH)/cosh(E*(TH-t))*E;
  else          for(int t=0;t<=TH;t++) temp[t]=corr[t]/exp(-E*TH)/sinh(E*(TH-t))*E;
  Z2=constant_fit(temp,tmin,tmax,path2);
  
  if(debug_fit) cout<<"E: "<<E<<endl;
}

//fit the mass and the matrix element in SS and SL combo
double *c_two_pts_fit,*e_two_pts_fit;
int TH_two_pts_fit;
int tmin_two_pts_fit;
int tmax_two_pts_fit;

template <class T> T fun_two_pts_migrad_fit(T Z2,T M,double t)
{return Z2*exp(-M*TH_two_pts_fit)*cosh(M*(TH_two_pts_fit-t))/M;}

void ch2_two_pts_migrad_fit(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  double M=p[0];
  double Z2=p[1];
  
  for(int t=tmin_two_pts_fit;t<=min(tmax_two_pts_fit,TH_two_pts_fit);t++)
    {
      double num=c_two_pts_fit[t];
      double teo=fun_two_pts_migrad_fit(Z2,M,t);
      double diff=num-teo;
      double err=e_two_pts_fit[t];
      double cont=sqr(diff/err);
      ch+=cont;
      if(flag==3)
	cout<<" Z2: "<<Z2<<", M: "<<M<<", t="<<t<<", diff=("<<num<<"-"<<teo<<")="<<diff<<" err="<<err<<" cont="<<cont<<endl;
    }
}

void two_pts_migrad_fit(jack &M,jack &Z2,jvec corr,int tmin,int tmax,const char *path=NULL)
{
  TMinuit minu;
  minu.SetPrintLevel(-1);
  minu.SetFCN(ch2_two_pts_migrad_fit);

  jvec ecorr=effective_mass(corr);
  M=constant_fit(ecorr,tmin,tmax,NULL);
  double M_med=M.med();
  if(std::isnan(M_med)) M_med=1;
  
  jvec temp(corr.nel,corr.njack);
  int TH=temp.nel-1;
  for(int t=0;t<=TH;t++)
    temp[t]=corr[t]/exp(-M*TH)/cosh(M*(TH-t))*M;
  
  Z2=constant_fit(temp,tmin,tmax,NULL);
  double Z2_med=Z2.med();
  if(std::isnan(Z2_med)) Z2_med=M_med=1;
  minu.DefineParameter(0,"M",M_med,0.001,0,0);
  minu.DefineParameter(1,"Z2",Z2_med,0.001,0,0);
  
  int njack=M.njack;
  c_two_pts_fit=new double[TH+1];
  e_two_pts_fit=new double[TH+1];
  
  TH_two_pts_fit=TH;
  tmin_two_pts_fit=tmin;
  tmax_two_pts_fit=tmax;
  
  for(int iel=0;iel<=TH;iel++)
      e_two_pts_fit[iel]=corr[iel].err();
  
  for(int ijack_fit=0;ijack_fit<=njack;ijack_fit++)
    {
      //minu.FixParameter(0);
      for(int iel=0;iel<=TH;iel++)
	  c_two_pts_fit[iel]=corr[iel][ijack_fit];
      minu.Migrad();
      double dum;
      minu.GetParameter(0,M.data[ijack_fit],dum);
      minu.GetParameter(1,Z2.data[ijack_fit],dum);
    }
  
  double ch2,grad[2],par[2]={M[njack],Z2[njack]};
  minu.Eval(2,grad,ch2,par,3);
  if(debug_fit) cout<<"M: "<<M<<", ch2: "<<ch2<<endl;
  
  if(path!=NULL) write_constant_fit_plot(path,ecorr,M,tmin,tmax);
}


//fit the mass and the matrix element in SS and SL combo
double *c_two_pts_SL_fit[2],*e_two_pts_SL_fit[2];
double *c_two_pts_SL_fit_teo[2];
double *c_two_pts_SL_fit_ch2_contr[2];
int TH_two_pts_SL_fit;
int tmin_two_pts_SL_fit[2];
int tmax_two_pts_SL_fit[2];

double fun_two_pts_SL_fit(double Z1,double Z2,double M,double t)
{return Z1*Z2*exp(-M*TH_two_pts_SL_fit)*cosh(M*(TH_two_pts_SL_fit-t))/M;}

void ch2_two_pts_SL_fit(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  double M=p[0];
  double ZL=p[1];
  double ZS=p[2];
  
  for(int t=tmin_two_pts_SL_fit[0];t<=min(tmax_two_pts_SL_fit[0],TH_two_pts_SL_fit);t++)
    {
      double num=c_two_pts_SL_fit[0][t];
      double teo=c_two_pts_SL_fit_teo[0][t]=fun_two_pts_SL_fit(ZL,ZS,M,t);
      double diff=num-teo;
      double err=e_two_pts_SL_fit[0][t];
      double cont=c_two_pts_SL_fit_ch2_contr[0][t]=sqr(diff/err);
      ch+=cont;
      if(flag==3) cout<<"SL, t="<<t<<", diff=("<<num<<"-"<<teo<<")="<<diff<<" err="<<err<<" cont="<<cont<<endl;
    }
  
  for(int t=tmin_two_pts_SL_fit[1];t<=min(tmax_two_pts_SL_fit[1],TH_two_pts_SL_fit);t++)
    {
      double num=c_two_pts_SL_fit[1][t];
      double teo=c_two_pts_SL_fit_teo[1][t]=fun_two_pts_SL_fit(ZS,ZS,M,t);
      double diff=num-teo;
      double err=e_two_pts_SL_fit[1][t];
      double cont=c_two_pts_SL_fit_ch2_contr[1][t]=sqr(diff/err);
      ch+=cont;
      if(flag==3) cout<<"SS, t="<<t<<", diff="<<diff<<" err="<<err<<" cont="<<cont<<endl;
    }
}

void two_pts_SL_fit(jack &M,jack &ZL,jack &ZS,jvec corrSL,jvec corrSS,int tminSL,int tmaxSL,int tminSS,int tmaxSS,const char *path1=NULL,const char *path2=NULL,const char *path_cls=NULL)
{
  //perform the fit and set an initial estimate
  jack MSL,MSS,ZSL,ZSS;
  two_pts_fit(MSL,ZSL,corrSL,tminSL,tmaxSL);
  two_pts_fit(MSS,ZSS,corrSS,tminSS,tmaxSS);
  
  //get estimates for ZS and ZL
  M=MSL;
  ZS=sqrt(ZSS);
  ZL=ZSL/ZS;
  
  //define minimzer
  TMinuit minu;
  minu.SetPrintLevel(-1);
  minu.SetFCN(ch2_two_pts_SL_fit);
  
  //copy parameters
  int njack=MSL.njack;
  int TH=TH_two_pts_SL_fit=corrSL.nel-1;
  tmin_two_pts_SL_fit[0]=tminSL;
  tmin_two_pts_SL_fit[1]=tminSS;
  tmax_two_pts_SL_fit[0]=tmaxSL;
  tmax_two_pts_SL_fit[1]=tmaxSS;
  
  //define temporary structures
  c_two_pts_SL_fit[0]=new double[TH+1];
  c_two_pts_SL_fit[1]=new double[TH+1];
  c_two_pts_SL_fit_ch2_contr[0]=new double[TH+1];
  c_two_pts_SL_fit_ch2_contr[1]=new double[TH+1];
  c_two_pts_SL_fit_teo[0]=new double[TH+1];
  c_two_pts_SL_fit_teo[1]=new double[TH+1];
  e_two_pts_SL_fit[0]=new double[TH+1];
  e_two_pts_SL_fit[1]=new double[TH+1];
  
  //copy errors
  for(int iel=0;iel<=TH;iel++)
    {
      e_two_pts_SL_fit[0][iel]=corrSL[iel].err();
      e_two_pts_SL_fit[1][iel]=corrSS[iel].err();
    }
  
  //loop over jacknife
  for(int ijack_fit=0;ijack_fit<=njack;ijack_fit++)
    {
      //set pars
      minu.DefineParameter(0,"M",MSL[ijack_fit],MSL.err(),0,0);
      minu.DefineParameter(1,"ZL",ZL[ijack_fit],ZL.err(),0,0);
      minu.DefineParameter(2,"ZS",ZS[ijack_fit],ZS.err(),0,0);
      
      //copy elements
      for(int iel=0;iel<=TH;iel++)
	{
	  c_two_pts_SL_fit[0][iel]=corrSL[iel][ijack_fit];
	  c_two_pts_SL_fit[1][iel]=corrSS[iel][ijack_fit];
	}
      
      //minimize
      minu.Migrad();
      
      //get back parameters
      double dum;
      minu.GetParameter(0,M[ijack_fit],dum);
      minu.GetParameter(1,ZL[ijack_fit],dum);
      minu.GetParameter(2,ZS[ijack_fit],dum);
      
      //write down chi2
      if(ijack_fit==njack)
	if(path_cls!=NULL)
	  {
	    ofstream out(path_cls);
	    double tot=0;
	    int ndof=0;
	    out<<"SS/SL t ((teo-data)/err)^2=chi2_contr"<<endl;
	    for(int isl=0;isl<2;isl++)
	      {
		out<<"================================="<<endl;
		const char SLSS[2][4]={"SL","SS"};
		for(int it=tmin_two_pts_SL_fit[isl];it<=min(tmax_two_pts_SL_fit[isl],TH_two_pts_SL_fit);it++)
		  {
		    
		    out<<SLSS[isl]<<" "<<it<<" (("<<c_two_pts_SL_fit[isl][it]<<"-"<<c_two_pts_SL_fit_teo[isl][it]<<")/"<<
		      e_two_pts_SL_fit[isl][it]<<")^2="<<c_two_pts_SL_fit_ch2_contr[isl][it]<<endl;
		    ndof++;
		    tot+=c_two_pts_SL_fit_ch2_contr[isl][it];
		  }
	      }
	    out<<"================================="<<endl;
	    out<<"Total chi2: "<<tot<<"/"<<ndof-3<<"="<<tot/(ndof-3)<<endl;
	  }
    }
  
  //double ch2,grad[3],par[3]={M[njack],ZL[njack],ZS[njack]};
  //minu.Eval(3,grad,ch2,par,3);
  //cout<<"M: "<<smart_print(M)<<", ch2: "<<ch2<<endl;
  
  //write plots
  if(path1!=NULL)
    {
      write_constant_fit_plot(path1,effective_mass(corrSL),M,tminSL,tmaxSL);
      if(path2==NULL) append_constant_fit_plot(path1,effective_mass(corrSS),M,tminSS,tmaxSS,3);
    }
  if(path2!=NULL) write_constant_fit_plot(path2,effective_mass(corrSS),M,tminSS,tmaxSS);
  
  //delete
  delete[] c_two_pts_SL_fit[0];
  delete[] c_two_pts_SL_fit[1];
  delete[] c_two_pts_SL_fit_ch2_contr[0];
  delete[] c_two_pts_SL_fit_ch2_contr[1];
  delete[] c_two_pts_SL_fit_teo[0];
  delete[] c_two_pts_SL_fit_teo[1];
  delete[] e_two_pts_SL_fit[0];
  delete[] e_two_pts_SL_fit[1];
}

///////////////////////////////////////////////////////////////////////////////////////////

double *c_two_states_fit,*e_two_states_fit;
int TH_two_states_fit;
int tmin_two_states_fit;
int tmax_two_states_fit;

template <class Ti> Ti fun_two_states_migrad_fit(Ti MA,Ti MB,Ti ZA,Ti ZB,int t)
{return fun_two_pts_migrad_fit(ZA,MA,t)+fun_two_pts_migrad_fit(ZB,MB,t);}

void ch2_two_states_migrad_fit(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  double MA=p[0];
  double MB=p[1];
  double ZA=p[2];
  double ZB=p[3];
  
  for(int t=tmin_two_states_fit;t<=min(tmax_two_states_fit,TH_two_states_fit);t++)
    {
      double num=c_two_states_fit[t];
      double teo=fun_two_states_migrad_fit(MA,MB,ZA,ZB,t);
      double diff=num-teo;
      double err=e_two_states_fit[t];
      double cont=sqr(diff/err);
      ch+=cont;
      if(flag==3)
	cout<<" ZA: "<<ZA<<", MA: "<<MA<<"; ZB: "<<ZB<<", MB: "<<MB<<
	  "; t="<<t<<", diff=("<<num<<"-"<<teo<<")="<<diff<<" err="<<err<<" cont="<<cont<<endl;
    }
}

void two_states_fit(jack &MA,jack &MB,jack &ZA,jack &ZB,int tmin,int TH,jvec corr,const char *path)
{
  TMinuit minu;
  minu.SetPrintLevel(-1);
  minu.SetFCN(ch2_two_states_migrad_fit);
  minu.DefineParameter(0,"MA",MA.med(),MA.err(),MA.med()-3*MA.err(),MA.med()+3*MA.err());
  minu.DefineParameter(1,"MB",MB.med(),MB.err(),MB.med()-3*MB.err(),MB.med()+3*MB.err());
  minu.DefineParameter(2,"ZA",ZA.med(),ZA.err(),ZA.med()-3*ZA.err(),ZA.med()+3*ZA.err());
  minu.DefineParameter(3,"ZB",ZB.med(),ZB.err(),ZA.med()-3*ZB.err(),ZB.med()+3*ZB.err());
  
  TH_two_pts_fit=TH_two_states_fit=TH;
  c_two_states_fit=new double[TH+1];
  e_two_states_fit=new double[TH+1];
  
  tmin_two_states_fit=tmin;
  tmax_two_states_fit=TH;
  
  for(int iel=0;iel<=TH;iel++) e_two_states_fit[iel]=corr[iel].err();
  int njack=corr[0].njack;

  for(int ijack_fit=0;ijack_fit<=njack;ijack_fit++)
    {
      for(int iel=0;iel<=TH;iel++) c_two_states_fit[iel]=corr[iel][ijack_fit];
      minu.Migrad();
      double dum;
      minu.GetParameter(0,MA.data[ijack_fit],dum);
      minu.GetParameter(1,MB.data[ijack_fit],dum);
      minu.GetParameter(2,ZA.data[ijack_fit],dum);
      minu.GetParameter(3,ZB.data[ijack_fit],dum);
    }
  
  //separate plots
  //summ of exponentials
  jvec temp_fun(TH+1,njack);
  for(int t=0;t<=TH;t++)
    temp_fun[t]=fun_two_states_migrad_fit(MA,MB,ZA,ZB,t);
  temp_fun=effective_mass(temp_fun);
  
  ofstream foutB(path);
  //single mass
  foutB<<"@s0 line type 1"<<endl;      
  foutB<<"@s0 line color 8"<<endl;
  foutB<<"@s0 fill color 8"<<endl;
  foutB<<"@s0 fill type 1"<<endl;
  foutB<<"@type xy"<<endl;
  for(int i=tmin;i<TH;i++) foutB<<i<<" "<<MA.med()+MA.err()<<endl;
  for(int i=TH-1;i>=tmin;i--) foutB<<i<<" "<<MA.med()-MA.err()<<endl;
  foutB<<"&"<<endl;

  //summ of expo
  foutB<<"@s1 line type 1"<<endl;      
  foutB<<"@s1 line color 7"<<endl;
  foutB<<"@s1 fill color 7"<<endl;
  foutB<<"@s1 fill type 1"<<endl;
  foutB<<"@type xy"<<endl;
  for(int t=tmin;t<TH;t++) foutB<<t<<" "<<temp_fun[t].med()+temp_fun[t].err()<<endl;
  for(int t=TH-1;t>=tmin;t--) foutB<<t<<" "<<temp_fun[t].med()-temp_fun[t].err()<<endl;
  foutB<<"&"<<endl;
  //central line
  foutB<<"@s2 line color 1"<<endl;
  foutB<<"@type xy"<<endl;      
  for(int i=tmin;i<TH;i++) foutB<<i<<" "<<temp_fun[i].med()<<endl;


  //central line
  foutB<<"@s3 line color 3"<<endl;
  foutB<<"@type xy"<<endl;      
  for(int i=tmin;i<TH;i++) foutB<<i<<" "<<MA.med()<<endl;
  //plot the original data with error  
  foutB<<"&"<<endl;
  foutB<<"@type xydy"<<endl;      
  foutB<<"@s4 line type 0"<<endl;      
  foutB<<"@s4 symbol color 1"<<endl;
  foutB<<"@s4 errorbar color 1"<<endl;
  foutB<<"@s4 symbol 1"<<endl;
  foutB<<effective_mass(corr);
  foutB<<"&"<<endl;
}

//////////////////////////////////////////////////////////////////////////////////////

bvec lin_solve(double *A,bvec b)
{
  int d=b.nel;
  int nboot=b.nboot;
  int njack=b.njack;

  bvec x(d,nboot,njack);
  
  for(int i=0;i<d;i++)
    {
      double C=A[i*d+i];
      for(int j=i;j<d;j++) A[i*d+j]/=C;
      b[i]/=C;
      
      for(int k=i+1;k<d;k++)
        {
          double C=A[k*d+i];
          for(int j=i;j<d;j++) A[k*d+j]-=A[i*d+j]*C;
          b[k]-=C*b[i];
        }
    }
  
  for(int k=d-1;k>=0;k--)
    {
      boot S(nboot,njack);
      S=0;
      for(int i=k+1;i<d;i++) S+=A[k*d+i]*x[i];
      x[k]=b[k]-S;
    }
  
  return x;
}

bvec poly_fit(double *x,bvec y,int d,double xmin,double xmax)
{
  int np=y.nel;
  int nboot=y.nboot;
  int njack=y.njack;
  
  double Al[2*d+1];memset(Al,0,sizeof(double)*(2*d+1));
  bvec c(d+1,nboot,njack);c=0;

  for(int p=0;p<np;p++)
    if(x[p]<=xmax &&x[p]>=xmin)
      {
	//calculate the weight
	double w=pow(y[p].err(),-2);
	//compute Al and c
	for(int f=0;f<=2*d;f++)
	  {
	    Al[f]+=w;
	    if(f<=d) c[f]+=y[p]*w;
	    w*=x[p];
	  }
      }
  
  double A[(d+1)*(d+1)];
  for(int i=0;i<=d;i++)
    for(int j=0;j<=d;j++)
      A[i*(d+1)+j]=Al[i+j];
  
  return lin_solve(A,c);
}      

boot chi2_poly_fit(double *x,bvec y,int d,double xmin,double xmax,bvec pars)
{
  int np=y.nel;
  int nboot=y.nboot;
  int njack=y.njack;
  
  boot ch2(nboot,njack);
  int ndof=-pars.nel;
  ch2*=0;
  for(int ip=0;ip<np;ip++)
    if(x[ip]<=xmax && x[ip]>=xmin)
      {
	boot t=pars[0];
	double R=x[ip];
	for(int ipow=1;ipow<pars.nel;ipow++)
	  {
	    t+=pars[ipow]*R;
	    R*=x[ip];
	  }

	ch2+=sqr((t-y[ip])/y[ip].err());
	ndof++;
      }
  
  return ch2/ndof;
}      
