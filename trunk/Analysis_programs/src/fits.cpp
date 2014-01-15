#pragma once

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
	{
	  double m=miniz;
	  double e=einiz;
	  
	  double targ=a[t+1][ijack]/a[t][ijack];
	  
	  double yl;
	  double yr;
	  
	  //increment the range up to reaching opposite sign
	  int q;
	  do
	    {
	      if(par==1)
		{
		  yl=cosh((m-e)*(TH-(t+1)))/cosh((m-e)*(TH-t))-targ;
		  yr=cosh((m+e)*(TH-(t+1)))/cosh((m+e)*(TH-t))-targ;
		}
	      else
		{
		  yl=sinh((m-e)*(TH-(t+1)))/sinh((m-e)*(TH-t))-targ;
		  yr=sinh((m+e)*(TH-(t+1)))/sinh((m+e)*(TH-t))-targ;
		}
	      q=((yl<0 && yr<0) || (yl>=0 && yr>=0));
	      //cout<<t<<" "<<ijack<<" "<<yl<<" "<<yr<<" "<<e<<endl;
	      if(q)
		{
		  e*=2;
		  if(m<=e) m+=(e-m);
		}
	    }
	  while(q);
	  
	  //bisect
	  double xl=m-e,xr=m+e,ym;
	  do
	    {
	      m=(xl+xr)/2;
	      if(par==1) ym=cosh(m*(TH-(t+1)))/cosh(m*(TH-t))-targ;
	      else       ym=sinh(m*(TH-(t+1)))/sinh(m*(TH-t))-targ;
	      if((yl<0 && ym<0) || (yl>0 && ym>0))
		xl=m;
	      else
		xr=m;
	      //cout<<t<<" "<<ijack<<" "<<m<<endl;
	    }
	  while(fabs(ym)>1.e-14);
	  b.data[t].data[ijack]=m;
	}
    }
  return b;
}

jvec numerical_derivative(jvec a)
{
  jvec b(a.nel-1,a.njack);
  
  for(int iel=0;iel<a.nel-1;iel++)
    b[iel]=a[iel+1]-a[iel];
  
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
  
  cout<<"E: "<<E<<endl;
}

//fit the mass and the matrix element in SS and SL combo
double *c_two_pts_fit,*e_two_pts_fit;
int TH_two_pts_fit;
int tmin_two_pts_fit;
int tmax_two_pts_fit;

double fun_two_pts_migrad_fit(double Z2,double M,double t)
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
      if(flag==3) cout<<" t="<<t<<", diff=("<<num<<"-"<<teo<<")="<<diff<<" err="<<err<<" cont="<<cont<<endl;
    }
}

void two_pts_migrad_fit(jack &M,jack &Z2,jvec corr,int tmin,int tmax,const char *path=NULL)
{
  jvec ecorr=effective_mass(corr);
  
  M=constant_fit(ecorr,tmin,tmax,NULL);
  jvec temp(corr.nel,corr.njack);
  int TH=temp.nel-1;
  for(int t=0;t<=TH;t++)
    temp[t]=corr[t]/exp(-M*TH)/cosh(M*(TH-t))*M;
  
  Z2=constant_fit(temp,tmin,tmax,NULL);
  
  TMinuit minu;
  minu.SetPrintLevel(-1);
  minu.SetFCN(ch2_two_pts_migrad_fit);
  minu.DefineParameter(1,"Z2",Z2.med(),0.001,0,0);
  
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
      minu.DefineParameter(0,"M",M[ijack_fit],0.001,0,0);
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
  cout<<"M: "<<M<<", ch2: "<<ch2<<endl;
  
  if(path!=NULL) write_constant_fit_plot(path,ecorr,M,tmin,tmax);
}


//fit the mass and the matrix element in SS and SL combo
double *c_two_pts_SL_fit[2],*e_two_pts_SL_fit[2];
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
      double teo=fun_two_pts_SL_fit(ZL,ZS,M,t);
      double diff=num-teo;
      double err=e_two_pts_SL_fit[0][t];
      double cont=sqr(diff/err);
      ch+=cont;
      if(flag==3) cout<<"SL, t="<<t<<", diff=("<<num<<"-"<<teo<<")="<<diff<<" err="<<err<<" cont="<<cont<<endl;
    }
  
  for(int t=tmin_two_pts_SL_fit[1];t<=min(tmax_two_pts_SL_fit[1],TH_two_pts_SL_fit);t++)
    {
      double diff=c_two_pts_SL_fit[1][t]-fun_two_pts_SL_fit(ZS,ZS,M,t);
      double err=e_two_pts_SL_fit[1][t];
      double cont=sqr(diff/err);
      ch+=cont;
      if(flag==3) cout<<"SS, t="<<t<<", diff="<<diff<<" err="<<err<<" cont="<<cont<<endl;
    }
}

void two_pts_SL_fit(jack &M,jack &ZL,jack &ZS,jvec corrSL,jvec corrSS,int tminL,int tmaxL,int tminS,int tmaxS,const char *path1=NULL,const char *path2=NULL)
{
  jvec ecorrSL=effective_mass(corrSL);
  jvec ecorrSS=effective_mass(corrSS);
  
  jack ML=constant_fit(ecorrSL,tminL,tmaxL,NULL);
  //jack MS=constant_fit(ecorrSS,tminS,tmaxS,NULL);
  M=ML;//jack_weighted_average(ML,MS);
  jvec tempSL(corrSS.nel,corrSS.njack),tempSS(corrSS.nel,corrSS.njack);
  int TH=tempSS.nel-1;
  for(int t=0;t<=TH;t++)
    {
      tempSL[t]=corrSL[t]/exp(-M*TH)/cosh(M*(TH-t))*M;
      tempSS[t]=corrSS[t]/exp(-M*TH)/cosh(M*(TH-t))*M;
    }
  
  ZS=sqrt(constant_fit(tempSS,tminS,tmaxS,NULL));
  ZL=constant_fit(tempSL,tminL,tmaxL,NULL)/ZS;
  
  TMinuit minu;
  minu.SetPrintLevel(-1);
  minu.SetFCN(ch2_two_pts_SL_fit);
  minu.DefineParameter(1,"ZL",ZL.med(),0.001,0,0);
  minu.DefineParameter(2,"ZS",ZS.med(),0.001,0,0);
  
  int njack=ML.njack;
  c_two_pts_SL_fit[0]=new double[TH+1];
  c_two_pts_SL_fit[1]=new double[TH+1];
  e_two_pts_SL_fit[0]=new double[TH+1];
  e_two_pts_SL_fit[1]=new double[TH+1];
  
  TH_two_pts_SL_fit=TH;
  tmin_two_pts_SL_fit[0]=tminL;
  tmin_two_pts_SL_fit[1]=tminS;
  tmax_two_pts_SL_fit[0]=tmaxL;
  tmax_two_pts_SL_fit[1]=tmaxS;
  
  for(int iel=0;iel<=TH;iel++)
    {
      e_two_pts_SL_fit[0][iel]=corrSL[iel].err();
      e_two_pts_SL_fit[1][iel]=corrSS[iel].err();
    }
  
  for(int ijack_fit=0;ijack_fit<=njack;ijack_fit++)
    {
      minu.DefineParameter(0,"M",M[ijack_fit],0.001,0,0);
      minu.FixParameter(0);
      for(int iel=0;iel<=TH;iel++)
	{
	  c_two_pts_SL_fit[0][iel]=corrSL[iel][ijack_fit];
	  c_two_pts_SL_fit[1][iel]=corrSS[iel][ijack_fit];
	}
      minu.Migrad();
      double dum;
      minu.GetParameter(0,M.data[ijack_fit],dum);
      minu.GetParameter(1,ZL.data[ijack_fit],dum);
      minu.GetParameter(2,ZS.data[ijack_fit],dum);
    }
  
  double ch2,grad[3],par[3]={M[njack],ZL[njack],ZS[njack]};
  minu.Eval(3,grad,ch2,par,3);
  cout<<"ML: "<<smart_print(ML)<<", ch2: "<<ch2<<endl;
  
  if(path1!=NULL) write_constant_fit_plot(path1,ecorrSL,M,tminL,tmaxL);
  if(path2!=NULL) write_constant_fit_plot(path2,ecorrSS,M,tminS,tmaxS);
}

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
