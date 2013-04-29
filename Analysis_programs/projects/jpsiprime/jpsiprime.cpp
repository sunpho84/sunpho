#include "include.h"

const int njacks=16;
const int T=48;
const int L=24;
const int nlevls=2;
const int nlevls_sto=2;
jvec data[nlevls][nlevls];

template<class to> to two(to Z1,to Z2,to M,int t)
{return Z1*Z2*exp(-M*L)*cosh(M*(L-t))/M;}

jvec load(const char *path,int ism_so_lv,int ism_si_lv)
{
  jvec a(T,njacks);
  jvec b(T,njacks);
  
  //load the charged corrs
  a.load(path,0+2*(0+4*(ism_si_lv+nlevls_sto*ism_so_lv)));
  b.load(path,0+2*(3+4*(ism_si_lv+nlevls_sto*ism_so_lv)));
  
  return (a+b)/2;
}

jack find_tmin(int *tmin,jvec *buf)
{
 jvec Mcho(nlevls*nlevls,njacks);
 
 for(int ism_so=0;ism_so<nlevls;ism_so++)
    for(int ism_si=0;ism_si<nlevls;ism_si++)
      {
	int ttest=0;
	double ch;
	do
	  {
	    jvec corr=effective_mass(buf[ism_so*nlevls+ism_si],L);

	    Mcho[ism_so*nlevls+ism_si]=constant_fit(corr,ttest,L,combine("plots/fit_choose_%02d_%02d.xmg",ism_so,ism_si).c_str());
	    int n=-1;
	    ch=0;
	    for(int t=ttest;t<corr.nel-1;t+=2)
	      if(!isnan(corr[t].err()))
		{
		  double ch_cont=sqr((corr[t].med()-Mcho[ism_so*nlevls+ism_si].med())/corr[t].err());
		  ch+=ch_cont;
		  n++;
		  //cout<<"tmin: "<<ttest<<", t: "<<t<<", ch: "<<ch_cont<<", n: "<<n<<endl;
		}
	    ch/=n;
	    
	    if(ch>2) ttest++;
	  }
	while(ch>2);
	cout<<"Ism_so: "<<ism_so<<", Ism_si: "<<ism_si<<", Ch2: "<<ch<<", t: "<<ttest<<endl;
	
	tmin[ism_so*nlevls+ism_si]=ttest;
      }
  cout<<Mcho<<endl;
  
  return constant_fit(Mcho,0,nlevls*nlevls);
}

double ch2_multi_buf[nlevls*nlevls][L+1];
double ch2_multi_err[nlevls*nlevls][L+1];
int *ch2_multi_tmin;
int ch2_multi_tmax[nlevls*nlevls];

void ch2_multi_lev(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  for(int ism_so=0;ism_so<nlevls;ism_so++)
    for(int ism_si=0;ism_si<nlevls;ism_si++)
      {
	int i=ism_so*nlevls+ism_si;
	for(int t=ch2_multi_tmin[i];t<ch2_multi_tmax[i];t++)
	  {
	    double yteo=two(p[ism_so],p[ism_si],p[nlevls],t);
	    double yspe=ch2_multi_buf[i][t];
	    double yerr=ch2_multi_err[i][t];
	    double ch_contr=sqr((yteo-yspe)/yerr);
	    ch+=ch_contr;
	    
	    //cout<<t<<" "<<ch_contr<<"=(("<<yteo<<"-"<<yspe<<")/"<<yerr<<")^2"<<endl;
	  }
      }
}

void fit_single_stat(jvec &Z,jack &M,jvec *buf,int *tmin)
{
  TMinuit minu;
  minu.SetFCN(ch2_multi_lev);
  
  //link tmin
  ch2_multi_tmin=tmin;
  for(int i=0;i<nlevls*nlevls;i++) ch2_multi_tmax[i]=buf[i].nel;
  
  //prepare Z estimates
  for(int ilev=0;ilev<nlevls;ilev++)
    {
      jvec test=buf[ilev*nlevls+ilev];
      for(int ijack_fit=0;ijack_fit<=njacks;ijack_fit++)
	for(int t=0;t<test.nel;t++)
	  test[t].data[ijack_fit]/=two(1.0,1.0,M[ijack_fit],t);

      Z[ilev]=sqrt(constant_fit(test,0,tmin[ilev*nlevls+ilev]));
    }
  
  //prepare the data error
  for(int ism_so=0;ism_so<nlevls;ism_so++)
    for(int ism_si=0;ism_si<nlevls;ism_si++)
      {
	int i=ism_so*nlevls+ism_si;
	for(int t=0;t<buf[i].nel;t++)
	  ch2_multi_err[i][t]=buf[i][t].err();
      }
  
  minu.SetPrintLevel(1);
  for(int ijack_fit=0;ijack_fit<=njacks;ijack_fit++)
    {
      //define parameters
      for(int ilevl=0;ilevl<nlevls;ilevl++)
	if(!isnan(Z[ilevl][ijack_fit]))
	  minu.DefineParameter(ilevl,combine("Z_%02d",ilevl).c_str(),Z[ilevl][ijack_fit],Z[ilevl].err(),0,0);
	else
	  minu.DefineParameter(ilevl,combine("Z_%02d",ilevl).c_str(),0.001,0.001,0,0);
      minu.DefineParameter( nlevls,"M",M[ijack_fit],M.err(),0,0);
      
      //copy data
      for(int ism_so=0;ism_so<nlevls;ism_so++)
	for(int ism_si=0;ism_si<nlevls;ism_si++)
	  {
	    int i=ism_so*nlevls+ism_si;
	    for(int t=0;t<buf[i].nel;t++)
	      ch2_multi_buf[i][t]=buf[i][t][ijack_fit];
	  }
      
      //minimze
      minu.Migrad();
      minu.SetPrintLevel(-1);
      
      //get pars back
      double dum;
      for(int ilevl=0;ilevl<nlevls;ilevl++)
	minu.GetParameter(ilevl,Z[ilevl].data[ijack_fit],dum);
      minu.GetParameter(nlevls,M.data[ijack_fit],dum);
    }
  
  //remove ground state
  for(int ism_so=0;ism_so<nlevls;ism_so++)
    for(int ism_si=0;ism_si<nlevls;ism_si++)
      {
	int i=ism_so*nlevls+ism_si;
	buf[i]=buf[i].subset(0,tmin[i]);
	for(int t=0;t<tmin[i];t++)
	  buf[i][t]-=two(Z[ism_so],Z[ism_si],M,t);
      }
}

void estimate_M_and_Z(jack &M0,jack &M1,jvec &Z0,jvec &Z1,int *tmin1,jvec *buf)
{
  //first of all finds for all the combo the plateaux region for ground state
  int tmin0[nlevls*nlevls];
  M0=find_tmin(tmin0,buf);
  cout<<"tmin first state: "<<tmin0[0]<<endl;
  //fit all
  fit_single_stat(Z0,M0,buf,tmin0);
  
  //redo
  M1=find_tmin(tmin1,buf);
  cout<<"tmin second state: "<<tmin1[0]<<endl;
  //redo again
  jvec Zfit1(nlevls,njacks);
  fit_single_stat(Z1,M1,buf,tmin1);
}

double ch2_two_buf[nlevls*nlevls][L+1];
double ch2_two_err[nlevls*nlevls][L+1];
int *ch2_two_tmin;
double ch2_two_ext;

void ch2_two_states_multi_lev(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  
  for(int ism_so=0;ism_so<nlevls;ism_so++)
    for(int ism_si=0;ism_si<nlevls;ism_si++)
      {
        int i=ism_so*nlevls+ism_si;
	for(int t=ch2_two_tmin[i];t<=L;t++)
	  {
	    double yspe=ch2_two_buf[i][t];
	    double yteo=two(p[ism_so],p[ism_si],p[nlevls],t)+two(p[nlevls+1+ism_so],p[nlevls+1+ism_si],p[nlevls+1+nlevls],t);
	    double yerr=ch2_two_err[i][t];
	    
            double ch_contr=sqr((yteo-yspe)/yerr);
            ch+=ch_contr;
	    
	    //cout<<t<<" "<<ch_contr<<"=(("<<yteo<<"-"<<yspe<<")/"<<yerr<<")^2"<<endl;	  
	  }
      }
  ch2_two_ext=ch;
}

void fit_M_and_Z(jvec &Z0,jvec &Z1,jack &M0,jack &M1,jvec *buf,int *tmin)
{
  TMinuit minu;
  minu.SetFCN(ch2_two_states_multi_lev);
  
  //link tmin
  ch2_two_tmin=tmin;
  
  //prepare the data error
  for(int ism_so=0;ism_so<nlevls;ism_so++)
    for(int ism_si=0;ism_si<nlevls;ism_si++)
      {
	int i=ism_so*nlevls+ism_si;
	for(int t=0;t<buf[i].nel;t++)
	  ch2_two_err[i][t]=buf[i][t].err();
      }
  
  minu.SetPrintLevel(1);
  for(int ijack_fit=0;ijack_fit<=njacks;ijack_fit++)
    {
      //define parameters for first state
      for(int ilevl=0;ilevl<nlevls;ilevl++)
	  minu.DefineParameter(ilevl,combine("Z0_%02d",ilevl).c_str(),Z0[ilevl][ijack_fit],Z0[ilevl].err(),0,0);
      minu.DefineParameter( nlevls,"M0",M0[ijack_fit],M0.err(),0,0);
      //define parameters for excited state
      for(int ilevl=0;ilevl<nlevls;ilevl++)
	  minu.DefineParameter(ilevl+1+nlevls,combine("Z1_%02d",ilevl).c_str(),Z1[ilevl][ijack_fit],Z1[ilevl].err(),0,0);
      minu.DefineParameter(nlevls+1+nlevls,"M1",M1[ijack_fit],M1.err(),0,0);
      
      //copy data
      for(int ism_so=0;ism_so<nlevls;ism_so++)
	for(int ism_si=0;ism_si<nlevls;ism_si++)
	  {
	    int i=ism_so*nlevls+ism_si;
	    for(int t=0;t<buf[i].nel;t++)
	      ch2_two_buf[i][t]=buf[i][t][ijack_fit];
	  }
      
      //minimze
      minu.Migrad();
      minu.SetPrintLevel(-1);
      
      //get pars back
      double dum;
      for(int ilevl=0;ilevl<nlevls;ilevl++)
	minu.GetParameter(ilevl,Z0[ilevl].data[ijack_fit],dum);
      minu.GetParameter(nlevls,M0.data[ijack_fit],dum);
      for(int ilevl=0;ilevl<nlevls;ilevl++)
	minu.GetParameter(ilevl+nlevls+1,Z1[ilevl].data[ijack_fit],dum);
      minu.GetParameter(nlevls+1+nlevls,M1.data[ijack_fit],dum);
    }
  
  cout<<"Chi2: "<<ch2_two_ext/(nlevls*nlevls-2*nlevls-2)<<endl;
}

int main()
{
  //buffer where to copy the data
  jvec buf[nlevls*nlevls];
  
  //load data
  for(int ism_so=0;ism_so<nlevls;ism_so++)
    for(int ism_si=0;ism_si<nlevls;ism_si++)
      buf[ism_so*nlevls+ism_si]=data[ism_so][ism_si]=load("correlators/2pts_P5P5",ism_so,ism_si).simmetrized(1);

  //perform the estimate of Z,M and tmin
  int tmin[nlevls*nlevls];
  jvec Z0(nlevls,njacks),Z1(nlevls,njacks);
  jack M0(njacks),M1(njacks);
  estimate_M_and_Z(M0,M1,Z0,Z1,tmin,buf);
  cout<<"M0: "<<M0<<endl;
  cout<<"M1: "<<M1<<endl;

  for(int i=0;i<nlevls*nlevls;i++)
    {
      if(tmin[i]<7) tmin[i]=7;
      cout<<i<<" "<<tmin[i]<<endl;
    }
  
  //now make the true fit
  fit_M_and_Z(Z0,Z1,M0,M1,(jvec*)data,tmin);
  
  cout<<"M0: "<<M0<<endl;
  cout<<"M1: "<<M1<<endl;
  cout<<"M1/M0-1: "<<M1/M0-1<<endl;

  return 0;
}
