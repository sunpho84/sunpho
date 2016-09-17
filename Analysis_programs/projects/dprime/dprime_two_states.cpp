#include <math.h>
#include "include.h"
#include "eff_bis.cpp"

const int deb=0;

int clov_tm;
const int njacks=16;
int T,TH;
int L;
int nml,iml,imc;
int **tmin,**tmax;
int nst;
double ml,mc;
int nsm,*sm;
int **use,**plot;

double **corr_two_pts_fit,**err_two_pts_fit;
int *tmin_fit;
int *tmax_fit;
int *use_fit;
int *sign,*ri;

void read_input(const char *path)
{
  FILE *fin=open_file(path,"r");
  
  //read T
  read_formatted_from_file_expecting((char*)&T,fin,"%d","T");
  L=TH=T/2;
  
  //read tm or wclov
  read_formatted_from_file_expecting((char*)&clov_tm,fin,"%d","clov_tm");
  
  //read nml
  read_formatted_from_file_expecting((char*)&nml,fin,"%d","nml");
  read_formatted_from_file_expecting((char*)&iml,fin,"%d","iml");
  read_formatted_from_file_expecting((char*)&imc,fin,"%d","imc");
  
  //read ml,mc
  read_formatted_from_file_expecting((char*)&ml,fin,"%lg","ml");
  read_formatted_from_file_expecting((char*)&mc,fin,"%lg","mc");
  
  //read nsm
  read_formatted_from_file_expecting((char*)&nsm,fin,"%d","nsm");
  sm=(int*)malloc(sizeof(int)*nsm);
  for(int ism=0;ism<nsm;ism++) read_formatted_from_file((char*)(sm+ism),fin,"%d","sm");
  
  //read sign and parity
  sign=(int*)malloc(sizeof(int)*nsm*nsm);
  ri=(int*)malloc(sizeof(int)*nsm*nsm);
  expect_string_from_file(fin,"sign_ri");
  for(int ism=0;ism<nsm;ism++)
    for(int jsm=0;jsm<nsm;jsm++)
      {
	read_formatted_from_file((char*)(sign+ism*nsm+jsm),fin,"%d","sign");
	read_formatted_from_file((char*)(ri+ism*nsm+jsm),fin,"%d","ri");
      }
  
  //read nst
  read_formatted_from_file_expecting((char*)&nst,fin,"%d","nst");
  tmin=(int**)malloc(sizeof(int*)*nst);
  tmax=(int**)malloc(sizeof(int*)*nst);
  for(int ist=0;ist<nst;ist++)
    {
      tmin[ist]=(int*)malloc(sizeof(int)*nsm);
      tmax[ist]=(int*)malloc(sizeof(int)*nsm);
      
      read_formatted_from_file_expecting((char*)tmin[ist],fin,"%d",combine("tint_st_%d",ist).c_str());
      read_formatted_from_file((char*)tmax[ist],fin,"%d","tint");
      
      for(int ism=1;ism<nsm;ism++)
	{
	  read_formatted_from_file((char*)(tmin[ist]+ism),fin,"%d","tint");
	  read_formatted_from_file((char*)(tmax[ist]+ism),fin,"%d","tint");
	}
    }
  
  //decide whether to use or not the level, and to plot it
  use=(int**)malloc(sizeof(int*)*nst);
  plot=(int**)malloc(sizeof(int*)*nst);
  for(int ist=0;ist<nst;ist++)
    {
      use[ist]=(int*)malloc(sizeof(int)*nsm);
      plot[ist]=(int*)malloc(sizeof(int)*nsm);
      for(int ism=0;ism<nsm;ism++)
	{
	  use[ist][ism]=(tmin[ist][ism]!=tmax[ist][ism]);
	  plot[ist][ism]=(ist==0||use[ist-1][ism]);
	}
    }
  
  fclose(fin);
}

jvec load(const char *path,int im1,int im2,int ri)
{
  if(clov_tm)
    {
      int r=0;
      jvec a(T,njacks),b(T,njacks);
      
      a.load(path,ri+2*( r+2*(im1+nml*( r+2*im2))));
      b.load(path,ri+2*(!r+2*(im1+nml*(!r+2*im2))));
      
      if(ri==0) return ((a+b)/2).simmetrized(1);
      else      return ((b-a)/2).simmetrized(1);
    }
  else
    {
      jvec a(T,njacks);
      
      a.load(path,ri+2*(im1+nml*im2));
      
      return a.simmetrized(1);
    }
}

template <class T> T fun_migrad_A_fit(T Za,T Zb,T M,double t)
{return Za*Zb*exp(-M*TH)*cosh(M*(TH-t))/M;}

template <class T> T fun_migrad_B_fit(T *Za,T *Zb,T *M,double t)
{
  T out=Za[0]*0;
  
  for(int ist=0;ist<nst;ist++) out+=fun_migrad_A_fit(Za[ist],Zb[ist],M[ist],t);
  
  return out;
}

template <class T> T fD(T Z,T M)
{return (ml+mc)*Z/M/sinh(M);}

void chi2_stage_A(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  double M=p[nsm];
  double *Z=p;
  
  for(int ism=0;ism<nsm;ism++)
    for(int jsm=0;jsm<nsm;jsm++)
      if(use_fit[ism]&&use_fit[jsm])
	for(int t=tmin_fit[ism*nsm+jsm];t<=tmax_fit[ism*nsm+jsm];t++)
	  {
	    double num=corr_two_pts_fit[ism*nsm+jsm][t];
	    double teo=fun_migrad_A_fit(Z[ism],Z[jsm],M,t);
	    double diff=num-teo;
	    double err=err_two_pts_fit[ism][t];
	    double cont=sqr(diff/err);
	    ch+=cont;
	  }
}

void stage_A_fit(jack *M,jack *Z,jvec *corr,int ist)
{
  //minimizator
  TMinuit minu;
  minu.SetPrintLevel(-1);
  minu.SetFCN(chi2_stage_A);

  //decide the intervals
  use_fit=use[ist];
  tmin_fit=new int[nsm*nsm];
  tmax_fit=new int[nsm*nsm];
  for(int ism=0;ism<nsm;ism++)
    for(int jsm=0;jsm<nsm;jsm++)
      {
	tmin_fit[ism*nsm+jsm]=min(tmin[ist][ism],tmin[ist][jsm]);
	tmax_fit[ism*nsm+jsm]=min(tmax[ist][ism],tmax[ist][jsm]);
      }
  
  //copy the data
  corr_two_pts_fit=new double*[nsm*nsm];
  err_two_pts_fit=new double*[nsm*nsm];
  for(int ism=0;ism<nsm;ism++)
    for(int jsm=0;jsm<nsm;jsm++)
      {
	err_two_pts_fit [ism*nsm+jsm]=new double[TH+1];
	corr_two_pts_fit[ism*nsm+jsm]=new double[TH+1];

	for(int t=0;t<=TH;t++) err_two_pts_fit[ism*nsm+jsm][t]=corr[ism*nsm+jsm][t].err();
    }
  
  //define the parameters
  for(int ism=0;ism<nsm;ism++)
    {
      if(use[ist][ism])
	{
	  two_pts_fit(M[ist],Z[ist*nsm+ism],corr[ism*nsm+ism],tmin[ist][ism],tmax[ist][ism],combine("plots/pre_fit_st_%d_sm_%02d.xmg",ist,sm[ism]).c_str());
	  Z[ist*nsm+ism]=sqrt(Z[ist*nsm+ism]);
	  minu.DefineParameter(ism,combine("Z_%d",ism).c_str(),Z[ist*nsm+ism].med(),Z[ist*nsm+ism].err(),0,0);
	}
      else
	{
	  minu.DefineParameter(ism,combine("Z_%d",ism).c_str(),0,0,0,0);
	  minu.FixParameter(ism);
	}
    }
  minu.DefineParameter(nsm,"M",M[ist].med(),M[ist].err(),0,0);
  
  for(int ijack_fit=0;ijack_fit<=njacks;ijack_fit++)
    {
      //copy jack
      for(int ism=0;ism<nsm;ism++)
	for(int jsm=0;jsm<nsm;jsm++)
	  for(int t=0;t<=TH;t++)
	    corr_two_pts_fit[ism*nsm+jsm][t]=corr[ism*nsm+jsm][t][ijack_fit];
      
      //minimize
      minu.Migrad();
      
      //get back the parameters
      double dum;
      minu.GetParameter(nsm,M[ist].data[ijack_fit],dum);
      for(int ism=0;ism<nsm;ism++)
	if(use[ist][ism])
	  minu.GetParameter(ism,Z[ist*nsm+ism].data[ijack_fit],dum);
    }
  
  //get pars back
  double par[nsm+1];
  for(int ism=0;ism<=nsm;ism++)
    {
      double p,dum;
      minu.GetParameter(ism,p,dum);
      par[ism]=p;
    }
  
  //compute chi2
  double chi2;
  minu.Eval(nsm+1,NULL,chi2,par,2);
  cout<<"M: "<<smart_print(M[ist])<<", Z: "<<smart_print(Z[ist*nsm+0])<<", ch2: "<<chi2<<endl;
  
  //write the plots
  for(int ism=0;ism<nsm;ism++)
    for(int jsm=0;jsm<nsm;jsm++)
      if(plot[ism]&&plot[jsm])
	write_constant_fit_plot(combine("plots/fit_A_st_%d_sm_%02d_%02d.xmg",ist,sm[ism],sm[jsm]).c_str(),
				effective_mass(corr[ism*nsm+jsm],TH),M[ist],max(tmin[ist][ism],tmin[ist][jsm]),min(tmax[ist][ism],tmax[ist][jsm]));
}

void chi2_stage_B(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  double M[nst],Z[nsm][nst];
  
  for(int ist=0;ist<nst;ist++)
    {
      M[ist]=p[ist*(nsm+1)+nsm];
      for(int ism=0;ism<nsm;ism++) Z[ism][ist]=p[ist*(nsm+1)+ism];
    }
  
  for(int ism=0;ism<nsm;ism++)
    for(int jsm=0;jsm<nsm;jsm++)
      for(int t=tmin_fit[ism*nsm+jsm];t<=tmax_fit[ism*nsm+jsm];t++)
	{
	  double num=corr_two_pts_fit[ism*nsm+jsm][t];
	  double teo=fun_migrad_B_fit(Z[ism],Z[jsm],M,t);
	  double diff=num-teo;
	  double err=err_two_pts_fit[ism][t];
	  double cont=sqr(diff/err);
	  ch+=cont;
	}
}

void stage_B_fit(jack *M,jack *Z,jvec *corr)
{
  //minimizator
  TMinuit minu;
  minu.SetPrintLevel(-1);
  minu.SetFCN(chi2_stage_B);

  //decide the intervals
  tmin_fit=new int[nsm*nsm];
  tmax_fit=new int[nsm*nsm];
  for(int ism=0;ism<nsm;ism++)
    for(int jsm=0;jsm<nsm;jsm++)
      {
	tmin_fit[ism*nsm+jsm]=min(tmin[0][ism],tmin[0][jsm]);
	tmax_fit[ism*nsm+jsm]=max(tmax[0][ism],tmax[0][jsm]);
	for(int ist=1;ist<nst;ist++)
	  {
	    tmin_fit[ism*nsm+jsm]=min(tmin_fit[ism*nsm+jsm],min(tmin[ist][ism],tmin[ist][jsm]));
	    tmax_fit[ism*nsm+jsm]=max(tmax_fit[ism*nsm+jsm],max(tmax[ist][ism],tmax[ist][jsm]));
	  }
      }
  
  //copy the data
  corr_two_pts_fit=new double*[nsm*nsm];
  err_two_pts_fit=new double*[nsm*nsm];
  for(int ism=0;ism<nsm;ism++)
    for(int jsm=0;jsm<nsm;jsm++)
      {
	err_two_pts_fit [ism*nsm+jsm]=new double[TH+1];
	corr_two_pts_fit[ism*nsm+jsm]=new double[TH+1];
	for(int t=0;t<=TH;t++) err_two_pts_fit[ism*nsm+jsm][t]=corr[ism*nsm+jsm][t].err();
    }
  
  //define the parameters
  for(int ist=0;ist<nst;ist++)
    {
      for(int ism=0;ism<nsm;ism++)
	{
	  if(use[ist][ism]) minu.DefineParameter(ist*(nsm+1)+ism,combine("Z_st_%d_sm_%02d",ist,sm[ism]).c_str(),Z[ist*nsm+ism].med(),Z[ist*nsm+ism].err(),0,0);
	  else
	    {
	      minu.DefineParameter(ist*(nsm+1)+ism,combine("Z_st_%d_sm_%02d",ist,sm[ism]).c_str(),0,0,0,0);
	      minu.FixParameter(ist*(nsm+1)+ism);
	    }
	}
      minu.DefineParameter(ist*(nsm+1)+nsm,combine("M_st_%d",ist).c_str(),M[ist].med(),M[ist].err(),0,0);
    }
  
  for(int ijack_fit=0;ijack_fit<=njacks;ijack_fit++)
    {
      //copy jack
      for(int ism=0;ism<nsm;ism++)
	for(int jsm=0;jsm<nsm;jsm++)
	  for(int t=0;t<=TH;t++)
	    corr_two_pts_fit[ism*nsm+jsm][t]=corr[ism*nsm+jsm][t][ijack_fit];
      
      //minimize
      minu.Migrad();
      
      //get back the parameters
      double dum;
      for(int ist=0;ist<nst;ist++)
	{
	  minu.GetParameter(ist*(nsm+1)+nsm,M[ist].data[ijack_fit],dum);
	  for(int ism=0;ism<nsm;ism++)
	    if(use[ist][ism])
	      minu.GetParameter(ist*(nsm+1)+ism,Z[ist*nsm+ism].data[ijack_fit],dum);
	}
    }
  
  //get pars back
  int npar=(nsm+1)*nst;
  double par[npar];
  for(int ipar=0;ipar<npar;ipar++)
    {
      double p,dum;
      minu.GetParameter(ipar,p,dum);
      par[ipar]=p;
    }
  
  //compute chi2
  double chi2;
  minu.Eval(npar,NULL,chi2,par,2);
  cout<<"ch2: "<<chi2<<endl;
  for(int ist=0;ist<nst;ist++)
    cout<<"M: "<<smart_print(M[ist])<<", Z: "<<smart_print(Z[ist*nsm+0])<<endl;
  
  //write the plots
  for(int ist=0;ist<=nst;ist++)
    for(int ism=0;ism<nsm;ism++)
      for(int jsm=0;jsm<nsm;jsm++)
	if(plot[ism]&&plot[jsm])
	  {
	    //subtracted corr
	    jvec temp_corr=corr[ism*nsm+jsm];	    
	    for(int jst=0;jst<ist;jst++)
	      for(int t=0;t<=TH;t++)
		temp_corr[t]-=fun_migrad_A_fit(Z[jst*nsm+ism],Z[jst*nsm+jsm],M[jst],t);
	    
	    //summ of exponentials
	    int npoints=tmax_fit[ism*nsm+jsm]-tmin_fit[ism*nsm+jsm];
	    double dx=((double)tmax_fit[ism*nsm+jsm]-tmin_fit[ism*nsm+jsm])/(npoints-1);
	    double x[npoints];
	    jvec temp_fun(npoints,njacks);
	    for(int i=0;i<npoints;i++)
	      {
		double x0=x[i]=tmin_fit[ism*nsm+jsm]+i*dx;
		double x1=x0+dx;
		jack y0(njacks),y1(njacks);
		y0=y1=0;
		for(int jst=ist;jst<nst;jst++)
		  {
		    y0+=fun_migrad_A_fit(Z[jst*nsm+ism],Z[jst*nsm+jsm],M[jst],x0);
		    y1+=fun_migrad_A_fit(Z[jst*nsm+ism],Z[jst*nsm+jsm],M[jst],x1);
		  }
		temp_fun[i]=y0;//effective_masx(x,;//-log(y1/y0)/dx;
	      }
	    temp_fun=effective_mass(x,temp_fun,TH);
	    
	    ofstream out(combine("plots/fit_B_st_%d_sm_%02d_%02d.xmg",ist,sm[ism],sm[jsm]).c_str());
	    //error of the line
	    out<<"@s0 line type 1"<<endl;      
	    out<<"@s0 line color 7"<<endl;
	    out<<"@s0 fill color 7"<<endl;
	    out<<"@s0 fill type 1"<<endl;
	    out<<"@type xy"<<endl;
	    for(int i=0;i<npoints-1;i++) out<<x[i]<<" "<<temp_fun[i].med()+temp_fun[i].err()<<endl;
	    for(int i=npoints-2;i>=0;i--) out<<x[i]<<" "<<temp_fun[i].med()-temp_fun[i].err()<<endl;
	    out<<"&"<<endl;
	    //central line
	    out<<"@s1 line color 1"<<endl;
	    out<<"@type xy"<<endl;      
	    for(int i=0;i<npoints-1;i++) out<<x[i]<<" "<<temp_fun[i].med()<<endl;
	    //plot the original data with error  
	    out<<"&"<<endl;
	    out<<"@type xydy"<<endl;      
	    out<<"@s2 line type 0"<<endl;      
	    out<<"@s2 symbol color 1"<<endl;
	    out<<"@s2 errorbar color 1"<<endl;
	    out<<"@s2 symbol 1"<<endl;
	    out<<effective_mass(temp_corr);
	    out<<"&"<<endl;
	    out.close();
	  }
}

int main(int narg,char **arg)
{
  if(narg<2) crash("use %s input",arg[0]);
  read_input(arg[1]);
  
  jvec corr[nst][nsm*nsm];
  jack M[nst],Z[nst*nsm];
  
  cout<<"Stage A"<<endl;
  for(int ist=0;ist<nst;ist++)
    {
      for(int ism=0;ism<nsm;ism++)
	for(int jsm=0;jsm<nsm;jsm++)
	  if(ist==0)
	    {
	      corr[ist][ism*nsm+jsm]=sign[ism*nsm+jsm]*load(combine("corrs/2pts_P5P5_%02d_%02d",sm[ism],sm[jsm]).c_str(),imc,iml,ri[ism*nsm+jsm]);
	      corr[ist][ism*nsm+jsm].print_to_file(combine("plots/raw_data_%02d_%02d.xmg",sm[ism],sm[jsm]).c_str());
	    }
	  else
	    if(use[ist-1][ism]&&use[ist-1][jsm])
	      {
		corr[ist][ism*nsm+jsm]=corr[ist-1][ism*nsm+jsm];
		for(int t=0;t<=TH;t++)
		  corr[ist][ism*nsm+jsm][t]-=
		    fun_migrad_A_fit(Z[(ist-1)*nsm+ism],Z[nsm*(ist-1)+jsm],M[ist-1],t);
	      }
      
      stage_A_fit(M,Z,corr[ist],ist);
    }
  
  cout<<"MD'/MD: "<<smart_print(M[1]/M[0])<<endl;
  cout<<"fD'/fD: "<<fD(Z[1*nsm+0],M[1])/fD(Z[0*nsm+0],M[0])<<endl;
  
  cout<<endl<<"Stage B"<<endl;
  stage_B_fit(M,Z,corr[0]);
  cout<<"MD'/MD: "<<smart_print(M[1]/M[0])<<endl;
  cout<<"fD'/fD: "<<fD(Z[1*nsm+0],M[1])/fD(Z[0*nsm+0],M[0])<<endl;
    
  return 0;
}
