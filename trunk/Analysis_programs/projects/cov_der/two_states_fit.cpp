#include <math.h>
#include "include.h"
#include "../jpsiprime/eff_bis.cpp"
#include "../nf2/common_pars.cpp"

const int deb=0;

int itheta;
int njacks;
int only_charged,tm_run;
int T,TH,L;
int ibeta;
int **tmin,**tmax;
int nst;
double mc;
int nsm;
char sm[100][10];
int **use,**plot;

char infile[100],plotpath[100],current[100];
double **corr_two_pts_fit,**err_two_pts_fit;
int *tmin_fit;
int *tmax_fit;
int *use_fit;
int *sign,*parity;

char outfile[100];

void read_input(const char *path)
{
  FILE *fin=open_file(path,"r");
  
  read_formatted_from_file_expecting(infile,fin,"%s","infile");
  read_formatted_from_file_expecting(plotpath,fin,"%s","plotpath");
  read_formatted_from_file_expecting(outfile,fin,"%s","outfile");  
  
  //read itheta
  read_formatted_from_file_expecting((char*)&itheta,fin,"%d","itheta");
  
  //read T
  read_formatted_from_file_expecting((char*)&T,fin,"%d","T");
  L=TH=T/2;
  
  //read ibeta and correlator name
  read_formatted_from_file_expecting((char*)&ibeta,fin,"%d","Beta"); 
  read_formatted_from_file_expecting(current,fin,"%s","Current");
  
  read_formatted_from_file_expecting((char*)&njacks,fin,"%d","NJacks");  
  read_formatted_from_file_expecting((char*)&only_charged,fin,"%d","OnlyCharged");
  read_formatted_from_file_expecting((char*)&tm_run,fin,"%d","TwistedMass");
  
  //read mc
  read_formatted_from_file_expecting((char*)&mc,fin,"%lg","mc");
  
  //read nsm
  read_formatted_from_file_expecting((char*)&nsm,fin,"%d","nsm");
  for(int ism=0;ism<nsm;ism++) read_formatted_from_file(sm[ism],fin,"%s","sm");
  
  //read sign and ri
  sign=(int*)malloc(sizeof(int)*nsm*nsm);
  parity=(int*)malloc(sizeof(int)*nsm*nsm);
  expect_string_from_file(fin,"sign");
  for(int ism=0;ism<nsm;ism++)
    for(int jsm=0;jsm<nsm;jsm++)
      read_formatted_from_file((char*)(sign+ism*nsm+jsm),fin,"%d","sign");
  
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
	  use[ist][ism]=(tmin[ist][ism]<tmax[ist][ism]);
	  plot[ist][ism]=(ist==0||use[ist-1][ism]);
	}
    }
  
  fclose(fin);
}

template <class T> T fun_migrad_A_fit(T Za,T Zb,T M,double t,int par)
{
  if(par==1) return Za*Zb*exp(-M*TH)*cosh(M*(TH-t))/M;
  else return Za*Zb*exp(-M*TH)*sinh(M*(TH-t))/M;
}

template <class T> T fun_migrad_B_fit(T *Za,T *Zb,T *M,double t,int par)
{
  T out=Za[0]*0;
  
  for(int ist=0;ist<nst;ist++) out+=fun_migrad_A_fit(Za[ist],Zb[ist],M[ist],t,par);
  
  return out;
}

template <class T> T f(T Z,T M)
{
  if(string(current)==string("P")) return 2*mc*Z/M/sinh(M);
  else
    if(string(current)==string("V")) return Z/M*Za_med[ibeta];
    else
      {
	crash("Unknown current: %s",current);
	return M;
      }
}

void chi2_stage_A(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  double M=p[nsm];
  double *Z=p;
  
  for(int ism=0;ism<nsm;ism++)
    for(int jsm=0;jsm<nsm;jsm++)
      if(sign[ism*nsm+jsm]!=0)
	if(use_fit[ism]&&use_fit[jsm])
	  for(int t=tmin_fit[ism*nsm+jsm];t<=tmax_fit[ism*nsm+jsm];t++)
	    {
	      double num=corr_two_pts_fit[ism*nsm+jsm][t];
	      double teo=fun_migrad_A_fit(Z[ism],Z[jsm],M,t,parity[ism*nsm+jsm]);
	      double diff=num-teo;
	      double err=err_two_pts_fit[ism*nsm+jsm][t];
	      double cont=sqr(diff/err);
	      ch+=cont;
	      
	      if(flag==3) cout<<ism<<" "<<jsm<<" "<<t<<" "<<parity[ism*nsm+jsm]<<": [("<<num<<"-"<<teo<<")/"<<err<<"]^2="<<cont<<endl;
	    }
}

void stage_A_fit(jvec &M,jvec &Z,jvec *corr,int ist)
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
	tmin_fit[ism*nsm+jsm]=(tmin[ist][ism]+tmin[ist][jsm])/2.0;
	tmax_fit[ism*nsm+jsm]=(tmax[ist][ism]+tmax[ist][jsm])/2.0;
	cout<<"Using for combo "<<ism<<" "<<jsm<<" fit interval ["<<tmin_fit[ism*nsm+jsm]<<":"
	    <<tmax_fit[ism*nsm+jsm]<<"]"<<endl; 
      }
  
  //copy the data
  corr_two_pts_fit=new double*[nsm*nsm];
  err_two_pts_fit=new double*[nsm*nsm];
  for(int ism=0;ism<nsm;ism++)
    for(int jsm=0;jsm<nsm;jsm++)
      if(sign[ism*nsm+jsm]!=0)
	{
	  err_two_pts_fit [ism*nsm+jsm]=new double[TH+1];
	  corr_two_pts_fit[ism*nsm+jsm]=new double[TH+1];
	  
	  for(int t=0;t<=TH;t++) err_two_pts_fit[ism*nsm+jsm][t]=corr[ism*nsm+jsm][t].err();
	}
  
  cout<<"PREFIT "<<ist<<endl;
  
  jvec testM(nsm,njacks);
  //define the parameters of diagonal-defined operators
  for(int ism=0;ism<nsm;ism++)
    if(use[ist][ism] && sign[ism*nsm+ism]!=0)
      {
	cout<<" prefitting "<<ism<<" between "<<tmin_fit[ism*nsm+ism]<<" & "<<tmax_fit[ism*nsm+ism]<<endl;
	two_pts_fit(testM[ism],Z[ist*nsm+ism],corr[ism*nsm+ism],tmin_fit[ism*nsm+ism],tmax_fit[ism*nsm+ism],
		    combine("%s/pre_fit_st_%d_sm_%s.xmg",plotpath,ist,sm[ism]).c_str(),NULL,TH,parity[ism*nsm+ism]);
	Z[ist*nsm+ism]=sqrt(Z[ist*nsm+ism]);
	cout<<" guess for sm "<<ism<<" ("<<sm[ism]<<"): "<<Z[ist*nsm+ism]<<endl;
	minu.DefineParameter(ism,combine("Z_%d",ism).c_str(),Z[ist*nsm+ism].med(),Z[ist*nsm+ism].err(),0,0);
      }
  
  cout<<endl;
  
  //fix them using non-diagonal correlator
  for(int ism=0;ism<nsm;ism++)
    if(!(use[ist][ism] && sign[ism*nsm+ism]!=0)||isnan(Z[ist*nsm+ism].err()))
      {
	//cout<<use[ist][ism]<<" "<<(sign[ism*nsm+ism]!=0)<<" "<<isnan(Z[ist*nsm+ism].err())<<" "<<Z[ist*nsm+ism]<<" sm: "<<ism<<" "<<nst<<endl;
	//search a non-diagonal correlator
	int jsm=0;
	while(jsm<nsm && (sign[ism*nsm+jsm]==0 || ism==jsm)) jsm++;
	if(jsm==nsm) crash("operator %d is not fixed anyway",ism);
	two_pts_fit(testM[ism],Z[ist*nsm+ism],corr[ism*nsm+jsm],tmin_fit[ism*nsm+jsm],tmax_fit[ism*nsm+jsm],
		    combine("%s/pre_fit_st_%d_sm_%s.xmg",plotpath,ist,sm[ism]).c_str());
	cout<<" refixing "<<ism<<" using "<<jsm<<": "<<Z[ist*nsm+ism]<<" / "<<Z[ist*nsm+jsm]<<endl;
	Z[ist*nsm+ism]/=Z[ist*nsm+jsm];
	minu.DefineParameter(ism,combine("Z_%d",ism).c_str(),Z[ist*nsm+ism].med(),Z[ist*nsm+ism].err(),0,0);
      }
  
  //print the correlator for check
  for(int ism=0;ism<nsm;ism++)
    for(int jsm=0;jsm<nsm;jsm++)
      if(sign[ism*nsm+jsm]!=0)
	cout<<"ism "<<ism<<" jsm "<<jsm<<": "<<corr[ism*nsm+jsm][1]<<endl;
  
  cout<<"FIT"<<endl;
  
  //compute average M
  M[ist]=0;
  double wtot=0;
  for(int ism=0;ism<nsm;ism++)
    {
      double w=sqr(1/testM[ism].err());
      M[ist]+=testM[ism]*w;
      wtot+=w;
    }
  M[ist]/=wtot;
  
  minu.DefineParameter(nsm,"M",M[ist].med(),M[ist].err(),0,0);
  
  for(int ijack_fit=0;ijack_fit<=njacks;ijack_fit++)
    {
      //copy jack
      for(int ism=0;ism<nsm;ism++)
	for(int jsm=0;jsm<nsm;jsm++)
	  if(sign[ism*nsm+jsm]!=0)
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
  minu.Eval(nsm+1,NULL,chi2,par,3);
  cout<<"M: "<<smart_print(M[ist])<<", Z: "<<smart_print(Z[ist*nsm+0])<<", ch2: "<<chi2<<endl;
  
  //write the plots
  for(int ism=0;ism<nsm;ism++)
    for(int jsm=0;jsm<nsm;jsm++)
      if(sign[ism*nsm+jsm]!=0)	  
	if(plot[ist][ism]&&plot[ist][jsm])
	  write_constant_fit_plot(combine("%s/fit_A_st_%d_sm_%s_%s.xmg",plotpath,ist,sm[ism],sm[jsm]).c_str(),
				  effective_mass(corr[ism*nsm+jsm],TH,parity[ism*nsm+jsm]),M[ist],tmin_fit[ism*nsm+jsm],tmax_fit[ism*nsm+jsm]);
}

void chi2_stage_B(int &npar,double *fuf,double &ch,double *p,int flag)
{
  int npoints=-nst*(1+nsm); //M+nsm*Z
  ch=0;
  double M[nst],Z[nsm][nst];
  
  for(int ist=0;ist<nst;ist++)
    {
      M[ist]=p[ist*(nsm+1)+nsm];
      for(int ism=0;ism<nsm;ism++) Z[ism][ist]=p[ist*(nsm+1)+ism];
    }
  
  for(int ism=0;ism<nsm;ism++)
    for(int jsm=0;jsm<nsm;jsm++)
      if(sign[ism*nsm+jsm]!=0)
	for(int t=tmin_fit[ism*nsm+jsm];t<=tmax_fit[ism*nsm+jsm];t++)
	  {
	    double num=corr_two_pts_fit[ism*nsm+jsm][t];
	    double teo=fun_migrad_B_fit(Z[ism],Z[jsm],M,t,parity[ism*nsm+jsm]);
	    double diff=num-teo;
	    double err=err_two_pts_fit[ism*nsm+jsm][t];
	    double cont=sqr(diff/err);
	    
	    npoints++;
	    ch+=cont;
	  }
  
  ch/=npoints;
}

void stage_B_fit(jvec &M,jvec &Z,jvec *corr)
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
	tmin_fit[ism*nsm+jsm]=max(tmin[0][ism],tmin[0][jsm]);
	tmax_fit[ism*nsm+jsm]=max(tmax[0][ism],tmax[0][jsm]);
	for(int ist=1;ist<nst;ist++)
	  {
	    tmin_fit[ism*nsm+jsm]=min(tmin_fit[ism*nsm+jsm],(tmin[ist][ism]+tmin[ist][jsm])/2);
	    tmax_fit[ism*nsm+jsm]=max(tmax_fit[ism*nsm+jsm],max(tmax[ist][ism],tmax[ist][jsm]));
	  }
	
	cout<<"stage B fit intervals, ism="<<ism<<", jsm="<<jsm<<", ["<<tmin_fit[ism*nsm+jsm]<<";"<<tmax_fit[ism*nsm+jsm]<<"]"<<endl;
      }  

  corr_two_pts_fit=new double*[nsm*nsm];
  err_two_pts_fit=new double*[nsm*nsm];
  for(int ism=0;ism<nsm;ism++)
    for(int jsm=0;jsm<nsm;jsm++)
      if(sign[ism*nsm+jsm]!=0)
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
	  //if(use[ist][ism])
	    minu.DefineParameter(ist*(nsm+1)+ism,combine("Z_st_%d_sm_%s",ist,sm[ism]).c_str(),Z[ist*nsm+ism].med(),Z[ist*nsm+ism].err(),0,00);
	    //else
	    //{
	    //minu.DefineParameter(ist*(nsm+1)+ism,combine("Z_st_%d_sm_%s",ist,sm[ism]).c_str(),0,0,0,0);
	    //minu.FixParameter(ist*(nsm+1)+ism);
	    //}
	}
      minu.DefineParameter(ist*(nsm+1)+nsm,combine("M_st_%d",ist).c_str(),M[ist].med(),M[ist].err(),0,0);
    }
  
  for(int ijack_fit=0;ijack_fit<=njacks;ijack_fit++)
    {
      //copy jack
      for(int ism=0;ism<nsm;ism++)
	for(int jsm=0;jsm<nsm;jsm++)
	  if(sign[ism*nsm+jsm]!=0)
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
	    //if(use[ist][ism])
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
  
  ///////////////// write the plots //////////////
  
  //separate plots
  for(int ist=0;ist<=nst;ist++)
    for(int ism=0;ism<nsm;ism++)
      for(int jsm=0;jsm<nsm;jsm++)
	if(sign[ism*nsm+jsm]!=0)
	  if(ist==nst||(plot[ist][ism]&&plot[ist][jsm]))
	    {
	      //subtracted corr
	      jvec temp_corr=corr[ism*nsm+jsm];	    
	      for(int jst=0;jst<ist;jst++)
		for(int t=0;t<=TH;t++)
		  temp_corr[t]-=fun_migrad_A_fit(Z[jst*nsm+ism],Z[jst*nsm+jsm],M[jst],t,parity[ism*nsm+jsm]);
	      
	      //summ of exponentials
	      int npoints=tmax_fit[ism*nsm+jsm]-tmin_fit[ism*nsm+jsm];
	      if(npoints!=0)
		{
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
			  y0+=fun_migrad_A_fit(Z[jst*nsm+ism],Z[jst*nsm+jsm],M[jst],x0,parity[ism*nsm+jsm]);
			  y1+=fun_migrad_A_fit(Z[jst*nsm+ism],Z[jst*nsm+jsm],M[jst],x1,parity[ism*nsm+jsm]);
			}
		      temp_fun[i]=y0;//effective_masx(x,;//-log(y1/y0)/dx;
		    }
		  temp_fun=effective_mass(x,temp_fun,TH,parity[ism*nsm+jsm]);
		  
		  ofstream out(combine("%s/fit_B_st_%d_sm_%s_%s.xmg",plotpath,ist,sm[ism],sm[jsm]).c_str());
		  //error of the line
		  if(ist!=nst)
		    {
		      //summ of expo
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
		      
		      if(ist!=nst-1)
			{
			  //single mass
			  out<<"@s2 line type 1"<<endl;      
			  out<<"@s2 line color 8"<<endl;
			  out<<"@s2 fill color 8"<<endl;
			  out<<"@s2 fill type 1"<<endl;
			  out<<"@type xy"<<endl;
			  for(int i=0;i<npoints-1;i++) out<<x[i]<<" "<<M[ist].med()+M[ist].err()<<endl;
			  for(int i=npoints-2;i>=0;i--) out<<x[i]<<" "<<M[ist].med()-M[ist].err()<<endl;
			  out<<"&"<<endl;
			  //central line
			  out<<"@s3 line color 3"<<endl;
			  out<<"@type xy"<<endl;      
			  for(int i=0;i<npoints-1;i++) out<<x[i]<<" "<<M[ist].med()<<endl;
			}
		    }
		  //plot the original data with error  
		  out<<"&"<<endl;
		  out<<"@type xydy"<<endl;      
		  out<<"@s4 line type 0"<<endl;      
		  out<<"@s4 symbol color 1"<<endl;
		  out<<"@s4 errorbar color 1"<<endl;
		  out<<"@s4 symbol 1"<<endl;
		  out<<effective_mass(temp_corr,TH,parity[ism*nsm+jsm]);
		  out<<"&"<<endl;
		  out.close();
		}
	    }
  
  //plot together all the different sink smearing for a fixed source level
  for(int ist=0;ist<=nst;ist++)
    for(int ism=0;ism<nsm;ism++)
      {
	ofstream out(combine("%s/fit_B_st_%d_sm_%s_XX.xmg",plotpath,ist,sm[ism]).c_str());
	
	for(int jsm=0;jsm<nsm;jsm++)
	  if(sign[ism*nsm+jsm]!=0)
	    if(ist==nst||(plot[ist][ism]&&plot[ist][jsm]))
	      {
		//subtracted corr
		jvec temp_corr=corr[ism*nsm+jsm];	    
		for(int jst=0;jst<ist;jst++)
		  for(int t=0;t<=TH;t++)
		    temp_corr[t]-=fun_migrad_A_fit(Z[jst*nsm+ism],Z[jst*nsm+jsm],M[jst],t,parity[ist*nsm+jsm]);
		
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
			y0+=fun_migrad_A_fit(Z[jst*nsm+ism],Z[jst*nsm+jsm],M[jst],x0,parity[ist*nsm+jsm]);
			y1+=fun_migrad_A_fit(Z[jst*nsm+ism],Z[jst*nsm+jsm],M[jst],x1,parity[ist*nsm+jsm]);
		      }
		    temp_fun[i]=y0;//effective_masx(x,;//-log(y1/y0)/dx;
		  }
		temp_fun=effective_mass(x,temp_fun,TH,parity[ism*nsm+jsm]);
		
		//error of the line
		out<<"@s"<<jsm*3+0<<" line type 1"<<endl;      
		out<<"@s"<<jsm*3+0<<" line color 7"<<endl;
		out<<"@s"<<jsm*3+0<<" fill color 7"<<endl;
		out<<"@s"<<jsm*3+0<<" fill type 1"<<endl;
		out<<"@type xy"<<endl;
		for(int i=0;i<npoints-1;i++) out<<x[i]<<" "<<temp_fun[i].med()+temp_fun[i].err()<<endl;
		for(int i=npoints-2;i>=0;i--) out<<x[i]<<" "<<temp_fun[i].med()-temp_fun[i].err()<<endl;
		out<<"&"<<endl;
		//central line
		out<<"@s"<<jsm*3+1<<" line color 1"<<endl;
		out<<"@type xy"<<endl;      
		for(int i=0;i<npoints-1;i++) out<<x[i]<<" "<<temp_fun[i].med()<<endl;
		//plot the original data with error  
		out<<"&"<<endl;
		out<<"@type xydy"<<endl;      
		out<<"@s"<<jsm*3+2<<" line type 0"<<endl;      
		out<<"@s"<<jsm*3+2<<" symbol color 1"<<endl;
		out<<"@s"<<jsm*3+2<<" errorbar color 1"<<endl;
		out<<"@s"<<jsm*3+2<<" symbol 1"<<endl;
		out<<effective_mass(temp_corr,TH,parity[ism*nsm+jsm]);
		out<<"&"<<endl;
	      }
	out.close();
      }
}

int main(int narg,char **arg)
{
  if(narg<2) crash("use %s input",arg[0]);
  read_input(arg[1]);
  
  jvec corr[nst][nsm*nsm];
  jvec M(nst,njacks),Z(nst*nsm,njacks);

  ////////////////////// Stage A, loop over different states fitting them one by one //////////////////////

  cout<<"Stage A"<<endl;
  for(int ist=0;ist<nst;ist++)
    {
      for(int ism=0;ism<nsm;ism++)
	for(int jsm=0;jsm<nsm;jsm++)
	  if(sign[ism*nsm+jsm]!=0) //exclude if sign is 0
	    //if it is the first state, load, otherwise subtract previously fitted state
	    if(ist==0)
	      {
		//load and put the correct sign
		jvec tem=jvec_load(infile,T,njacks,ism*nsm+jsm);
		int loc_sign=(tem[1].med()>0)*2-1;
		tem=loc_sign*tem;
		
		//take correct parity
		parity[ism*nsm+jsm]=2*((tem[1]*tem[T-2]).med()>0)-1;
		corr[ist][ism*nsm+jsm]=tem.simmetrized(parity[ism*nsm+jsm]);
		cout<<"(parity,sign)["<<ism<<"]["<<jsm<<"]: "<<parity[ism*nsm+jsm]<<" "<<loc_sign<<endl;
		
		corr[ist][ism*nsm+jsm].print_to_file(combine("%s/raw_data_%s_%s.xmg",plotpath,sm[ism],sm[jsm]).c_str());
	      }
	    else
	      if(use[ist-1][ism]&&use[ist-1][jsm])
		{
		  corr[ist][ism*nsm+jsm]=corr[ist-1][ism*nsm+jsm];
		  for(int t=0;t<=TH;t++)
		    corr[ist][ism*nsm+jsm][t]-=
		      fun_migrad_A_fit(Z[(ist-1)*nsm+ism],Z[(ist-1)*nsm+jsm],M[ist-1],t,parity[ist*nsm+jsm]);
		}
      
      //perform the fit
      stage_A_fit(M,Z,corr[ist],ist);
    }
  
  if(nst>=2)
    {
      cout<<"M'/M: "<<smart_print(M[1]/M[0])<<endl;
      cout<<"f'/f: "<<smart_print(f(Z[1*nsm+0],M[1])/f(Z[0*nsm+0],M[0]))<<endl;
    }
  
  if(nst>=3)
    {
      cout<<"M''/M: "<<smart_print(M[2]/M[0])<<endl;
      cout<<"f''/f: "<<smart_print(f(Z[2*nsm+0],M[2])/f(Z[0*nsm+0],M[0]))<<endl;
    }
  
  ////////////////// Stage B, fit globally all the states on the relevant interval /////////////////
  
  cout<<endl<<"Stage B"<<endl;
  
  stage_B_fit(M,Z,corr[0]);
  
  if(nst>=2)
    {
      cout<<"M'/M: "<<smart_print(M[1]/M[0])<<endl;
      cout<<"f'/f: "<<smart_print(f(Z[1*nsm+0],M[1])/f(Z[0*nsm+0],M[0]))<<endl;
    }
  cout<<"(f: "<<f(Z[0*nsm+0],M[0])<<")"<<endl;
  cout<<"(f: "<<smart_print(f(Z[0*nsm+0],M[0])/lat_med[ibeta])<<" GeV )"<<endl;
  
  if(nst>=3)
    {
      cout<<"M''/M: "<<smart_print(M[2]/M[0])<<endl;
      cout<<"f''/f: "<<smart_print(f(Z[2*nsm+0],M[2])/f(Z[0*nsm+0],M[0]))<<endl;
    }
  
  //write the effect of the smearing
  cout<<"smearing effect on states:"<<endl;
  for(int ist=0;ist<nst;ist++)
    {
      cout<<" 1";
      for(int ism=1;ism<nsm;ism++) cout<<"\t"<<smart_print(Z[ist*nsm+ism]/Z[ist*nsm+0]);
      cout<<endl;
    }
  cout<<endl;
  
  //write the composition of the interpolating operators
  cout<<"relative composition of interpolating operators:"<<endl;
  for(int ism=0;ism<nsm;ism++)
    {
      jack tot=sqr(Z[0*nsm+ism]);
      for(int ist=1;ist<nst;ist++) tot+=sqr(Z[ist*nsm+ism]);
      
      cout<<" op "<<ism<<": 0";
      for(int ist=0;ist<nst;ist++)
	{
	  jack compo_sq=sqr(Z[ist*nsm+ism])/tot;
	  jack compo_sqrt=sqrt(compo_sq);
	  cout<<"\t"<<smart_print(Z[ist*nsm+ism]/Z[0*nsm+ism]);
	}
      cout<<endl;
    }
  cout<<endl;
  
  //print the relative error of the suppression of the ground state in the 3pts
  if(nst>=2)
    {
      int ist=0,jst=1;
      
      cout<<"Relative error between state "<<ist<<" and "<<jst<<": "<<endl;
      for(int ism=0;ism<nsm;ism++)
	for(int jsm=0;jsm<ism;jsm++)
	  {
	    jack C=(Z[ist*nsm+ism]*Z[jst*nsm+jsm]-Z[jst*nsm+ism]*Z[ist*nsm+jsm]);
	    jack D=Z[ist*nsm+ism]*Z[ist*nsm+jsm];
	    jack P1=C/D.err();
	    jack P2=C/D;
	    cout<<" using operators "<<ism<<" and "<<jsm<<": "<<smart_print(P1)<<", "<<smart_print(P2)<<endl;
	  }
    }
  
  //save M and Z
  M.write_to_binfile(outfile);
  Z.append_to_binfile(outfile);
  
  return 0;
}
