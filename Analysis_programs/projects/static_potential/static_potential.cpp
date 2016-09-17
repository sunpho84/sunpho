#include "include.h"

int njacks;
int term_conf;
char data_path[200];
int Tmin,Tmax,Dmin,Dmax,n_sm,conf_min,conf_max;
jvec *data;
int dmin_fit,dmax_fit;
int *tmin_fit,*tmax_fit;

void read_input(const char *input_path)
{
  FILE *fin=open_file(input_path,"r");
  
  read_formatted_from_file_expecting((char*)&njacks,fin,"%d","njacks");
  read_formatted_from_file_expecting(data_path,fin,"%s","data_path");
  read_formatted_from_file_expecting((char*)&term_conf,fin,"%d","term_conf");
  read_formatted_from_file_expecting((char*)&dmin_fit,fin,"%d","dfit_int");
  read_formatted_from_file((char*)&dmax_fit,fin,"%d","dfit_int");
  tmin_fit=(int*)malloc(sizeof(int)*dmax_fit);
  tmax_fit=(int*)malloc(sizeof(int)*dmax_fit);
  for(int d=dmin_fit;d<dmax_fit;d++)
    {
      read_formatted_from_file_expecting((char*)&(tmin_fit[d]),fin,"%d",combine("d%d_tfit_int",d).c_str());
      read_formatted_from_file((char*)&(tmax_fit[d]),fin,"%d",combine("d%d_tfit_int",d).c_str());
    }
  
  fclose(fin);
}

void load_data()
{
  FILE *fin=open_file(data_path,"r");
  
  int nr=0;

  char *line=NULL;
  int rc;
  size_t nc=0;
  do
    {
      rc=getline(&line,&nc,fin);
      if(rc>0)
	{
	  double fuf;
	  int iconf,ism,t,d;
	  sscanf(line,"%d %d %d %d %lg",&iconf,&ism,&t,&d,&fuf);
	  if(nr==0||t>Tmax) Tmax=t;
	  if(nr==0||d>Dmax) Dmax=d;
	  if(nr==0||t<Tmin) Tmin=t;
	  if(nr==0||d<Dmin) Dmin=d;
	  if(nr==0||ism>n_sm) n_sm=ism;
	  if(nr==0||iconf<conf_min) conf_min=iconf;
	  if(nr==0||iconf>conf_max) conf_max=iconf;
	  nr++;
	}
    }
  while(rc>0);

  n_sm++;
  Dmax++;
  Tmax++;
  
  printf("conf available: %d %d\n",conf_min,conf_max);
  printf("T: %d %d\n",Tmin,Tmax);
  printf("D: %d %d\n",Dmin,Dmax);
  printf("n_sm: %d\n",n_sm);
  
  //get back
  fseek(fin,0,SEEK_SET);
  
  conf_min=max(term_conf,conf_min);
  conf_max=max(term_conf,conf_max);
  int nconf_ava=conf_max-conf_min+1;
  int clust_size=nconf_ava/njacks;
  int nconf=clust_size*njacks;
  conf_max=conf_min+nconf;
  
  printf("conf used: %d %d\n",conf_min,conf_max);
  
  //allocate data
  data=(jvec*)malloc(Dmax*sizeof(jvec));
  for(int d=Dmin;d<Dmax;d++)
    {
      data[d]=jvec(Tmax,njacks);
      data[d]*=0;
    }
  
  //clusterize
  do
    {
      rc=getline(&line,&nc,fin);
      if(rc>0)
	{
	  double fuf;
	  int iconf,ism,t,d;
	  sscanf(line,"%d %d %d %d %lg",&iconf,&ism,&t,&d,&fuf);
	  
	  if(iconf>=conf_min && iconf<conf_max)
	    {
	      int iclust=(iconf-conf_min)/clust_size;
	      data[d].data[t].data[iclust]+=fuf;
	    }
	}
    }
  while(rc>0);
  
  fclose(fin);
  
  //jack
  for(int d=Dmin;d<Dmax;d++)
    for(int t=Tmin;t<Tmax;t++)
      {
	
	for(int iclust=0;iclust<njacks;iclust++)
	  data[d].data[t].data[njacks]+=data[d].data[t].data[iclust];
	
	for(int iclust=0;iclust<njacks;iclust++)
	  data[d].data[t].data[iclust]=(data[d].data[t].data[njacks]-data[d].data[t].data[iclust])/(nconf-clust_size);
	
	data[d].data[t].data[njacks]/=nconf;
      } 
}

void extract_potential_from_correlators(jvec &V)
{
  ofstream out("/tmp/potential.xmg");
  for(int d=dmin_fit;d<dmax_fit;d++)
    {
      jvec corr=aperiodic_effective_mass(data[d]);
      V[d]=constant_fit(corr,tmin_fit[d],tmax_fit[d]);
      
      cout<<d<<" "<<V[d]<<endl;
      
      out<<"@type xydy"<<endl;
      out<<corr;
      out<<"&"<<endl;
      out<<write_constant_with_error(V[d],tmin_fit[d],tmax_fit[d]);
      out<<"&"<<endl;
    }
}

double *c_r0_fr_a_fit,*e_r0_fr_a_fit;
int dmin_r0_fr_a_fit;
int dmax_r0_fr_a_fit;

template <class T> T fun_r0_fr_a_migrad_fit(T A,T B,T sigma,double d)
{return A+B/(d+1)+sigma*(d+1);}

void ch2_r0_fr_a_migrad_fit(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  double A=p[0];
  double B=p[1];
  double sigma=p[2];
  
  for(int d=dmin_r0_fr_a_fit;d<dmax_r0_fr_a_fit;d++)
    {
      double num=c_r0_fr_a_fit[d];
      double teo=fun_r0_fr_a_migrad_fit(A,B,sigma,d);
      double diff=num-teo;
      double err=e_r0_fr_a_fit[d];
      double cont=sqr(diff/err);
      ch+=cont;
      if(flag==3) cout<<" d="<<d<<", diff=("<<num<<"-"<<teo<<")="<<diff<<" err="<<err<<" cont="<<cont<<endl;
    }
}

void r0_fr_a_migrad_fit(jack &A,jack &B,jack &sigma,jvec corr,int dmin,int dmax,const char *path=NULL)
{
  jvec ecorr=effective_mass(corr);
  
  A=B=sigma=jack(njacks);
  A=B=sigma=0;
  
  TMinuit minu;
  minu.SetPrintLevel(-1);
  minu.SetFCN(ch2_r0_fr_a_migrad_fit);
  minu.DefineParameter(0,"A",0.1,0.001,0,0);
  minu.DefineParameter(1,"B",0.1,0.001,0,0);
  minu.DefineParameter(2,"sigma",0.1,0.001,0,0);
  
  c_r0_fr_a_fit=new double[dmax];
  e_r0_fr_a_fit=new double[dmax];
  
  dmin_r0_fr_a_fit=dmin;
  dmax_r0_fr_a_fit=dmax;
  
  for(int iel=0;iel<dmax;iel++)
    e_r0_fr_a_fit[iel]=corr[iel].err();
  
  for(int ijack_fit=0;ijack_fit<=njacks;ijack_fit++)
    {
      for(int iel=0;iel<dmax;iel++) c_r0_fr_a_fit[iel]=corr[iel][ijack_fit];
      
      minu.Migrad();
      double dum;
      minu.GetParameter(0,A.data[ijack_fit],dum);
      minu.GetParameter(1,B.data[ijack_fit],dum);
      minu.GetParameter(2,sigma.data[ijack_fit],dum);
    }

  double ch2,grad[3],par[3]={A[njacks],B[njacks],sigma[njacks]};
  minu.Eval(3,grad,ch2,par,3);
  cout<<"A: "<<A<<", B: "<<B<<", sigma: "<<sigma<<", ch2: "<<ch2<<" / "<<dmax-dmin-3<<" = "<<ch2/(dmax-dmin-3)<<endl;
  
  if(path!=NULL)
    {
      ofstream out(path);

      //fit_function
      int npoints=100;
      double dx=((double)dmax-1-dmin)/(npoints-1);
      double x[npoints];
      jvec temp_fun(npoints,njacks);
      for(int i=0;i<npoints;i++)
	{
	  x[i]=dmin+i*dx;
	  temp_fun[i]=fun_r0_fr_a_migrad_fit(A,B,sigma,x[i]);
	}
      
      //errr band
      out<<"@s0 line type 1"<<endl;      
      out<<"@s0 line color 7"<<endl;
      out<<"@s0 fill color 7"<<endl;
      out<<"@s0 fill type 1"<<endl;
      out<<"@type xy"<<endl;
      for(int i=0;i<npoints-1;i++) out<<x[i]+1<<" "<<temp_fun[i].med()+temp_fun[i].err()<<endl;
      for(int i=npoints-2;i>=0;i--) out<<x[i]+1<<" "<<temp_fun[i].med()-temp_fun[i].err()<<endl;
      out<<"&"<<endl;
      //central line
      out<<"@s1 line color 1"<<endl;
      out<<"@type xy"<<endl;      
      for(int i=0;i<npoints-1;i++) out<<x[i]+1<<" "<<temp_fun[i].med()<<endl;
      //original data
      out<<"&"<<endl;
      out<<"@type xydy"<<endl;      
      out<<"@s4 line type 0"<<endl;      
      out<<"@s4 symbol color 1"<<endl;
      out<<"@s4 errorbar color 1"<<endl;
      out<<"@s4 symbol 1"<<endl;
      for(int d=dmin;d<dmax;d++)
	out<<d+1<<" "<<corr[d]<<endl;
      out.close();
    }  
}

jack determine_r0_fr_a(jvec V)
{
  jack out;
  
  jack A,B,sigma;
  r0_fr_a_migrad_fit(A,B,sigma,V,dmin_fit,dmax_fit,"potential_d_fit.xmg");
  
  return sqrt((1.65+B)/sigma);
}

int main(int narg,char **arg)
{
  if(narg<2) crash("use %s input",arg[0]);
  read_input(arg[1]);

  load_data();
  
  //determine the potential from correlators
  jvec V(Dmax,njacks);
  extract_potential_from_correlators(V);
  
  //fit the potential to determine r0/a
  jack r0_fr_a=determine_r0_fr_a(V);
  
  cout<<"r0/a: "<<r0_fr_a<<endl;
  cout<<"a: "<<0.44/r0_fr_a<<" fm"<<endl;
  
  return 0;
}
