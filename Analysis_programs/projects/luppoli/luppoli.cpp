#include "load.cpp"

const int L=48,LH=L/2,VH=(LH+1)*(LH+1)*(LH+1);
int clust_size;
int njacks,nterm;

int dmin_termoli_fit,dmax_termoli_fit;
double *c_termoli_fit,*e_termoli_fit;

template <class T> T fun_termoli_migrad_fit(T A,T SIG,T C,double d)
{return -A/d+SIG*d+C;}

void ch2_termoli_migrad_fit(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  double A=p[0];
  double SIG=p[1];
  double C=p[2];

  for(int d=dmin_termoli_fit;d<=dmax_termoli_fit;d++)
    {
      double num=c_termoli_fit[d];
      double teo=fun_termoli_migrad_fit(A,SIG,C,d);
      double diff=num-teo;
      double err=e_termoli_fit[d];
      double cont=sqr(diff/err);
      if(!std::isnan(err)) ch+=cont;
      //if(flag==3)
      //cout<<" A: "<<A<<", SIG: "<<SIG<<", C: "<<C<<", d="<<d<<", diff=("<<num<<"-"<<teo<<")="<<diff<<" err="<<err<<" cont="<<cont<<endl;
    }
}

void termoli_migrad_fit(jack &A,jack &SIG,jack &C,jvec corr,int dmin,int dmax,const char *path=NULL)
{
  TMinuit minu;
  minu.SetPrintLevel(-1);
  minu.SetFCN(ch2_termoli_migrad_fit);

  minu.DefineParameter(0,"A",1,0.001,0,0);
  minu.DefineParameter(1,"SIG",1,0.001,0,0);
  minu.DefineParameter(2,"C",1,0.001,0,0);

  c_termoli_fit=new double[LH+1];
  e_termoli_fit=new double[LH+1];

  dmin_termoli_fit=dmin;
  dmax_termoli_fit=dmax;

  for(int iel=0;iel<LH+1;iel++)
    e_termoli_fit[iel]=corr[iel].err();

  for(int ijack_fit=0;ijack_fit<=njacks;ijack_fit++)
    {
      for(int iel=0;iel<LH+1;iel++) c_termoli_fit[iel]=corr[iel][ijack_fit];
      minu.Migrad();
      double dum;
      minu.GetParameter(0,A.data[ijack_fit],dum);
      minu.GetParameter(1,SIG.data[ijack_fit],dum);
      minu.GetParameter(2,C.data[ijack_fit],dum);
    }
  int ndof=dmax-dmin+1-3;
  
  double ch2,grad[3],par[3]={A[njacks],SIG[njacks],C[njacks]};
  minu.Eval(3,grad,ch2,par,3);
  if(debug_fit) cout<<"A: "<<smart_print(A)<<", SIG: "<<smart_print(SIG)<<", C: "<<smart_print(C)<<", ch2: "<<ch2<<"/"<<ndof<<"="<<ch2/ndof<<endl;

  if(path!=NULL)
    {
      int npoints=100;
      double x[npoints];

      jvec y(npoints,njacks);

      for(int ip=0;ip<npoints;ip++)
	{
	  x[ip]=dmin+(double)(dmax-dmin)/(npoints-1)*ip;
	  y[ip]=-A/x[ip]+SIG*x[ip]+C;
	}

      ofstream out(path);
      out<<write_polygon(x,y);
      out<<write_ave_line(x,y);
      out<<"@type xydy"<<endl;
      out<<"@s2 line type 0"<<endl;
      out<<corr<<endl;
    }
}

void read(int dmin,int dmax)
{
  //read data
  double *ll_dag_raw,*ll_raw,*l_raw;
  int nconfs;
  read(nconfs,ll_dag_raw,ll_raw,l_raw,L,"Luppoli");
  
  //find cluster size
  njacks=(nconfs-nterm)/clust_size;
  int nconfs_temp=clust_size*njacks+nterm;
  if(nconfs_temp!=nconfs) cout<<"Changing nconfs from "<<nconfs<<" to "<<nconfs_temp<<endl;
  nconfs=nconfs_temp;
  
  //clusterize
  jvec ll(2*VH,njacks),ll_dag(2*VH,njacks),l(2,njacks);
  ll=ll_dag=0;
  l=0;
  for(int ijack=0;ijack<njacks;ijack++)
    for(int iconf=nterm+ijack*clust_size;iconf<nterm+(ijack+1)*clust_size;iconf++)
      {
	for(int i=0;i<2*VH;i++)
	  {
	    ll[i][ijack]+=ll_raw[i+2*VH*iconf];
	    ll_dag[i][ijack]+=ll_dag_raw[i+2*VH*iconf];
	  }
	for(int i=0;i<2;i++) l[i][ijack]+=l_raw[i+2*iconf];
      }
  
  //free
  delete[] ll_raw;
  delete[] ll_dag_raw;
  delete[] l_raw;
  
  //clusterize
  for(int i=0;i<2*VH;i++)
    {
      ll[i].clusterize(clust_size);
      ll_dag[i].clusterize(clust_size);
    }
  for(int i=0;i<2;i++) l[i].clusterize(clust_size);
  
  cout<<"Corr_ll[0]: \t"<<smart_print(ll[0])<<","<<smart_print(ll[1])<<endl;
  cout<<"Corr_ll+[0]:\t"<<smart_print(ll_dag[0])<<","<<smart_print(ll_dag[1])<<endl;
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //build magnetic correlators
  
  //ouptut and filter
  jvec corr_orth_mag[3],corr_orth_ele[3];
  for(int mu=0;mu<3;mu++) corr_orth_mag[mu]=corr_orth_ele[mu]=jvec(LH+1,njacks);
  ofstream allpoints_plot("plots/all_points.xmg");
  allpoints_plot<<"@type xydy"<<endl;
  allpoints_plot<<"@s0 line type 0"<<endl;
  for(int i=0;i<2*VH;i+=2)
    {
      //compute distance d2
      int c[3];
      int j=i/2;
      int d2=0;
      for(int mu=2;mu>=0;mu--)
	{
	  c[mu]=j%(LH+1);
	  j/=LH+1;
	  
	  d2+=c[mu]*c[mu];
	}
      
      //subtract
      jack Mp=0.5*(ll_dag[i]+ll[i])-(l[0]*l[0]+l[1]*l[1]);
      //jack Mm=0.5*(out_dag[i+1]+out[i+1]);
      jack Em=0.5*(ll_dag[i]-ll[i]);
      //jack Ep=0.5*(out_dag[i+1]-out[i+1]);
      Mp*=sqrt(d2);
      //Mm*=sqrt(d2);
      Em*=sqrt(d2);
      //Ep*=sqrt(d2);
      
      for(int mu=0;mu<3;mu++)
	if(c[(mu+1)%3]==0 && c[(mu+2)%3]==0)
	  {
	    corr_orth_mag[mu][c[mu]]=Mp;
	    corr_orth_ele[mu][c[mu]]=Em;
	  }
      if(!std::isnan(Em.err())) allpoints_plot<<sqrt(d2)<<" "<<Em<<endl;
    }
  
  ofstream orth_dirs_plot("plots/orth_dirs.xmg");
  orth_dirs_plot<<"@type xydy"<<endl;
  for(int mu=0;mu<3;mu++)
    orth_dirs_plot<<corr_orth_mag[mu]<<"&"<<endl;
  orth_dirs_plot<<"&"<<endl;
  for(int mu=0;mu<3;mu++)
    orth_dirs_plot<<corr_orth_ele[mu]<<"&"<<endl;
  orth_dirs_plot<<"&"<<endl;
  
  //jvec A(3,njacks),SIG(3,njacks),C(3,njacks);
  //for(int mu=0;mu<3;mu++) termoli_migrad_fit(A[mu],SIG[mu],C[mu],-log(corr_orth[mu]),dmin,dmax,combine("plots/termoli_fit_%d.xmg",mu).c_str());
  
  /*
  cout<<"----------------------------"<<endl;
  cout<<"A[1]: "<<smart_print(A[1]-A[0])<<endl;
  cout<<"SIGMA[1]: "<<smart_print(SIG[1]-SIG[0])<<endl;
  cout<<"C[1]: "<<smart_print(C[1]-C[0])<<endl;
  cout<<"----------------------------"<<endl;
  cout<<"A[2]: "<<smart_print(A[2]-(A[1]+A[0])/2)<<endl;
  cout<<"SIGMA[2]: "<<smart_print(SIG[2]-(SIG[1]+SIG[0])/2)<<endl;
  cout<<"C[2]: "<<smart_print(C[2]-(C[1]+C[0])/2)<<endl;
  cout<<"----------------------------"<<endl;
  */
}

int main(int narg,char **arg)
{
  if(narg<5) crash("use: %s clust_size term dmin dmax",arg[0]);
  clust_size=atoi(arg[1]);
  nterm=atoi(arg[2]);
  int dmin=atoi(arg[3]);
  int dmax=atoi(arg[4]);
  cout<<"dmin: "<<dmin<<", dmax: "<<dmax<<endl;
  read(dmin,dmax);
  
  return 0;
}
