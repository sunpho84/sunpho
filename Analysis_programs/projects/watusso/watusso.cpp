#include "../../src/include.h"
#include <TMath.h>

using namespace TMath;

int nmu=3;
int nrho=2;
int njacks;
int nape;
int dmax;
vector<int> sizes;
const int nind_planes=1;

const int debug=0;

TMinuit minu(3);

void discard(ifstream &fin,int n)
{
  if(debug) cout<<"Discharging "<<n<<endl;
  char head[1000];
  for(int idisc=0;idisc<n;idisc++)
    {
      fin>>head;
      //if(!fin.good()) crash()
      if(debug) cout<<head<<" ";
    }
  if(debug) cout<<endl;
}

int ind_ind(int iape,int isize,int iind)
{return iape+nape*(isize+sizes.size()*iind);}

int ind(int iape,int isize,int imu)
{return imu+nmu*(isize+sizes.size()*iape);}

int ind(int iape,int isize,int imu,int irho)
{return irho+nrho*ind(iape,isize,imu);}

jack BesselK0(jack a)
{jack out(njacks);for(int ijack=0;ijack<=njacks;ijack++) out[ijack]=BesselK0(a[ijack]);return out;}
jack BesselK1(jack a)
{jack out(njacks);for(int ijack=0;ijack<=njacks;ijack++) out[ijack]=BesselK1(a[ijack]);return out;}

//eq. 3.80 saldatore
template <class T1,class T2> T1 fun_fit(T1 phi,T1 mu,T1 alpha,T2 x)
{return phi/(2*M_PI)*sqr(mu)/alpha*BesselK0(sqrt(sqr(mu*x)+sqr(alpha)))/BesselK1(alpha);}

//fit the mass and the matrix element in SS and SL combo
vector<double> data,err;
int dmin_fit=0;
int dmax_fit;
double plot_ymax,plot_ymin;

void ch2_fit(int &npar,double *fuf,double &ch,double *p,int flag)
{
  double phi=p[0];
  double mu=p[1];
  double alpha=p[2];
  
  ch=0;
  for(int x=dmin_fit;x<=dmax_fit;x++)
    {
      double num=data[x];
      double teo=fun_fit(phi,mu,alpha,x);
      double diff=num-teo;
      double cont=sqr(diff/err[x]);
      ch+=cont;
    }
}

jvec fit(const char *path,jvec &corr,int isize,int iape)
{
  //find max range
  int nexcl=0;//2
  dmax_fit=dmin_fit;
  while(corr[dmax_fit].med()>corr[dmax_fit].err() && dmax_fit<dmax) dmax_fit++;
  dmax_fit=max(dmax_fit-nexcl,dmin_fit);
  
  int ndof=dmax_fit-3-dmin_fit;
  
  //copy err
  for(int i=0;i<=dmax;i++) err[i]=corr[i].err();
  
  //fit
  jack phi(njacks),mu(njacks),alpha(njacks),ch2(njacks);
  for(int ijack=0;ijack<=njacks;ijack++)
    {
      //if(ijack==njacks) minu.SetPrintLevel(1);
      
      for(int i=0;i<=dmax;i++) data[i]=corr[i][ijack];
      minu.Migrad();
      double dum;
      minu.GetParameter(0,phi[ijack],dum);
      minu.GetParameter(1,mu[ijack],dum);
      minu.GetParameter(2,alpha[ijack],dum);
      
      double pars[3]={phi[ijack],mu[ijack],alpha[ijack]};
      minu.Eval(3,NULL,ch2[ijack],pars,0);
    }
  
  if(debug) cout<<"Phi: "<<smart_print(phi)<<endl;
  if(debug) cout<<"Mu: "<<smart_print(mu)<<endl;
  if(debug) cout<<"Alpha: "<<smart_print(alpha)<<endl;
  
  int npoints=200;
  jvec fun_teo(npoints,njacks);
  double x[npoints];
  for(int i=0;i<npoints;i++)
    {
      double dx=(double)(dmax_fit-dmin_fit+0.5)/(npoints-1);
      x[i]=dmin_fit+dx*i;
      fun_teo[i]=fun_fit(phi,mu,alpha,x[i]);
    }
  
  //preamble
  ofstream out(path);
  out<<
    "@ world 0, "<<plot_ymin*0.5<<", "<<dmax+0.5<<", "<<plot_ymax*2<<endl<<
    "@ yaxes scale Logarithmic"<<endl<<
    "@ yaxis  tick major 10"<<endl<<
    "@ yaxis  tick minor ticks 9"<<endl<<
    "@ s0 symbol color 10"<<endl<<
    "@ s0 symbol fill color 10"<<endl<<
    "@ s0 line color 10"<<endl<<
    "@ s0 fill type 1"<<endl<<
    "@ s0 fill color 10"<<endl<<
    "@ s0 errorbar color 10"<<endl<<
    "@ s1 symbol color 12"<<endl<<
    "@ s1 symbol fill color 12"<<endl<<
    "@ s1 line linewidth 2.0"<<endl<<
    "@ s1 line color 12"<<endl<<
    "@ s1 errorbar color 12"<<endl<<
    "@ s2 symbol 1"<<endl<<
    "@ s2 symbol color 12"<<endl<<
    "@ s2 symbol fill color 12"<<endl<<
    "@ s2 symbol fill pattern 1"<<endl<<
    "@ s2 symbol linewidth 2.0"<<endl<<
    "@ s2 line type 0"<<endl<<
    "@ s2 line color 12"<<endl<<
    "@ s2 errorbar color 12"<<endl<<
    "@ s2 errorbar linewidth 2.0"<<endl<<
    "@ s2 errorbar riser linewidth 2.0"<<endl<<
    "@ s3 symbol 1"<<endl<<
    "@ s3 symbol color 12"<<endl<<
    "@ s3 symbol fill color 12"<<endl<<
    "@ s3 symbol fill pattern 2"<<endl<<
    "@ s3 line type 0"<<endl<<
    "@ s3 line color 12"<<endl<<
    "@ s3 errorbar color 12"<<endl<<
    "@ s3 errorbar linewidth 2.0"<<endl<<
    "@ s3 errorbar linestyle 2"<<endl<<
    "@ s3 errorbar riser linewidth 2.0"<<endl<<
    "@ s3 errorbar riser linestyle 2"<<endl<<
    "@with string"<<endl<<
    "@    string on"<<endl<<
    "@    string loctype view"<<endl<<
    "@    string 0.267678100264, 0.417308707124"<<endl<<
    "@    string color 1"<<endl<<
    "@    string char size 1.700000"<<endl<<
    "@    string def \""
    "\\xf\\0="<<smart_print(phi)<<"\\n"<<
    "\\xm\\0="<<smart_print(mu)<<"\\n"<<
    "\\xa\\0="<<smart_print(alpha)<<"\\n\\n"<<
    "\\xc\\0\\S2\\N/d.o.f="<<smart_print(ch2)<<"/"<<ndof<<"\\n\""<<endl<<
    "@with string"<<endl<<
    "@    string on"<<endl<<
    "@    string loctype view"<<endl<<
    "@    string 0.846965699208, 0.72364116095"<<endl<<
    "@    string color 1"<<endl<<
    "@    string char size 1.500000"<<endl<<
    "@    string def \""
    "Size="<<sizes[isize]<<"\\n"<<
    "N\\sAPE\\N="<<iape*10<<"\\n\""<<endl;
  
  for(int i=0;i<npoints;i++)          out<<x[i]<<" "<<fun_teo[i].med()-fun_teo[i].err()<<endl;
  for(int i=npoints-1;i>=0;i--)       out<<x[i]<<" "<<fun_teo[i].med()+fun_teo[i].err()<<endl;
  out<<"&"<<endl;
  for(int i=0;i<npoints;i++)          out<<x[i]<<" "<<fun_teo[i].med()<<endl;
  out<<"&"<<endl<<"@type xydy"<<endl;
  for(int i=0;i<=dmax_fit;i++)        out<<i<<" "<<corr[i]<<endl;
  out<<"&"<<endl;
  for(int i=dmax_fit+1;i<=dmax;i++)   out<<i<<" "<<corr[i]<<endl;
  out<<"&"<<endl;
  
  jvec pars(3,njacks);
  pars[0]=phi;
  pars[1]=mu;
  pars[2]=alpha;
  
  return pars;
}

int main()
{
  const char path[]="watusso";
  ifstream fin(path);
  if(!fin.good()) crash("opening the file \"%s\"",path);
  
  FILE *infile=open_file("input","r");
  int sizemin,sizemax,sizestep;
  int nconfs;
  read_formatted_from_file_expecting((char*)&nconfs,infile,"%d","NConfs");
  cout<<"Number configurations: "<<nconfs<<endl;
  njacks= sqrt(nconfs);
  //njacks=nconfs;
  int clust_size=nconfs/njacks;
  cout<<"Cluster size: "<<clust_size<<endl;
  cout<<"Njacks: "<<njacks<<endl;
  nconfs=njacks*clust_size;
  cout<<"Adjusted NConfs: "<<nconfs<<endl;
  read_formatted_from_file_expecting((char*)&nape,infile,"%d","NApe");
  read_formatted_from_file_expecting((char*)&dmax,infile,"%d","DMax");
  read_formatted_from_file_expecting((char*)&sizemin,infile,"%d","SizeMin");
  read_formatted_from_file_expecting((char*)&sizestep,infile,"%d","SizeStep");
  read_formatted_from_file_expecting((char*)&sizemax,infile,"%d","SizeMax");
  fclose(infile);
  
  for(int size=sizemin;size<=sizemax;size+=sizestep) sizes.push_back(size);
  
  // vector<ofstream> cali_out(nape);
  // for(int iape=0;iape<nape;iape++)
  //   {
  //     cali_out[iape].open(combine("/tmp/wils_sm%02d.dat",iape*10));
  //     cali_out[iape].precision(14);
  // }
  
  const char perp_dir[]="xyz";
  const char perp2_dir[3][3]={"yz","xz","xy"};
  
  //define the average
  vector<jvec> ave(ind_ind(nape-1,sizes.size()-1,nind_planes-1)+1,jvec(dmax+1,njacks));
  for(size_t i=0;i<ave.size();i++) ave[i]=0;
  
  int deg[2];
  if(nind_planes==2)
    {
      deg[0]=2;
      deg[1]=4;
    }
  else deg[0]=6;
  
  if(!file_exists("watusso_bin"))
    {
      jvec w(nmu*sizes.size()*nape,njacks);
      jvec t(nmu*sizes.size()*nape,njacks);
      w=t=0;
      
      //init output
      int ncorrs=ind(nape-1,sizes.size()-1,nmu-1,nrho-1)+1;
      vector<jvec> corr(ncorrs,jvec(2*dmax+1,njacks));
      for(int icorr=0;icorr<ncorrs;icorr++) corr[icorr]=0;
      
      for(int iconf=0;iconf<nconfs;iconf++)
	{
	  int ijack=iconf/clust_size;
	  
	  discard(fin,3);
	  int jconf;
	  fin>>jconf;
	  for(int iape=0;iape<nape;iape++)
	    {
	      for(int imu=0;imu<nmu;imu++)
		{
		  discard(fin,6);
		  int nape;
		  fin>>nape;
		  discard(fin,8);
		  int mu;
		  fin>>mu;
		  discard(fin,3);
		  double plaqr,plaqi;
		  fin>>plaqr>>plaqi;
		  
		  for(size_t isize=0;isize<sizes.size();isize++)
		    {
		      discard(fin,3);
		      int size;
		      fin>>size;
		      discard(fin,3);
		      
		      //load w
		      double wr,wi;
		      fin>>wr>>wi;
		      w[ind(iape,isize,imu)][ijack]+=wr;
		      t[ind(iape,isize,imu)][ijack]+=wi;
		      
		      for(int irho=0;irho<nrho;irho++)
			{
			  discard(fin,3);
			  double rho;
			  fin>>rho;
			  
			  for(int d=0;d<2*dmax+1;d++)
			    {
			      int e;
			      double wpr,wpi,pr,pi;
			      fin>>e>>wpr>>wpi>>pr>>pi;
			      if(!fin.good()) crash("reading conf %d: %d %lg %lg",iconf,e,wpr,wpi);
			      corr[ind(iape,isize,imu,irho)][d][ijack]+=wpr-pr/3;
			      
			      // cali_out[iape]<<iconf<<" "<<imu<<" "<<irho<<" "<<d-dmax<<" "<<size<<" "<<size<<" "<<wpr<<" "<<wpi<<" "<<pr<<" "<<pi<<" "<<wr<<" "<<wi<<" "<<plaqr<<" "<<plaqi<<endl;
			    }
			}
		    }
		}
	    }
	}
      cout<<"Finished reading txt"<<endl;
      
      //clusterize w and t
      w.clusterize(clust_size);
      t.clusterize(clust_size);
      
      for(int irho=0;irho<nrho;irho++)
	for(int imu=0;imu<nmu;imu++)
	  for(size_t isize=0;isize<sizes.size();isize++)
	    {
	      //clusterize
	      for(int iape=0;iape<nape;iape++) corr[ind(iape,isize,imu,irho)].clusterize(clust_size);
	      
	      //normalize
	      for(int iape=0;iape<nape;iape++) corr[ind(iape,isize,imu,irho)]/=w[ind(iape,isize,imu)];
	      
	      //find independent plane
	      int iind;
	      if(nind_planes==2) iind=((perp_dir[imu]=='z')||(perp2_dir[imu][irho]=='z'));
	      else iind=0;
	      
	      //symmetrize
	      for(int iape=0;iape<nape;iape++)
		for(int d=0;d<=dmax;d++)
		  ave[ind_ind(iape,isize,iind)][d]+=(corr[ind(iape,isize,imu,irho)][dmax+d]+corr[ind(iape,isize,imu,irho)][dmax-d])/(2*deg[iind]);
	    }
      
      //write
      ofstream wat_bin("watusso_bin");
      for(auto &it : ave) for(int i=0;i<it.nel;i++) wat_bin.write((char*)it[i].data,sizeof(double)*(njacks+1));
    }
  else
    {
      ifstream wat_bin("watusso_bin");
      for(auto &it : ave) for(int i=0;i<it.nel;i++) wat_bin.read((char*)it[i].data,sizeof(double)*(njacks+1));
      cout<<"Finished reading bin"<<endl;
    }
  
  minu.SetFCN(ch2_fit);
  minu.SetPrintLevel(-1);
  
  minu.DefineParameter(0,"phi",1,0.001,0,40);
  minu.DefineParameter(1,"mu",1,0.001,0,40);
  minu.DefineParameter(2,"alpha",1,0.001,0,40);
  
  data.resize(dmax+1);
  err.resize(dmax+1);
  
  for(int iind=0;iind<nind_planes;iind++)
    for(size_t isize=0;isize<sizes.size();isize++)
      {
	//fix min and max
	plot_ymax=-1e300;
	plot_ymin=+1e300;
	for(int iape=0;iape<nape;iape++)
	  for(int i=dmin_fit;i<=dmax;i++)
	    {
	      double m=ave[ind_ind(iape,isize,iind)][i].med();
	      double e=ave[ind_ind(iape,isize,iind)][i].err();
	      if(m>e)
		{
		  plot_ymax=max(plot_ymax,m+e);
		  plot_ymin=min(plot_ymin,m-e);
		}
	    }
	
	//fit
	vector<jvec> pars(nape,jvec(3,njacks));
	for(int iape=0;iape<nape;iape++)
	  pars[iape]=fit(combine("plots/fitted_size%02d_ape%02d_plane%d.xmg",sizes[isize],iape,iind).c_str(),ave[ind_ind(iape,isize,iind)],isize,iape);
	
	//find max and min for each par
	double min_pars[3],max_pars[3];
	for(int i=0;i<3;i++)
	  {
	    min_pars[i]=+1e300;
	    max_pars[i]=-1e300;
	    for(int iape=0;iape<nape;iape++)
	      {
		min_pars[i]=min(min_pars[i],pars[iape][i].med()-2*pars[iape][i].err());
		max_pars[i]=max(max_pars[i],pars[iape][i].med()+2*pars[iape][i].err());
	      }
	  }
	
	ofstream out_pars[3];
	const char par_name[3][10]={"phi","mu","alpha"},par_initial[4]="fma";
	for(int i=0;i<3;i++)
	  {
	    out_pars[i].open(combine("pars/%s_size%02d_plane%d.xmg",par_name[i],sizes[isize],iind).c_str());
	    if(!out_pars[i].good()) crash("check pars folder");
	    out_pars[i]<<
	      "@with string"<<endl<<
	      "@  string on"<<endl<<
	      "@  string loctype view"<<endl<<
	      "@  string 0.660817941953, 0.51715039578"<<endl<<
	      "@  string color 1"<<endl<<
	      "@  string rot 0"<<endl<<
	      "@  string font 0"<<endl<<
	      "@  string just 0"<<endl<<
	      "@  string char size 4.617346"<<endl<<
	      "@  string def \"\\x"<<par_initial[i]<<"\\N\""<<endl<<
	      "@world -10, "<<min_pars[i]<<", "<<nape*10<<", "<<max_pars[i]<<""<<endl<<
	      "@s0 symbol 1"<<endl<<
	      "@s0 symbol color 12"<<endl<<
	      "@s0 symbol fill color 12"<<endl<<
	      "@s0 symbol fill pattern 1"<<endl<<
	      "@s0 line type 0"<<endl<<
	      "@s0 line color 12"<<endl<<
	      "@s0 errorbar color 12"<<endl<<
	      "@s0 errorbar linewidth 2.0"<<endl<<
	      "@s0 errorbar riser linewidth 2.0"<<endl<<
	      "@ type xydy"<<endl<<
	      endl;
	  }
	
	//print
	for(int iape=0;iape<nape;iape++)
	  for(int i=0;i<3;i++) out_pars[i]<<iape*10<<" "<<pars[iape][i]<<endl;
      }
  
  return 0;
}
