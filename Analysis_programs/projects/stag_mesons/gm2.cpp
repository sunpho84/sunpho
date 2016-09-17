#include "include.h"

int nflavs;
int njacks;
int nconfs;
int nop;
int T;
double a;

//////////////////////////////////// to find phys point /////////////////////////////

int glb_tmin,glb_tmax;

void ignore_until_eol(ifstream &in)
{in.ignore(10000,'\n');}

int icorr(int iflav,int iop)
{return iop+nop*iflav;}

jvec calc_HVP(jvec corr)
{
  jvec HVP(T,njacks);
  
  //point at q=0
  HVP[0]=0;
  
  //all other q
  for(int iq=1;iq<T;iq++)
    {
      double q0=2*M_PI*iq/T;
      HVP[iq]=0;
      for(int t=0;t<=T/2;t++) HVP[iq]+=2*((cos(q0*t)-1)/sqr(q0)+sqr(t)/2)*corr[t];
    }
  return HVP/sqr(137.0*M_PI);
}

void write_HVP(const char *path,jvec corr,double val_0=0)
{
  corr[0]=val_0;
  ofstream out(path);
  out<<"@type xydy"<<endl;
  for(int iq=0;iq<T;iq++) out<<sqr(2*M_PI*iq/T/a)<<" "<<corr[iq]<<endl;
}

int main(int narg,char **arg)
{
  FILE *fin=open_file("input","r");
  read_formatted_from_file_expecting((char*)&T,fin,"%d","T");
  read_formatted_from_file_expecting((char*)&nconfs,fin,"%d","NConfs");
  cout<<"Read NConfs: "<<nconfs<<endl;
  njacks=sqrt(nconfs);
  int clust_size=nconfs/njacks;
  cout<<"Cluster size: "<<clust_size<<endl;
  cout<<"Njacks: "<<njacks<<endl;
  nconfs=njacks*clust_size;
  cout<<"Adjusted NConfs: "<<nconfs<<endl;
  read_formatted_from_file_expecting((char*)&nflavs,fin,"%d","NFlavs");
  read_formatted_from_file_expecting((char*)&nop,fin,"%d","NOp");
  read_formatted_from_file_expecting((char*)&a,fin,"%lg","a");
  
  //allocate
  int ncorr=icorr(nflavs-1,nop-1)+1;
  vector<jvec> corr(ncorr);
  for(int i=0;i<ncorr;i++) corr[i]=0*jvec(T,njacks);
  
  //load
  std::ifstream in("gm2");
  for(int iconf=0;iconf<nconfs;iconf++)
    for(int iop=0;iop<nop;iop++)
      for(int iflav=0;iflav<nflavs;iflav++)
	{
	  ignore_until_eol(in);
	  for(int jt=0;jt<T;jt++)
	    {
	      int t;
	      double r,i;
	      if(!(in>>t>>r>>i)) crash("reading iconf %d , iop %d , flav %d , jt %d, eof: %d",iconf,iop,iflav,jt,in.eof());
	      if(t!=jt) crash("iconf %d expected %d obtained %d",iconf,jt,t);
	      //cout<<iconf<<" "<<iop<<" "<<jflav<<" "<<t<<" "<<r<<" "<<i<<endl;
	      
	      int ijack=iconf/clust_size;
	      corr[icorr(iflav,iop)][t][ijack]+=r;
	    }
	  ignore_until_eol(in);
	  ignore_until_eol(in);
	}
  
  jvec AP(nflavs*nop,njacks),MP(nflavs*nop,njacks),AM(nflavs*nop,njacks),MM(nflavs*nop,njacks);
  for(int iop=0;iop<nop;iop++)
    for(int iflav=0;iflav<nflavs;iflav++)
      {
	int ind=icorr(iflav,iop);
	corr[ind].clusterize(clust_size);
	corr[ind]=corr[ind].simmetrized(1);
	corr[ind].print_to_file(combine("plots/corr_flv%d_op%d.xmg",iflav,iop).c_str());
	effective_mass(corr[ind],T/2,1,2).print_to_file(combine("plots/eff_flv%d_op%d.xmg",iflav,iop).c_str());
      }

  for(int iflav=0;iflav<nflavs;iflav++)
    for(int iop=1;iop<nop;iop++)
      {
	double coef=1;
	if(iop==2) coef=-1;
	write_HVP(combine("plots/gm2_flv%d_op%d.xmg",iflav,iop).c_str(),calc_HVP(coef*corr[icorr(iflav,iop)]));
      }
  
  return 0;
}

