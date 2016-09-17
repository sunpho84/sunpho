#include "include.h"

int nflavs;
int njacks;
int nconfs;
int nop;
int T;
vector<int> tmin,tmax;
double am;
const int tav=3;

const double mpi_phys=0.135;
const double mss_phys=0.692;
//const double meta_phys=0.547862;
const double metap_phys=0.958;
const double fpi_phys=0.1302;
const double fact=metap_phys/mss_phys;

//////////////////////////////////// to find phys point /////////////////////////////

int glb_tmin,glb_tmax;

void ignore_until_eol(ifstream &in)
{in.ignore(10000,'\n');}

double *c,*e;

template <class _T> _T fun_fit_nosymm_noosc(_T &A,_T &M,int t)
{
  return A*exp(-M*t);
}

template <class _T> _T fun_fit_nosymm(_T &AP,_T &MP,_T &AM,_T &MM,int t)
{
  return fun_fit_nosymm_noosc(AP,MP,t)/(2*MP)+pow(-1,t)*fun_fit_nosymm_noosc(AM,MM,t);
}

template <class _T> _T fun_fit(_T &AP,_T &MP,_T &AM,_T &MM,int t)
{
  return fun_fit_nosymm(AP,MP,AM,MM,t)+fun_fit_nosymm(AP,MP,AM,MM,T-t);
}

void ch2(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  
  double MP=p[1];
  double MM=p[3];
  
  for(int t=glb_tmin;t<=glb_tmax;t++)
    {
      int tp;
      if(t<=T/2) tp=t;
      else tp=T-t;
      double teo=fun_fit(p[0],MP,p[2],MM,t);
      double cont=sqr((c[tp]-teo)/e[tp]);
      ch+=cont;
    }
}

void fit(jack &AP,jack &MP,jack &AM,jack &MM,jvec corr,int tmin,int tmax)
{
  glb_tmin=tmin;
  glb_tmax=tmax;
  MP=MM=0;
  
  TMinuit minu;
  minu.SetPrintLevel(-1);
  minu.SetFCN(ch2);
  //set tolerance
  double tol=1e-16;
  int iflag;
  minu.SetMaxIterations(10000000);
  minu.mnexcm("SET ERR",&tol,1,iflag);
  
  minu.DefineParameter(0,"AP",0,0.001,0,0);
  minu.DefineParameter(1,"MP",0.2,0.001,0,1);
  minu.DefineParameter(2,"AM",0,0.001,0,0);
  minu.DefineParameter(3,"MM",0,0.001,0,0);
  minu.FixParameter(3);
  for(int iel=0;iel<=T/2;iel++) e[iel]=corr[iel].err();
  
  for(int ijack=0;ijack<=corr.njack;ijack++)
    {
      //if(ijack==0) minu.SetPrintLevel(-1);
      //else         minu.SetPrintLevel(1);
      for(int iel=0;iel<=T/2;iel++) c[iel]=corr[iel][ijack];
      minu.Migrad();
      double dum;
      minu.GetParameter(0,AP[ijack],dum);
      minu.GetParameter(1,MP[ijack],dum);
      minu.GetParameter(2,AM[ijack],dum);
      minu.GetParameter(3,MM[ijack],dum);
    }
}

void fit(jack &AP,jack &MP,jack &AM,jack &MM,jvec corr,int tmin,int tmax,ofstream &out)
{
  fit(AP,MP,AM,MM,corr,tmin,tmax);
  for(int t=glb_tmin;t<=glb_tmax;t++) out<<t<<" "<<fun_fit(AP,MP,AM,MM,t)<<endl;
  out<<"&"<<endl;
}

int icorr(int iflav,int iop)
{return iop+nop*iflav;}

template<class T> T correction_without_f(T m5,T mss)
{
  T num=sqr(m5)/8;
  T den=1+sqr(m5)/(2*sqr(mss))+3*sqr(m5)/(2*sqr(mss*fact));
  return num/den;
}

template<class T> T correction(T fpi0,T m5,T mss)
{return correction_without_f(m5,mss)*sqr(fpi0);}

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
  
  //allocate
  int ncorr=icorr(nflavs-1,nop-1)+1;
  vector<jvec> corr(ncorr);
  tmin.resize(ncorr);
  tmax.resize(ncorr);
  for(int i=0;i<ncorr;i++) corr[i]=0*jvec(T,njacks);
  
  for(int iflav=0;iflav<nflavs;iflav++)
    for(int iop=0;iop<nop;iop++)
      {
	expect_string_from_file(fin,combine("TInt_flv%d_op%d",iflav,iop).c_str());
	int ind=icorr(iflav,iop);
	read_formatted_from_file((char*)&tmin[ind],fin,"%d","TMin");
	read_formatted_from_file((char*)&tmax[ind],fin,"%d","TMax");
      }
  
  read_formatted_from_file_expecting((char*)&am,fin,"%lg","am");
  cout<<"am: "<<am<<endl;
  
  c=new double[T];
  e=new double[T];
  
  //load
  std::ifstream in("stag_mesons");
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
  for(int iflav=0;iflav<nflavs;iflav++)
    {
      cout<<endl<<endl<<"================ flav "<<iflav<<" ================"<<endl;
      
      //print out
      ofstream corr_out(combine("corr_flv%d.xmg",iflav).c_str());
      corr_out<<"@type xydy"<<endl;
      for(int iop=0;iop<nop;iop++)
	{
	  int ind=icorr(iflav,iop);
	  corr[ind].clusterize(clust_size);
	  corr[ind]=corr[ind].simmetrized(1);
	  corr_out<<corr[ind]<<"&"<<endl;
	}
      
      //print effective mass
      for(int iop=0;iop<nop;iop++)
	{
	  int ind=icorr(iflav,iop);
	  constant_fit(effective_mass(corr[icorr(iflav,iop)],T/2,1,2),tmin[ind],tmax[ind],combine("effmass_flv%d_op%d.xmg",iflav,iop).c_str());
	}
      
      for(int iop=0;iop<nop;iop++)
	{
	  int ind=icorr(iflav,iop);
	  fit(AP[ind],MP[ind],AM[ind],MM[ind],corr[ind],tmin[ind],tmax[ind]);
	  ofstream outp(combine("fit_flv%d_op%d.xmg",iflav,iop).c_str());
	  outp<<"@type xydy"<<endl;
	  for(int t=tmin[ind];t<tmax[ind];t++) outp<<t<<" "<<(fun_fit_nosymm_noosc(AP[ind],MP[ind],t)+fun_fit_nosymm_noosc(AP[ind],MP[ind],T-t))/(2*MP[ind])<<endl;
	  outp<<"&"<<endl;
	  for(int t=0;t<T/2;t++) outp<<t<<" "<<corr[ind][t]<<endl;
	  
	  cout<<endl<<"======iop"<<iop<<"======"<<endl;
	  cout<<"AP "<<smart_print(AP[ind])<<endl;
	  cout<<"MP "<<smart_print(MP[ind])<<endl;
	  cout<<"AM "<<smart_print(AM[ind])<<endl;
	  cout<<"MM "<<smart_print(MM[ind])<<endl;
	}
    }
  
  jack fpi0=2*am*sqrt(AP[0])/(2*MP[0]*MP[0]);
  jack mpi0=MP[icorr(0,0)];
  jack mpi5=MP[icorr(0,1)];
  jack mss0=MP[icorr(1,0)];
  jack mss5=MP[icorr(1,1)];
  
  jack st_chpt_pred=correction(fpi0,mpi5,mss5);
  double st_chpt_pred_phys=correction(fpi_phys,mpi_phys,mss_phys);
  
  jack st_chpt_pred_without_f=correction_without_f(mpi5,mss5);
  double st_chpt_pred_phys_without_f=correction_without_f(mpi_phys,mss_phys);
  
  cout<<"==================================="<<endl;
  cout<<endl;
  cout<<" fpi: "<<smart_print(fpi0)<<endl;
  cout<<" mpi*T: "<<smart_print(mpi0*T)<<endl;
  cout<<" fpi/mpi: "<<smart_print(fpi0/mpi0)<<" phys one: "<<fpi_phys/mpi_phys<<endl;
  cout<<" M_eta_ss/M_pi: "<<smart_print(MP[icorr(1,0)]/MP[icorr(0,0)])<<" phys one: "<<mss_phys/mpi_phys<<endl;
  cout<<" M_5/M_pi: "<<mpi5/MP[0]<<endl;
  cout<<" sqrt(M_pi/M_5): "<<sqrt(MP[0]/mpi5)<<endl;
  cout<<" correction (M_pi/M_5)^2: "<<sqr(MP[0]/mpi5)<<endl;
  cout<<" st_chpt_pred: "<<smart_print(st_chpt_pred)<<endl;
  cout<<" physical value: "<<st_chpt_pred_phys<<", fourth root: "<<pow(st_chpt_pred_phys,0.25)<<endl;
  cout<<" correcting factor for fourth root: "<<pow(st_chpt_pred_phys/st_chpt_pred,0.25)<<endl;
  cout<<" correcting factor for susc: "<<st_chpt_pred_phys/st_chpt_pred<<endl;
  cout<<" correcting factor for susc (without f, units a^-2): "<<st_chpt_pred_phys_without_f/st_chpt_pred_without_f<<endl;
  cout<<" correcting factor for fourth root susc (without f, only latt, phys units): "<<pow(correction_without_f(mpi0,mss0)/correction_without_f(mpi5,mss5),0.25)<<endl;
  
  delete[] c;
  delete[] e;
  
  return 0;
}

