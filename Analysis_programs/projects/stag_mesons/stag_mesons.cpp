#include "include.h"

int nflavs;
int nconfs;
int nop;
int T;
int tmin,tmax;
double am;
const int tav=3;

const char op_tag[6][100]={
  "\\gamma_5 \\times \\gamma_5 (0 link)",
  "\\gamma_5 \\times \\gamma_1\\gamma_5 (2 links)",
  "\\gamma_5 \\times \\gamma_0\\gamma_1 (2 links)",
  "\\gamma_5 \\times \\gamma_0 (3 links)",
  "\\gamma_0 \\gamma_5 \\times 1 (3 links)",
  ""};

//////////////////////////////////// to find phys point /////////////////////////////

double db0=11.972809241766459,f0=0.26772510275074468;
double r0=2.1975;
double l3=1.0423672088359199,l4=2.1490196003108331;
double xi=db0/(16*M_PI*M_PI*f0*f0);
double ml_phys=3.6e-3;

template <class T> T m_fun(T m)
{
  m*=r0;
  
  T m2=db0*m*(1+  xi*m*(log(db0*m)+l3));
  
  return sqrt(m2)/r0;
}

template <class T> T f_fun(T m)
{
  m*=r0;
  
  T f=   f0* (1-2*xi*m*(log(db0*m)-l4));
  
  return f/r0;
}

template <class T> T ratio(T m)
{return m_fun(m)/f_fun(m);}


double find(double x)
{
  double m=1e-5;
  double s=1e-5;
  
  do
    if(ratio(m+s)<=x)
      {
	m+=s;
	s*=2;
      }
    else s/=2;
  while(s>1e-14);
  
  return m;
}

jack find(jack x)
{
  jack out(nconfs);
  for(int iconf=0;iconf<=nconfs;iconf++) out[iconf]=find(x[iconf]);
  return out;
}

///////////////////////////////////////////////////

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

void print_data(const char *path,vector<jvec> &data)
{
  ofstream out(path);
  int s=0;
  out<<"@type xydy"<<endl;
  for(int iop=0;iop<nop;iop++)
    out<<data[iop]<<"@s"<<s++<<" legend \"$$"<<op_tag[iop]<<"\"\n&"<<endl;
}

int main(int narg,char **arg)
{
  FILE *fin=open_file("input","r");
  read_formatted_from_file_expecting((char*)&T,fin,"%d","T");
  read_formatted_from_file_expecting((char*)&tmin,fin,"%d","TMin");
  read_formatted_from_file_expecting((char*)&tmax,fin,"%d","TMax");
  read_formatted_from_file_expecting((char*)&nconfs,fin,"%d","NConfs");
  read_formatted_from_file_expecting((char*)&nflavs,fin,"%d","NFlavs");
  read_formatted_from_file_expecting((char*)&nop,fin,"%d","NOp");
  read_formatted_from_file_expecting((char*)&am,fin,"%lg","am");
  cout<<"am: "<<am<<endl;
  
  c=new double[T];
  e=new double[T];
  
  for(int iflav=0;iflav<nflavs;iflav++)
    {
      vector<jvec> corr(nop);
      for(int iop=0;iop<nop;iop++) corr[iop]=jvec(T,nconfs);
      std::ifstream in("stag_mesons");
      if(nflavs>1) cout<<endl<<endl<<"================ flav "<<iflav<<" ================"<<endl;
      
      //load
      for(int iconf=0;iconf<nconfs;iconf++)
	{
	  for(int iop=0;iop<nop;iop++)
	    for(int jflav=0;jflav<nflavs;jflav++)
	      {
		ignore_until_eol(in);
		for(int jt=0;jt<T;jt++)
		  {
		    int t;
		    double r,i;
		    if(!(in>>t>>r>>i)) crash("reading iconf %d , iop %d , flav %d , jt %d",iconf,iop,jflav,jt);
		    if(t!=jt) crash("iconf %d expected %d obtained %d",iconf,jt,t);
		    //cout<<iconf<<" "<<iop<<" "<<jflav<<" "<<t<<" "<<r<<" "<<i<<endl;
		    
		    if(iflav==jflav) corr[iop][t][iconf]=r;
		  }
		ignore_until_eol(in);
		ignore_until_eol(in);
	      }
	}
      
      ofstream corr_out(combine("corr_flv%d.xmg",iflav).c_str());
      corr_out<<"@type xydy"<<endl;
      for(int iop=0;iop<nop;iop++)
	{
	  corr[iop].clusterize();
	  corr[iop]=corr[iop].simmetrized(1);
	  corr_out<<corr[iop]<<"&"<<endl;
	}
      
      double mpi_phys=0.135;
      double fpi_phys=0.1302;
      
      jvec fm(nop,nconfs);
      for(int iop=0;iop<nop;iop++)
	{
	  int dt=2;
	  fm[iop]=constant_fit(effective_mass(corr[iop],T/2,1,dt),tmin,T/2,combine("effmass%d_flv%d.xmg",iop,iflav).c_str());
	}
      vector<jvec> AP_corr(nop);
      vector<jvec> MP_corr(nop);
      vector<jvec> AM_corr(nop);
      vector<jvec> MM_corr(nop);
      int s=0;
      for(int iop=0;iop<nop;iop++)
	{
	  AP_corr[iop]=AM_corr[iop]=MP_corr[iop]=MM_corr[iop]=jvec(T/2+1,nconfs);
	  ofstream out(combine("corr%d_flv%d_fit.xmg",iop,iflav).c_str());
	  out<<"@type xydy"<<endl;
	  for(int tf=0;tf<=T/2;tf++) fit(AP_corr[iop][tf],MP_corr[iop][tf],AM_corr[iop][tf],MM_corr[iop][tf],corr[iop],tf,tf+tav,out);
	}
      
      ofstream out(combine("ops_flv%d.xmg",iflav).c_str());
      out<<"@type xydy"<<endl;
      for(int iop=0;iop<nop;iop++)
	out<<MP_corr[iop]<<"@s"<<s++<<" legend \"$$"<<op_tag[iop]<<"\"\n&"<<endl;
      
      print_data(combine("AP_flv%d.xmg",iflav).c_str(),AP_corr);
      print_data(combine("MP_flv%d.xmg",iflav).c_str(),MP_corr);
      print_data(combine("AM_flv%d.xmg",iflav).c_str(),AM_corr);
      print_data(combine("MM_flv%d.xmg",iflav).c_str(),MM_corr);
      
      jvec AP(nop,nconfs),MP(nop,nconfs),AM(nop,nconfs),MM(nop,nconfs);
      for(int iop=0;iop<nop;iop++)
	{
	  fit(AP[iop],MP[iop],AM[iop],MM[iop],corr[iop],tmin,T/2-1);
	  ofstream outp(combine("/tmp/p%d.xmg",iop).c_str());
	  outp<<"@type xydy"<<endl;
	  for(int t=tmin;t<T/2;t++) outp<<t<<" "<<(fun_fit_nosymm_noosc(AP[iop],MP[iop],t)+fun_fit_nosymm_noosc(AP[iop],MP[iop],T-t))/(2*MP[iop])<<endl;
	  outp<<"&"<<endl;
	  for(int t=0;t<T/2;t++) outp<<t<<" "<<corr[iop][t]<<endl;
	  
	  ofstream outm(combine("/tmp/m%d.xmg",iop).c_str());
	  outm<<"@type xydy"<<endl;
	  for(int t=tmin;t<T/2;t++) outm<<t<<" "<<pow(-1,t)*(corr[iop][t]-fun_fit_nosymm_noosc(AP[iop],MP[iop],t)-
							     fun_fit_nosymm_noosc(AP[iop],MP[iop],T-t))/(2*MP[iop])<<endl;
	  cout<<endl<<"======iop"<<iop<<"======"<<endl;
	  cout<<"AP "<<smart_print(AP[iop])<<endl;
	  cout<<"MP "<<smart_print(MP[iop])<<endl;
	  cout<<"AM "<<smart_print(AM[iop])<<endl;
	  cout<<"MM "<<smart_print(MM[iop])<<endl;
	  
	  if(iop==0 && iflav==0)
	    {
	      jack fpi_bis=2*am*sqrt(AP[iop])/(2*MP[iop]*MP[iop]);
	      jack rat=MP[iop]/fpi_bis;
	      cout<<"f: "<<smart_print(fpi_bis)<<endl;
	      cout<<"M/f, ops: "<<smart_print(rat)<<" exp: "<<mpi_phys/fpi_phys<<endl;
	      jack ml_reno=find(rat);
	      jack agevmuno=MP[iop]/m_fun(ml_reno);
	      cout<<"Renormalized ml: "<<smart_print(ml_reno)<<" = "<<smart_print(ml_reno/ml_phys)<<" physical mass"<<endl;
	      cout<<"Renormalization constant: "<<smart_print(ml_reno/(am/agevmuno))<<endl;
	      cout<<"Lattice spacing: "<<smart_print(0.197*agevmuno)<<" fm, "<<smart_print(agevmuno)<<" GeV^-1; inverse: "<<smart_print(1/agevmuno)<<" GeV"<<endl;
	    }
	  cout<<"ratio with taste 0: "<<smart_print(MP[iop]/MP[0])<<endl;
	}
	
      if(nop==5)
	{
	  double deg[5]={1,4,6,4,1};
	  jack invave(nconfs),ave(nconfs);
	  invave=ave=0;
	  for(int iop=0;iop<5;iop++)
	    {
	      ave+=deg[iop]*pow(MP[iop],2);
	      invave+=deg[iop]*pow(MP[iop],-2);
	    }
	  invave=1/sqrt(invave/16);
	  ave=sqrt(ave/16);
	  cout<<"inverse average: "<<smart_print(invave)<<endl;
	  cout<<"average: "<<smart_print(ave)<<endl;
	}
      
    }
  
  delete[] c;
  delete[] e;
  
  return 0;
}


