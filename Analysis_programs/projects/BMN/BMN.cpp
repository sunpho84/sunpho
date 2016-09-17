#include "include.h"

const int njacks=16;

//data
vector<double> px;
jvec py;
vector<double> pe;
TMatrixD cova,inv_cova;

double dtmeas;
int N;
size_t ixmin_fit=1;
int ixmax;
double xmax=-1e300;
double ymin=1e300,ymax=-1e300;

int nppm=5;
int nmodes=1,nmodes_used;

//get a single par
template <class T> void get_par(T &out,string path,T def_val)
{
  ifstream input(path);
  if(!(input>>out)) out=def_val;
}

//load a single file
jack load_single_obs(const char *path)
{
  int clust_size;
  int ijack,it;
  double y;
  jack obs(njacks);
  FILE *fin_obs=open_file(path,"r");
  while(fscanf(fin_obs,"%d %d %d %lg",&ijack,&it,&clust_size,&y)==4)
    obs[ijack]=y;
  obs.clusterize(clust_size);
  return obs;
}

//load all the data
void load()
{
  get_par(dtmeas,"dt",0.1);
  get_par(N,"../../N",10);
  
  //open the file and load
  FILE *fin=open_file("obs_sq_Y_trace_ch2","r");
  map<int,array<double,njacks>> temp;
  int clust_size;
  {
    int ijack,it;
    double y;
    int each;
    get_par(each,"../each",3);
    cout<<"Each: "<<each<<endl;
    while(fscanf(fin,"%d %d %d %lg",&ijack,&it,&clust_size,&y)==4)
      if(it%each==0) temp[it][ijack]=y;
  }
  
  load_single_obs("obs_energy").write_to_file("energy");

  if(file_exists("obs_sq_Y_trace_ch_modulo"))
  {
    jack mod=load_single_obs("obs_sq_Y_trace_ch_modulo");
    mod.write_to_file("mod");
    jack sing=load_single_obs("obs_sq_Y_trace");
    sing.write_to_file("sing");
    cout<<(sqr(N)-1)*mod/(sqr(sing)/36)/8<<endl;
  }
  
  //resize using temp
  px.resize(temp.size());
  pe.resize(temp.size());
  py=jvec(temp.size(),njacks);
  
  //load all temp
  int it=0;
  double x0;
  for(auto d : temp)
    {
      for(int ijack=0;ijack<njacks;ijack++)
	{
	  px[it]=d.first*dtmeas;
	  py[it][ijack]=d.second[ijack];
	}
      py[it].clusterize(clust_size);
      pe[it]=py[it].err();
      
      //shift all times
      if(it==0) x0=px[it];
      px[it]-=x0;
      
      //take maximal values
      xmax=max(xmax,px[it]);
      ymin=min(ymin,py[it].med());
      ymax=max(ymax,py[it].med());
      
      it++;
    }
  
  //find xmax such that y is 2% of ymax
  ixmax=px.size()-1;
  if(file_exists("xmax"))
    {
      ifstream xmax_file("xmax");
      xmax_file>>xmax;
    }
  else
    {
      while(ixmax>=0 && fabs(py[ixmax].med()/ymax)<0.01) ixmax--;
      xmax=px[ixmax];
      ofstream xmax_file("xmax");
      xmax_file<<xmax<<endl;
    }
}

//fit function ansatz
double fun(double x,double *pars)
{
  double out=pars[0];
  for(int imode=0;imode<nmodes;imode++)
    {
      double *p=pars+nppm*imode;
            out+=p[2]*exp(-x*p[3]*p[5])*(cos((x-p[4])*p[5])+p[1]);
    }
  return out;
}

jack fun(double x,jvec p)
{
  jack out=p[0];
  for(int imode=0;imode<nmodes;imode++)
  out+=p[nppm*imode+2]*exp(-x*p[nppm*imode+3]*p[nppm*imode+5])*(cos((x-p[nppm*imode+4])*p[nppm*imode+5])+p[nppm*imode+1]);
  return out;
}

//fit data
int ijack_fit=0;
int ndof;
vector<double> diff(px.size());
int use_corr=0;
void uncorr_ch2(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  ndof=-nppm*nmodes+1;
  
  //compute and add to the ch2 if uncorr
  for(size_t i=ixmin_fit;i<px.size();i++)
    if(px[i]<xmax)
    {
      ndof++;
      double num=py[i][ijack_fit];
      double teo=fun(px[i],p);
      diff[i]=num-teo;
      ch+=sqr(diff[i]/pe[i]);
      //cout<<"uncorr "<<i<<1/sqr(py[i].err())<<" "<<inv_cova(i,i)<<endl;
    }
}

//compute the correlated ch2
void corr_ch2(int &npar,double *fuf,double &ch,double *p,int flag)
{
  uncorr_ch2(npar,fuf,ch,p,flag);
  
  //compute ch2 if correlated
  ch=0;
  for(size_t i=ixmin_fit;i<px.size();i++)
    for(size_t j=ixmin_fit;j<px.size();j++)
      {
	ch+=diff[i]*diff[j]*inv_cova(i,j);
	//cout<<"corr "<<i<<" "<<j<<" "<<cova(i,j)<<endl;
      }
}

//wrapper for jack
jack ch2(jvec &p,bool is_corr)
{
  jack ch(njacks);
  double pj[p.nel];
  
  for(ijack_fit=0;ijack_fit<=njacks;ijack_fit++)
    {
      for(int ip=0;ip<p.nel;ip++) pj[ip]=p[ip][ijack_fit];
      double *fuf=NULL;
      int flag=0;
      if(is_corr) corr_ch2(p.nel,fuf,ch[ijack_fit],pj,flag);
      else uncorr_ch2(p.nel,fuf,ch[ijack_fit],pj,flag);
    }
  return ch;
}

//super-wrappers
jack corr_ch2(jvec &p){return ch2(p,true);}
jack uncorr_ch2(jvec &p){return ch2(p,false);}

TMinuit minu;
double get_par(int i)
{
  double out,dum;
  minu.GetParameter(i,out,dum);
  return out;
}
void define_new_mode()
{
  cout<<"========================= adding new mode ("<<nmodes_used+1<<") ======================"<<endl;
  
  minu.DefineParameter(1+nppm*nmodes_used,"ConstDamp",(nmodes_used==0)?ymax*10:0.001,0.001,0,0);
  minu.DefineParameter(2+nppm*nmodes_used,"OscilDamp",(nmodes_used==0)?ymax:0.001,0.001,0,0);
  minu.DefineParameter(4+nppm*nmodes_used,"Offset-X",0,0.001,-100,100);
  if(nmodes_used==0)
    {
      minu.DefineParameter(3+nppm*nmodes_used,"W",0.1,0.001,0,100);
      minu.DefineParameter(5+nppm*nmodes_used,"wim",1.4,0.001,0,100);
    }
  else
    {
      minu.DefineParameter(3+nppm*nmodes_used,"W",get_par(3+nppm*(nmodes_used-1))*1.5,0.001,get_par(3+nppm*(nmodes_used-1)),100);
      minu.DefineParameter(5+nppm*nmodes_used,"wim",get_par(5+nppm*(nmodes_used-1))*1.5,0.001,get_par(5+nppm*(nmodes_used-1)),100);
    }
  
  nmodes_used++;
}

void setup_fit()
{
  //minu.SetPrintLevel(-1);
  minu.DefineParameter(0,"Y-Offset",0,0.001,0,0);
  
  //resize the covariance matrix
  diff.resize(px.size());
  cova.ResizeTo(px.size(),px.size());
  inv_cova.ResizeTo(cova);
  
  double tol=1e-16;
  int iflag;
  minu.mnexcm("SET ERR",&tol,1,iflag);
  minu.SetMaxIterations(10000000);
}

//compute the covariance matrix and its inverse
void prepare_cova()
{
  //compute covariance
  for(size_t i=0;i<px.size();i++)
    for(size_t j=i;j<px.size();j++)
      cova(j,i)=cova(i,j)=cov(py[i],py[j]);
  
  //put to 0 below a threshold
  for(size_t i=0;i<px.size();i++)
    for(size_t j=i+1,ex=false;j<px.size();j++)
      if(i<ixmin_fit||j<ixmin_fit||fabs(cova(i,j)/sqrt(cova(i,i)*cova(j,j)))<0.5||ex==true)
	{
	  cova(i,j)=cova(j,i)=0;
	  ex=true;
	}
  
  //invert
  inv_cova=cova;
  inv_cova.Invert();
  
  //check inversion
  if(0)
    {
      TMatrixD tempA=cova*inv_cova;
      TMatrixD tempB=inv_cova*cova;
      double diffA=0,diffB=0;
      for(size_t i=0;i<px.size();i++)
	for(size_t j=0;j<px.size();j++)
	  {
	    diffA+=sqr(tempA(i,j)-(i==j));
	    diffB+=sqr(tempB(i,j)-(i==j));
	  }
      cout<<"Diff of inverse: "<<diffA<<" "<<diffB<<endl;
    }
  
  //print the full matrix
  if(0)
    {
      cout<<"Matr"<<endl;
      for(size_t i=0;i<px.size();i++)
	for(size_t j=0;j<px.size();j++)
	  cout<<i<<" "<<j<<" "<<cova(i,j)<<endl;
    }
  
  //print the diagonal
  if(0)
    {
      cout<<"Diag"<<endl;
      for(size_t i=0;i<px.size();i++)
	cout<<i<<" "<<inv_cova(i,i)<<" "<<1/sqr(py[i].err())<<endl;
    }
}

//perform the fit
jvec fit()
{
  if(use_corr)
    {
      minu.SetFCN(corr_ch2);
      prepare_cova();
    }
  else         minu.SetFCN(uncorr_ch2);
  
  //the first time make a monte carlo minimization
  static int first_time=1;
  if(first_time)
    for(int imode=1;imode<=nmodes;imode++)
      {
	define_new_mode();
	minu.Migrad();
      }
  first_time=0;
  
  //loop over jacknives
  jvec pars(1+nmodes_used*nppm,njacks);
  for(ijack_fit=0;ijack_fit<=njacks;ijack_fit++)
    {
      minu.Migrad();
      
      //get parameters
      for(int i=0;i<1+nmodes_used*nppm;i++)
	pars[i][ijack_fit]=get_par(i);
    }
  
  return pars;
}

//print fit
void write_plot(string path,jvec &pars)
{
  ofstream out(path);
  //fit line
  out<<"@type xy"<<endl;
  for(double x=px[ixmin_fit];x<xmax;x+=0.01)
    out<<x<<" "<<fun(x,pars).med()<<endl;
  out<<"&"<<endl<<"@type xydy"<<endl;
  
  //data
  out<<"@s1 line type 0"<<endl;
  for(size_t i=0;i<px.size();i++)
    out<<px[i]<<" "<<py[i]<<endl;
}

//print fit
void write_exc_plot(string path,jvec pars)
{
  ofstream out(path);
  //data
  //pars[3]=0;
  
  out<<"@s0 line type 0"<<endl;
  out<<"@type xydy"<<endl;
  for(size_t i=0;i<px.size();i++)
    out<<px[i]<<" "<<py[i]-fun(px[i],pars)<<endl;
}

//fit moving xmin in the passed interval
jvec fit_in_the_interval(string path="o.xmg",size_t ixmin_fit_min=10,size_t ixmin_fit_max=px.size()/3*2)
{
  jvec pars;
  ofstream ruch2("ruch2.xmg");
  ofstream ruW("ruW.xmg");
  ofstream ruD("ruD.xmg");
  ofstream ruwre("ruwre.xmg");
  ofstream ruwim("ruwim.xmg");
  ruW<<"@type xydy"<<endl;
  ruD<<"@type xydy"<<endl;
  ruwre<<"@type xydy"<<endl;
  ruwim<<"@type xydy"<<endl;
  for(ixmin_fit=ixmin_fit_min;ixmin_fit<=ixmin_fit_max;ixmin_fit++)
    if(px[ixmin_fit+nmodes_used*nppm+1]<xmax)
    {
      pars=fit();
      //check
      jack c=ch2(pars,use_corr);
      //converged=(c.med()-0*c.err()<=ndof);
      
      cout<<" Ch2 "<<c<<endl;
      
      ixmin_fit++;
      
      jack W=pars[3];
      jack D=pars[1];
      jack wim=pars[5];
      jack wre=W*wim;;
      if(!isnan(W.err())) ruW<<px[ixmin_fit]<<" "<<W<<endl;
      if(!isnan(D.err())) ruD<<px[ixmin_fit]<<" "<<D<<endl;
      if(!isnan(wre.err())) ruwre<<px[ixmin_fit]<<" "<<wre<<endl;
      if(!isnan(wim.err())) ruwim<<px[ixmin_fit]<<" "<<wim<<endl;
      ruch2<<px[ixmin_fit]<<" "<<c/ndof<<endl;
    }
  
  write_plot(path,pars);
  write_exc_plot("oexc.xmg",pars);
  
  cout<<"Ch2: "<<uncorr_ch2(pars)<<" / "<<ndof;
  if(use_corr) cout<<", corr: "<<corr_ch2(pars)<<endl;
  else cout<<endl;
  
  pars[0].write_to_file("yoff");
  pars[1].write_to_file("cons_damp");
  pars[2].write_to_file("osci_damp");
  pars[3].write_to_file("W");
  pars[4].write_to_file("xoff");
  pars[5].write_to_file("wim");
  (pars[3]*pars[5]).write_to_file("wre");

  cout<<"Fit interval: "<<px[ixmin_fit]<<" "<<xmax<<endl;
  
  return pars;
}

void fit_up_to_x(double xmin_fit_ext)
{
  //search the minimum
  size_t ixmin_fit_ext=10;
  while(px[ixmin_fit_ext]<xmin_fit_ext) ixmin_fit_ext++;
  
  fit_in_the_interval("o.xmg",10,ixmin_fit_ext);
}

int main(int narg,char **arg)
{
  load();
  setup_fit();
  
  ofstream out("obs_resc_N.xmg");
  out<<"@type xydy"<<endl;
  for(size_t i=0;i<px.size();i++) out<<px[i]/pow(N,0.25)<<" "<<py[i]<<endl;
  
  //fit iteratively
  if(narg<2) fit_in_the_interval();
  else fit_up_to_x(stod(arg[1]));
  
  return 0;
}
