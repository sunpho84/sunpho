#include <algorithm>
#include <complex>
#include <iostream>
#include <fstream>
#include <random>

#define CRASH(...) internal_crash(__LINE__,__FILE__,__VA_ARGS__)

using namespace std;

//types and constants
#define ndims 2
typedef complex<double> dcomplex;
typedef int coords[ndims];

//random number generators
random_device *rd;
mt19937_64 *gen;
uniform_int_distribution<unsigned long long> *dis;

//parametes
int N=2;
double beta=1.1;
int L=36;
double g=1/(N*beta);

//geometry
int V;
int *neigh_data;

//system
dcomplex *zeta_data,*lambda_data;

///////////////////////////////////////////////////////////////////////////////////////////////

//crash
void internal_crash(int line,const char *file,const char *templ,...)
{
  //expand error message
  char mess[1024];
  va_list ap;
  va_start(ap,templ);
  vsprintf(mess,templ,ap);
  va_end(ap);
  
  fprintf(stderr,"ERROR on line %d of file \"%s\", message error: \"%s\".\n",line,file,mess);
  exit(1);
}

//compute the square
template <class T> T sqr(T in)
{return in*in;}

//return zeta and lambda
inline dcomplex *zeta(int site)
{return zeta_data+site*N;}
inline dcomplex *lambda(int site)
{return lambda_data+site*ndims;}

//return a scalar product between two zetas
double get_zeta_scalprod(dcomplex *a,dcomplex *b)
{
  double res=0;
  for(int n=0;n<N;n++) res+=(conj(b[n])*a[n]).real();
  return res;
}

//return the norm of a zeta
double get_zeta_norm(dcomplex *z)
{
  double norm2=0;
  for(int n=0;n<N;n++) norm2+=norm(z[n]);
  return sqrt(norm2);
}

//reunitarize a zeta
void zeta_unitarize(dcomplex *z)
{
  double zeta_norm_inv=1/get_zeta_norm(z);
  for(int n=0;n<N;n++) z[n]*=zeta_norm_inv;
}

//return the deviation from unitarity of a zeta
inline double check_zeta_unitarity(dcomplex *z)
{return fabs(get_zeta_norm(z)-1);}

//return the result of the scalar product of two lambda
double get_lambda_scalprod(dcomplex &a,dcomplex &b)
{return (conj(a)*b).real();}

//return the norm of a lambda
double get_lambda_norm(dcomplex &l)
{return sqrt(norm(l));}

//reunitarize a lambda
double lambda_unitarize(dcomplex &l)
{
  double n=get_lambda_norm(l);
  l*=1/n;
  return n;
}

//return the deviation from unitarity of a lambda
inline double check_lambda_unitarity(dcomplex &l)
{return fabs(get_lambda_norm(l)-1);}

//get coords of site
void coords_of_site(coords c,int site)
{
  for(int mu=ndims-1;mu>=0;mu--)
    {
      c[mu]=site%L;
      site/=L;
    }
}

//return the coordinate of a site
int site_of_coords(coords c)
{
  int site=0;
  for(int mu=0;mu<ndims;mu++) site=site*L+c[mu];
  return site;
}

//return the neighbors
inline int &neighdw(int site,int mu)
{return neigh_data[0+2*(mu+ndims*site)];}
inline int &neighup(int site,int mu)
{return neigh_data[1+2*(mu+ndims*site)];}

//return a double between [0,max)
double get_unif_double(double max,bool incl=false)
{
  double res;
  do res=max*(*dis)(*gen)/dis->b();
  while((!incl)&&(res==max));

  return res;
}

//compute p(theta)
inline double fun_ptheta(double theta,double a,int k)
{
  //compute the factor
  double f=sin(theta);
  f*=f;

  //compute first piece
  double p=1;
  for(int i=0;i<k-1;i++) p*=f;
  
  return p*exp(a*cos(theta));
}

//obtain theta
double get_theta(double a,int k)
{
  //compute parameters
  double zita=(k-1)/a;
  double theta0=acos(sqrt(1+zita*zita)-zita);
  double ctheta0=cos(theta0);
  double c=sqrt(2*(k-1)*(1-zita*ctheta0)/(1-ctheta0*ctheta0));
  double ptheta0=fun_ptheta(theta0,a,k);
  
  //extract theta
  double theta;
  double eta=0.99;
  bool acc;
  int no=0;
  do
    {
      double chi=get_unif_double(1,true);
      theta=theta0+tan(chi*atan(c*(M_PI-theta0))+(chi-1)*atan(c*theta0))/c;

      //reweighting
      double ptheta=fun_ptheta(theta,a,k);
      double pacc=ptheta/ptheta0*(1+sqr(c*(theta-theta0)))*eta;
      double extr=get_unif_double(1);
      if(pacc>1) CRASH("pacc: %lg",pacc);
      //cerr<<extr<<" "<<theta<<" "<<theta0<<" "<<zita<<endl;
      acc=(extr<pacc);
      if(!acc) no++;
      if(no>1000) CRASH("k: %lg, a: %d, zita: %lg, ptheta: %lg, ptheta0: %lg, c: %lg, theta: %lg, theta0: %lg, pacc: %lg",
			k,a,zita,ptheta,ptheta0,c,theta,theta0,pacc);
    }
  while(!acc);
  
  return theta;
}

//taken by appendix C of hep-lat/9210016
//correction done: h is returned in place of (wrong) g
double get_theta_1(const double a)
{
  double eps=0.001;
  double as=0.798953686083986;
  double dap=max(0.0,a-as);
  double del=0.35*dap+1.03*sqrt(dap);
  double alp=min(sqrt(a*(2-eps)),max(sqrt(eps*a),del));
  double bet=max(alp*alp/a,(cosh(M_PI*alp)-1)/(exp(2*a)-1))-1;
  double bt1=sqrt((1+bet)/(1-bet));
  
  double h; //result
  bool acc; //accepted or not
  int no=0;
  do
    {
      double r=get_unif_double(1);
      double h1=bt1*tan((2*r-1)*atan(tanh(M_PI*alp/2)/bt1));
      h=log((1+h1)/(1-h1))/alp;
      
      //decite if accept or reject
      double g=exp(-a*(1-cos(h)))*(cosh(alp*h)+bet)/(1+bet);
      double p=get_unif_double(1);
      acc=(p<g);
      if(!acc) no++;
      if(no>1000) CRASH("%d",no);
    }
  while(!acc);
  
  return h;
}

//set an U1 to random
void set_U1_to_rnd(dcomplex &U)
{
  //extract a phase
  double ph=get_unif_double(2*M_PI);
  U=dcomplex(cos(ph),sin(ph));
}

//set an O(N) to random respecting the bound of unitarity
void set_ON_to_rnd(dcomplex *O)
{
  //first of all extract their norm in such: ordering
  double w[N+1];
  w[0]=0;
  for(int i=1;i<N;i++) w[i]=get_unif_double(1);
  w[N]=1;
  sort(w,w+N+1);
  
  //extract a random complex number for each site using the extracted norm
  for(int i=0;i<N;i++)
    {
      double nor=sqrt(w[i+1]-w[i]);
      double the=get_unif_double(2*M_PI);
      O[i]=nor*dcomplex(cos(the),sin(the));
    }
}

//initialize the system
void init_system_to_hot()
{
  for(int site=0;site<V;site++)
    {
      //fill the lambda
      for(int mu=0;mu<ndims;mu++)
	set_U1_to_rnd(lambda(site)[mu]);
      
      //fill the Zeta
      set_ON_to_rnd(zeta(site));
    }
}

//control that each element is unitary
void check_system_unitarity(double res=1.e-13)
{
  //check all sites
  for(int site=0;site<V;site++)
    {
      double dev_zeta=check_zeta_unitarity(zeta(site));
      if(dev_zeta>res) CRASH("zeta norm for site %d deviates from 1 by %lg",site,dev_zeta);

      //check all links
      for(int mu=0;mu<ndims;mu++)
	{
	  double dev_lambda=check_lambda_unitarity(lambda(site)[mu]);
	  if(dev_lambda>res) CRASH("lambda norm for site %d mu %d deviates from 1 by %lg",site,mu,dev_lambda);
	}
    }
}

//compute the total energy or action
double energy()
{
  double res=0;
  for(int site=0;site<V;site++)
    for(int mu=0;mu<ndims;mu++)
      {
	int site_up=neighup(site,mu);
	for(int n=0;n<N;n++) res+=(conj(zeta(site_up)[n])*zeta(site)[n]*lambda(site)[mu]).real();
	}
  
  return -(2*res-2*V*ndims);
}
double action()
{return energy()/g;}

//return the energy/action of a single site
//NB: the total action will be half the sum of the energy of all sites!
double site_energy(int site)
{
  double res=0;
  for(int mu=0;mu<ndims;mu++)
    {
      int site_up=neighup(site,mu);
      for(int n=0;n<N;n++) res+=(conj(zeta(site_up)[n])*zeta(site)[n]*lambda(site)[mu]).real();
      int site_dw=neighdw(site,mu);
      for(int n=0;n<N;n++) res+=(conj(zeta(site_dw)[n])*zeta(site)[n]*conj(lambda(site_dw)[mu])).real();
    }
  
  return -(2*res-4*ndims);
}
inline double site_action(int site)
{return site_energy(site)/g;}

//return the energy/action of a single link
double link_energy(int site,int mu)
{
  double res=0;
  
  int site_up=neighup(site,mu);
  for(int n=0;n<N;n++) res+=(conj(zeta(site_up)[n])*zeta(site)[n]*lambda(site)[mu]).real();
  
  return -(2*res-2);
}
inline double link_action(int site,int mu)
{return link_energy(site,mu)/g;}

//compute the staple of zeta
void site_staple(dcomplex *staple,int site)
{
  for(int n=0;n<N;n++) staple[n]=0;

  for(int mu=0;mu<ndims;mu++)
    {
      int site_up=neighup(site,mu);
      for(int n=0;n<N;n++) staple[n]+=2.0*zeta(site_up)[n]*conj(lambda(site)[mu]);
      int site_dw=neighdw(site,mu);
      for(int n=0;n<N;n++) staple[n]+=2.0*zeta(site_dw)[n]*lambda(site_dw)[mu];
    }
}
inline double site_staple_energy(int site,dcomplex *staple)
{
  double res=0;
  for(int n=0;n<N;n++) res+=(staple[n]*conj(zeta(site)[n])).real();
  return -(res-4*ndims);
}
inline double site_staple_action(int site,dcomplex *staple)
{return site_staple_energy(site,staple)/g;}

//compute the staple of lambda
void link_staple(dcomplex &staple,int site,int mu)
{
  staple=0;
  int site_up=neighup(site,mu);
  for(int n=0;n<N;n++) staple+=2.0*conj(zeta(site)[n])*zeta(site_up)[n];
}
inline double link_staple_energy(int site,int mu)
{
  dcomplex staple;
  link_staple(staple,site,mu);
  return -((staple*conj(lambda(site)[mu])).real()-2);
}
inline double link_staple_action(int site,int mu)
{return link_staple_energy(site,mu)/g;}

//initialize the code
void init()
{
  //init the random generators
  rd=new random_device();
  gen=new mt19937_64((*rd)());
  gen->seed(100);
  dis=new uniform_int_distribution<unsigned long long>;
  
  //geometry
  V=1;
  for(int mu=0;mu<ndims;mu++) V*=L;
  cout<<"Volume: "<<V<<endl;
  neigh_data=new int[V*ndims*2];
  
  //loop over sites
  for(int site=0;site<V;site++)
    {
      //get the original coordinates
      coords c;
      coords_of_site(c,site);
      
      //loop over directions
      for(int mu=0;mu<ndims;mu++)
	{
	  //save original
	  int o=c[mu];
	  
	  //backward
	  c[mu]=(o+L-1)%L;
	  neighdw(site,mu)=site_of_coords(c);
	  
	  //forward
	  c[mu]=(o+L+1)%L;
	  neighup(site,mu)=site_of_coords(c);
	  
	  //restore original
	  c[mu]=o;
	}
    }
  
  //Zeta and Lambda
  zeta_data=new dcomplex[N*V];
  lambda_data=new dcomplex[V];

  //set the system to hot state
  init_system_to_hot();
}

//update zeta with metropolis
void metro_update_site(int site)
{
  //copy zeta
  dcomplex ori[N];
  for(int n=0;n<N;n++) ori[n]=zeta(site)[n];
  
  //change of action
  double ori_ac=site_action(site);
  set_ON_to_rnd(zeta(site));
  double fin_ac=site_action(site);
  
  //accept or not?
  double diff_ac=fin_ac-ori_ac;
  double t=exp(-diff_ac);
  double p=get_unif_double(1);
  if(p>t) for(int n=0;n<N;n++) zeta(site)[n]=ori[n];
}

//update lambda with metropolis
void metro_update_link(int site,int mu)
{
  //copy lambdaa
  dcomplex ori=lambda(site)[mu];
  
  //change of action
  double ori_ac=link_action(site,mu);
  set_U1_to_rnd(lambda(site)[mu]);
  double fin_ac=link_action(site,mu);
  
  //accept or not?
  double diff_ac=fin_ac-ori_ac;
  double t=exp(-diff_ac);
  double p=get_unif_double(1);
  if(p>t) lambda(site)[mu]=ori;
}

//update a site using overrelaxion/heatbath or microcanonical
void overheat_micro_update_site(int site,bool over)
{
  //get the staple, its norm and energy
  dcomplex staple[N];
  site_staple(staple,site);
  double staple_norm=get_zeta_norm(staple);
  if(isnan(staple_norm))
    {
      for(int mu=0;mu<ndims;mu++)
	{
	  int site_up=neighup(site,mu);
	  for(int n=0;n<N;n++) cout<<zeta(site_up)[n]<<" "<<conj(lambda(site)[mu])<<endl;
	  int site_dw=neighdw(site,mu);
	  for(int n=0;n<N;n++) cout<<zeta(site_dw)[n]<<" "<<lambda(site_dw)[mu]<<endl;
	}
      for(int n=0;n<N;n++) cout<<staple[n]<<endl;
      CRASH("%lg",staple_norm);
    }
  double staple_energy=get_zeta_scalprod(zeta(site),staple);
  
  //extract theta and compute its cos
  double ctheta_new=cos(get_theta(beta*N*staple_norm,N));
  
  //extract remaining components
  double ctheta_old=staple_energy/staple_norm;
  double a=ctheta_new/staple_norm,b=ctheta_old/staple_norm,c=sqrt((1-sqr(ctheta_new))/(1-sqr(ctheta_old)));
  for(int n=0;n<N;n++) zeta(site)[n]=a*staple[n]-(zeta(site)[n]-b*staple[n])*c;
  
  //reunitarize
  zeta_unitarize(zeta(site));
  //double dev_zeta=check_zeta_unitarity(zeta(site));
  //if(dev_zeta>1.e-10) CRASH("unitarity deviates by %lg for new zeta on site %d",dev_zeta,site);
}

//update a link using overrelaxion/heatbath or microcanonical
void overheat_micro_update_link(int site,int mu,bool over)
{
  static ofstream norm_file("/tmp/norm_file");
  static int nnn=0;
  nnn++;
  
  //get the staple
  dcomplex staple;
  link_staple(staple,site,mu);

  //compute the staple norm and energy
  double staple_norm=sqrt(norm(staple));
  double staple_energy=get_lambda_scalprod(lambda(site)[mu],staple);
  
  //cout<<"site: "<<site<<" staple: "<<staple<<" lambda(site)[mu]: "<<lambda(site)[mu]<<endl;
  
  //extract theta and compute its cos
  double ctheta_new=cos(get_theta_1(beta*N*staple_norm));
  
  //extract remaining components
  double ctheta_old=staple_energy/staple_norm;
  double a=ctheta_new/staple_norm,b=ctheta_old/staple_norm,c=sqrt((1-sqr(ctheta_new))/(1-sqr(ctheta_old)));
  auto term=(lambda(site)[mu]-b*staple);
  lambda(site)[mu]=a*staple-(lambda(site)[mu]-b*staple)*c;

  //reunitarize
  norm_file<<lambda_unitarize(lambda(site)[mu])-1<<endl;
  if(isnan(lambda(site)[mu].real()))
    {
      cout<<"a: "<<a<<" b: "<<b<<" c: "<<c<<" ctheta_old-1: "<<ctheta_old-1<<" chteta_new: "<<ctheta_new<<
	" staple_energy: "<<staple_energy<<" staple_norm: "<<staple_norm<<
	" term: "<<term<<" nnn: "<<nnn<<endl;
      CRASH("site: %d",site);
    }
  //double dev=check_lambda_unitarity(lambda(site)[mu]);
  //if(dev>1.e-15) CRASH("unitarity deviates by %lg for new lambda on site %d mu %d",dev,site,mu);
}

//sweep all the lattice
void metro_sweep()
{
  //loop over sites
  for(int site=0;site<V;site++)
    {
      metro_update_site(site);
      for(int mu=0;mu<ndims;mu++) metro_update_link(site,mu);
    }
}

//sweep all the lattice
void overheat_micro_sweep()
{
  //loop over sites
  for(int site=0;site<V;site++)
    {
      overheat_micro_update_site(site,false);
      for(int mu=0;mu<ndims;mu++) overheat_micro_update_link(site,mu,false);
    }
}

//close the code
void close()
{
  delete[] lambda_data;
  delete[] zeta_data;
  
  delete[] neigh_data;
  
  delete dis;
  delete gen;
  delete rd;
}

int main()
{
  //initialize
  init();
  
  ofstream energy_file("energy");
  
  //sweep with overheat-micro
  for(int i=0;i<20000;i++)
    {
      overheat_micro_sweep();
      //metro_sweep();
      //cout<<energy()/V/ndims<<endl;
    }
  
  //finalize
  close();
  
  return 0;
}
