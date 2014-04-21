#define HOT 1
#define COLD 0
#define TINY 1.e-13
#define HALFWAY_PRECISION 1.e-8

//#define USE_GOOD_GENERATOR
#define USE_SMART_EXTRACTION

//types and constants
#define ndims 2
typedef complex<double> dcomplex;
typedef int coords[ndims];

//random number generators
#ifdef USE_GOOD_GENERATOR
 random_device *rd;
 mt19937_64 *gen;
 uniform_real_distribution<double> *dis;
#else
 #define RAN2_NTAB 32
 //The structure for the random generator
 struct rnd_gen
 {
   int idum;
   int idum2;
   int iv[RAN2_NTAB];
   int iy;
 };
 rnd_gen gen;
#endif

//parameters
int N=21;
double beta=0.7;
int L=72;
double g=1/(N*beta);

//geometry
int V;
int *neigh_data;

//system
int init_time;
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
double get_lambda_scalprod(dcomplex a,dcomplex b)
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
#ifdef USE_GOOD_GENERATOR
  do res=max*(*dis)(*gen)/dis->b();
  while((!incl)&&(res==max));
#else
  const int im1=2147483563,im2=2147483399,imm1=im1-1,ia1=40014,ia2=40692;
  const int iq1=53668,iq2=52774,ir1=12211,ir2=3791,ndiv=1+imm1/RAN2_NTAB;
  const double am=1.0/im1,eps=1.2e-7,rnmx=1-eps;
  int j,k;
    
  k=gen.idum/iq1;
  gen.idum=ia1*(gen.idum-k*iq1)-k*ir1;
  if(gen.idum<0) gen.idum+=im1;
    
  k=gen.idum2/iq2;
  gen.idum2=ia2*(gen.idum2-k*iq2)-k*ir2;
  if(gen.idum2<0) gen.idum2+=im2;
    
  j=gen.iy/ndiv;
  gen.iy=gen.iv[j]-gen.idum2;
  gen.iv[j]=gen.idum;
  if(gen.iy<0) gen.iy+=imm1;
    
  res=max*std::min(am*gen.iy,rnmx);
#endif
  
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

#ifdef USE_SMART_EXTRACTION
//obtain theta
double get_theta(double a,int k)
{
  //compute parameters
  double zita=(k-1)/a;
  double theta0=acos(sqrt(1+zita*zita)-zita);
  double ctheta0=cos(theta0);
  double c=sqrt(2*(k-1)*(1-zita*ctheta0)/sqr(sin(theta0)));
  double ptheta0=fun_ptheta(theta0,a,k);
  
  //extract theta
  double theta;
  double eta=0.99;
  bool acc;
  int no=0;
  do
    {
      double chi=get_unif_double(1);
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
double get_theta_1(double a)
{
  double eps=0.001;
  double as=0.798953686083986;
  double dap=max((double)0.0,a-as);
  double del=0.35*dap+1.03*sqrt(dap);
  double alp=min(sqrt(a*(2-eps)),max(sqrt(eps*a),del));
  double bet=max(alp*alp/a,(double)((cosh(M_PI*alp)-1)/(exp(2*a)-1)))-1;
  double bt1=sqrt((1+bet)/(1-bet));
  
  double h; //result
  bool acc; //accepted or not
  int no=0;
  do
    {
      double r=get_unif_double(1);
      double h1=bt1*tan((2*r-1)*atan(tanh(M_PI*alp/2)/bt1));
      h=log((1+h1)/(1-h1))/alp;
      
      //decide if accept or reject
      double g=exp(-a*(1-cos(h)))*(cosh(alp*h)+bet)/(1+bet);
      if(g>1+TINY) CRASH("%lg",g-1);
      double p=get_unif_double(1);
      acc=(p<g);
      if(!acc) no++;
      if(no>1000) CRASH("%d",no);
    }
  while(!acc);
  
  return h;
}
#else
double get_theta(double a,int k)
{
  double th,no=exp(-fabs(a));
  
  int non=0;
  bool acc;
  do
    {
      th=get_unif_double(2*M_PI);
      double pacc=no*pow(sin(th),2*(k-1))*exp(a*cos(th));
      if(pacc>1) CRASH("a%lg",pacc);
      double ext=get_unif_double(1);
      acc=(ext<pacc);
      if(!acc) non++;
      if(non>10000) CRASH("non: %d, a: %lg",non,a);
    }
  while(!acc);
  
  return th;
}

double get_theta_1(const double a)
{
  double th,no=exp(-fabs(a));
  
  int non=0;
  bool acc;
  do
    {
      th=get_unif_double(2*M_PI);
      double pacc=no*exp(a*cos(th));
      if(pacc>1+TINY) CRASH("a%lg",pacc-1);
      double ext=get_unif_double(1);
      acc=(ext<pacc);
      if(!acc) non++;
      if(non>100) CRASH("non: %d",non);
    }
  while(!acc);
  
  return th;
}
#endif

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

//initialize the system to hot
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

//initialize to cold
void init_system_to_cold()
{
  for(int site=0;site<V;site++)
    {
      //fill the lambda
      for(int mu=0;mu<ndims;mu++) lambda(site)[mu]=1;
      
      //fill the Zeta
      for(int n=0;n<N;n++) zeta(site)[n]=(n==0);
    }
}

//switch
void init_system_to(int cond)
{
  if(cond==HOT) init_system_to_hot();
  else          init_system_to_cold();
}

//control that each element is unitary
void check_system_unitarity(double res=TINY)
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
      for(int n=0;n<N;n++) staple[n]+=(double)2.0*zeta(site_up)[n]*conj(lambda(site)[mu]);
      int site_dw=neighdw(site,mu);
      for(int n=0;n<N;n++) staple[n]+=(double)2.0*zeta(site_dw)[n]*lambda(site_dw)[mu];
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
  for(int n=0;n<N;n++) staple+=(double)2.0*conj(zeta(site)[n])*zeta(site_up)[n];
}
inline double link_staple_energy(int site,int mu)
{
  dcomplex staple;
  link_staple(staple,site,mu);
  return -((staple*conj(lambda(site)[mu])).real()-2);
}
inline double link_staple_action(int site,int mu)
{return link_staple_energy(site,mu)/g;}

//return inner product of zeta
void get_P(dcomplex *P,dcomplex *z)
{
  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      P[i*N+j]=conj(z[i])*z[j];
}

//return the angle of scalar product between two zetas
double arg_an(dcomplex *a,dcomplex *b)
{
  dcomplex res=0;
  for(int n=0;n<N;n++) res+=conj(a[n])*b[n];
  //cout<<log(res)<<" "<<log(conj(res))<<endl;
  return arg(res);
}
dcomplex sc(dcomplex *a,dcomplex *b)
{
  dcomplex res=0;
  for(int n=0;n<N;n++) res+=conj(a[n])*b[n];
  return res;
}

//return the geometric definition of topology
double geometric_topology_simplified()
{
  double topo=0;
  for(int n=0;n<V;n++)
    {
      int mu=0,nu=1;
      int nmu=neighup(n,mu);
      int nnu=neighup(n,nu);
      int nmu_nu=neighup(nmu,nu);
      
      topo+=arg(sc(zeta(nmu_nu),zeta(n))*sc(zeta(nmu),zeta(nmu_nu))*sc(zeta(n),zeta(nmu)))+
	arg(sc(zeta(nnu),zeta(n))*sc(zeta(nmu_nu),zeta(nnu))*sc(zeta(n),zeta(nmu_nu)));
    }
  
  return topo/(2*M_PI);
}
double geometric_topology()
{
  int mu=0,nu=1;
  double topo=0;
  for(int n=0;n<V;n++)
    {
      dcomplex P1[N*N],P2[N*N],P3[N*N];
      get_P(P1,zeta(n));
      get_P(P3,zeta(neighup(neighup(n,mu),nu)));

      dcomplex c;
      
      c=0;
      get_P(P2,zeta(neighup(n,mu)));
      for(int i=0;i<N;i++)
	for(int j=0;j<N;j++)
	  for(int k=0;k<N;k++)
	    c+=P3[i*N+j]*P2[j*N+k]*P1[k*N+i];
      topo+=log(c).imag();
      
      c=0;
      get_P(P2,zeta(neighup(n,nu)));
      for(int i=0;i<N;i++)
	for(int j=0;j<N;j++)
	  for(int k=0;k<N;k++)
	    c+=P2[i*N+j]*P3[j*N+k]*P1[k*N+i];
      topo+=log(c).imag();
    }
  
  return topo/(2*M_PI);
}

double topology()
{
  int mu=0,nu=1;
  double topo=0;
  for(int n=0;n<V;n++)
    topo+=(lambda(n)[mu]*lambda(neighup(n,mu))[nu]*conj(lambda(neighup(n,nu))[mu]*lambda(n)[nu])).imag();
  
  return topo;
}

//initialize the code
void init(int cond,int seed)
{
  init_time=time(0);
  
#ifdef USE_GOOD_GENERATOR
  //init the random generators
  rd=new random_device();
  gen=new mt19937_64((*rd)());
  gen->seed(seed);
  dis=new uniform_real_distribution<double>;
#else
  const int im1=2147483563,ia1=40014;
  const int iq1=53668,ir1=12211;
  int j,k;
  
  //initialization
  gen.idum=seed;
  gen.idum=std::max(gen.idum+1,1);
  gen.idum2=gen.idum;
  for(j=RAN2_NTAB+7;j>=0;j--)
    {
      k=gen.idum/iq1;
      gen.idum=ia1*(gen.idum-k*iq1)-k*ir1;
      if(gen.idum<0) gen.idum+=im1;
      if(j<RAN2_NTAB) gen.iv[j]=gen.idum;
    }
  gen.iy=gen.iv[0];
#endif
    
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
  lambda_data=new dcomplex[V*ndims];

  //set the system to hot state
  init_system_to(cond);
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

//update a site using microcanonical
void micro_update_site(int site)
{
  //get the staple, its norm and energy
  dcomplex staple[N];
  site_staple(staple,site);
  double staple_norm=get_zeta_norm(staple);
  double staple_energy=get_zeta_scalprod(zeta(site),staple);
  
  //extract site
  for(int n=0;n<N;n++) zeta(site)[n]=2*staple_energy/sqr(staple_norm)*staple[n]-zeta(site)[n];

  //reunitarize
  zeta_unitarize(zeta(site));
}

//update a site using overrelaxion/heatbath
void overheat_update_site(int site)
{
  //get the staple, its norm and energy
  dcomplex staple[N];
  site_staple(staple,site);
  double staple_norm=get_zeta_norm(staple);
  double staple_energy=get_zeta_scalprod(zeta(site),staple);
  
  //compute theta in the simple way
  double ctheta_old=staple_energy/staple_norm;
  double theta_old=acos(ctheta_old);
  
  //if theta is too small we switch to alternative method
  if(fabs(theta_old)<1.e-4) 
    {
      dcomplex diff[N];
      for(int n=0;n<N;n++) diff[n]=zeta(site)[n]-staple[n]/staple_norm;
      theta_old=asin(get_zeta_scalprod(diff,diff));
    }
  
  //it theta is too small the algorithm is undefined
  if(fabs(theta_old)>5.e-8)
    {
      //extract theta and compute its cos
      double theta_new=get_theta(beta*N*staple_norm,N);
      double ctheta_new=cos(theta_new);
      
      //extract remaining components
      double a=ctheta_new/staple_norm,b=ctheta_old/staple_norm,c=sin(theta_new)/sin(theta_old);
      for(int n=0;n<N;n++) zeta(site)[n]=a*staple[n]-(zeta(site)[n]-b*staple[n])*c;
      
      //reunitarize
      zeta_unitarize(zeta(site));
    }
  else cout<<"skipping site "<<site<<": "<<theta_old<<endl;
}

//update a link using microcanonical
void micro_update_link(int site,int mu)
{
  //get the staple
  dcomplex staple;
  link_staple(staple,site,mu);

  //compute the staple norm and energy
  double staple_norm=sqrt(norm(staple));
  double staple_energy=get_lambda_scalprod(lambda(site)[mu],staple);
  
  //extract link
  lambda(site)[mu]=2*staple_energy/sqr(staple_norm)*staple-lambda(site)[mu];

  //reunitarize
  lambda_unitarize(lambda(site)[mu]);
}

//update a link using overrelaxion/heatbath
void overheat_update_link(int site,int mu)
{
  //get the staple
  dcomplex staple;
  link_staple(staple,site,mu);

  //compute the staple norm and energy
  double staple_norm=sqrt(norm(staple));
  double staple_energy=get_lambda_scalprod(lambda(site)[mu],staple);
  
  //compute theta in the simple way
  double ctheta_old=staple_energy/staple_norm;
  double theta_old=acos(ctheta_old);
  
  //if theta is too small we switch to alternative method
  if(fabs(theta_old)<1.e-4) 
    {
      dcomplex diff=lambda(site)[mu]-staple/staple_norm;
      theta_old=asin(get_lambda_scalprod(diff,diff));
    }
  
  //it theta is too small the algorithm is undefined
  if(fabs(theta_old)>5.e-8)
    {
      //extract theta and compute its cos
      double theta_new=get_theta_1(beta*N*staple_norm);
      double ctheta_new=cos(theta_new);
      
      //extract remaining components  
      double a=ctheta_new/staple_norm,b=ctheta_old/staple_norm,c=sin(theta_new)/sin(theta_old);
      lambda(site)[mu]=a*staple-(lambda(site)[mu]-b*staple)*c;
      
      //reunitarize
      lambda_unitarize(lambda(site)[mu]);
    }
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

//sweep all the lattice with microcanonical
void micro_sweep()
{
  //loop over sites
  for(int site=0;site<V;site++)
    {
      micro_update_site(site);
      for(int mu=0;mu<ndims;mu++) micro_update_link(site,mu);
    }
}

//sweep all the lattice with overrelaxation/heatbath
void overheat_sweep()
{
  //loop over sites
  for(int site=0;site<V;site++)
    {
      overheat_update_site(site);
      for(int mu=0;mu<ndims;mu++) overheat_update_link(site,mu);
    }
}

//perform a hybid monte carlo update
void hmc_update()
{
  //allocate momenta
  dcomplex *pi=new dcomplex[V*N];
  double *omega=new double[V*ndims];
  
  //draw momenta
  for(int site=0;site<V;site++)
    {
      
    }
  
  delete[] pi;
  delete[] omega;
}

//close the code
void close()
{
  delete[] lambda_data;
  delete[] zeta_data;
  
  delete[] neigh_data;
  
#ifdef USE_GOOD_GENERATOR
  delete dis;
  delete gen;
  delete rd;
#endif
}

int main()
{
  //initialize
  init(HOT,100);

  ofstream energy_file("energy");
  energy_file.precision(16);
  ofstream topology_file("topology");
  
  //sweep with overheat-micro
  int nsweep=1000000;
  for(int isweep=1;isweep<=nsweep;isweep++)
    {
      //metro_sweep();
      for(int imicro=0;imicro<3;imicro++) micro_sweep();
      overheat_sweep();
      
      double topo_sim=geometric_topology_simplified();
      double topo_num=topology();
      
      energy_file<<energy()/V/ndims<<endl;
      topology_file<<topo_sim<<" "<<topo_num<<endl;
      
      //write time progress
      if(isweep%(nsweep/100)==0) cout<<isweep*100/nsweep<<"%, "<<time(0)-init_time<<" s"<<endl;
    }
  
  //finalize
  close();
  
  return 0;
}
