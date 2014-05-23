#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <cmath>
#include <fstream>
#include <iostream>
#include <string.h>
#include <stdlib.h>

using namespace std;

#define NDIMS 3
#define RAN2_NTAB 32
#define DEBUG 0

#define NTRAJ 20000
#define NSWEEP_PER_TRAJ 20

/*
const double tau=0.17;
const double kappa=0.05;
const int L0=24,L1=24;
*/

/*
const double tau=0.1782;
const double kappa=0.005;
const int L0=8,L1=8;
*/

const double tau=0.176;
const double kappa=0.01;
const int L0=16,L1=16;

const double pot_chim=0;
const int L[3]={L0,L0,L1};

const int mod[9]={2,0,1,
		  0,1,2,
		  1,2,0};

//////////////////////////////////////////////////////////

typedef int coords[NDIMS];
typedef short int z3_t;

//class to generate ranom numbers
class rnd_gen_t
{
public:
  int idum;
  int idum2;
  int iv[RAN2_NTAB];
  int iy;
};

//struct to hold ingredient needed to compute action
struct action_ingr_t
{
  int NS;
  int Ns[3];
  int Nb;
  action_ingr_t();
  double E();
  double U();
};

//struct to hold data
struct N_t
{
  int N;
  int N0;
};

//random generator
rnd_gen_t gen;

//geometry stuff
int V;
int *neigh_data;
z3_t *conf_b,*conf_s;

//probabilities and constants
double C,B,eta,eta_bar,Btilde;
double M[3],Mtilde[3];
double prob_B[9],prob_M[9];

//return the coordinate of a site
int site_of_coords(coords c)
{
  int site=0;
  for(int mu=0;mu<NDIMS;mu++) site=site*L[mu]+c[mu];
  return site;
}

//get coords of site
void coords_of_site(coords c,int site)
{
  for(int mu=NDIMS-1;mu>=0;mu--)
    {
      c[mu]=site%L[mu];
      site/=L[mu];
    }
}

//return the neighbors
inline int &neighdw(int site,int mu)
{return neigh_data[0+2*(mu+NDIMS*site)];}
inline int &neighup(int site,int mu)
{return neigh_data[1+2*(mu+NDIMS*site)];}

//initialize everything
void init_random_gen(int seed);
void init()
{
  init_random_gen(100);
  
  //geometry
  V=1;
  for(int mu=0;mu<NDIMS;mu++) V*=L[mu];
  cout<<"Volume: "<<V<<endl;
  neigh_data=new int[V*NDIMS*2];
  
  //loop over sites
  for(int site=0;site<V;site++)
    {
      //get the original coordinates
      coords c;
      coords_of_site(c,site);
      
      //loop over directions
      for(int mu=0;mu<NDIMS;mu++)
        {
          //save original
          int o=c[mu];
          
          //backward
          c[mu]=(o+L[mu]-1)%L[mu];
          neighdw(site,mu)=site_of_coords(c);
          
          //forward
          c[mu]=(o+L[mu]+1)%L[mu];
          neighup(site,mu)=site_of_coords(c);
          
          //restore original
          c[mu]=o;
        }
    }
  
  //allocate data
  conf_b=new z3_t[V*NDIMS];
  conf_s=new z3_t[V];
  
  //init to cold
  for(int i=0;i<V;i++)
    {
      conf_s[i]=1;
      for(int mu=0;mu<NDIMS;mu++) conf_b[i*NDIMS+mu]=1;
    }
  
  //define the constants used after
  C=(exp(2*tau)+2*exp(-tau))/3;
  B=(exp(2*tau)-exp(-tau))/(exp(2*tau)+2*exp(-tau));
  cout<<"B: "<<B<<endl;

  //lookup tables for prob_B
  for(int btilde=0;btilde<3;btilde++)
    for(int cur_b=0;cur_b<3;cur_b++)
      prob_B[btilde*3+cur_b]=pow(B,abs(btilde-1)-abs(cur_b-1));
  
  //hoppings
  eta=kappa*exp(+pot_chim);
  eta_bar=kappa*exp(-pot_chim);
  for(int s=0;s<3;s++)
    {
      M[s]=(exp(eta+eta_bar)+2*exp(-(eta+eta_bar)/2)*cos((eta-eta_bar)*sqrt(3)/2-(s-1)*2*M_PI/3))/3;
      cout<<"M["<<s<<"]="<<M[s]<<endl;
    }

  //lookup tables for prob_M
  for(int stilde=0;stilde<3;stilde++)
    for(int cur_s=0;cur_s<3;cur_s++)
      prob_M[stilde*3+cur_s]=M[stilde]/M[cur_s];
  
  //extra constants needed for observables
  Btilde=9*tau/(exp(3*tau)+1-exp(-3*tau));
  Mtilde[0]=(eta*M[2]+eta_bar*M[1])/M[0];
  Mtilde[1]=(eta*M[0]+eta_bar*M[2])/M[1];
  Mtilde[2]=(eta*M[1]+eta_bar*M[0])/M[2];
}

//deinitialize everything
void close()
{
  delete[] conf_b;
  delete[] conf_s;
}

//check that at every site the constraint on flux is satisfied
void check_constraint()
{
  //loop on every site
  for(int i=0;i<V;i++)
    {
      //start from value at the site
      int c=(conf_s[i]-1);
      for(int mu=0;mu<NDIMS;mu++)
	{
	  //add flux coming from neighbors
	  c+=(conf_b[i*NDIMS+mu]-1);
	  c-=(conf_b[neighdw(i,mu)*NDIMS+mu]-1);
	}
      //check
      if(((c>0) && ((c)%3))||((c<0) && ((-c)%3))) cout<<"failed constraint at: "<<i<<" "<<c<<endl;
    }
}

//return a double between [0,max)
double get_unif_double(double max,bool incl=false)
{
  double res;
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
  
  return res;
}

//generate
void init_random_gen(int seed)
{
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
}

//compute the ingredient needed to obtain action
action_ingr_t::action_ingr_t()
{
  //reset
  NS=Ns[0]=Ns[1]=Ns[2]=Nb=0;
  
  //count
  for(int i=0;i<V;i++)
    {
      //occupation per type
      Ns[conf_s[i]]++;
      //sites non zero
      NS+=(conf_s[i]!=1);
      //links non zero
      for(int mu=0;mu<NDIMS;mu++) Nb+=(conf_b[i*NDIMS+mu]!=1);
    }
}

//return the energy
double action_ingr_t::E()
{return Nb*Btilde/(V*tau);}

//return the total action
double action_ingr_t::U()
{return (2*NDIMS*V*tau*B+Nb*Btilde+Ns[0]*Mtilde[0]+Ns[1]*Mtilde[1]+Ns[2]*Mtilde[2])/((double)V);}

//////////////////////////////////////////////////////////////////////////////

struct worm_t
{
  int x,x0,worm_sign,worm_length;
  bool worm_complete;  
  void init();
  bool insert_monomer_first_tail();
  void insert_monomer_second_tail();
  void insert_dimer_pair();
  void move(int nu);
  void update();
};

//insert the first tail of a monomer pair
bool worm_t::insert_monomer_first_tail()
{
  int cur_s=conf_s[x];
  int stilde=mod[cur_s*3+worm_sign];
  double extr=get_unif_double(1);
  double prob=prob_M[stilde*3+cur_s];
  if(DEBUG)
    cout<<" trying to insert monomer "<<stilde<<" (first tail) at "<<x<<", (currently: "<<conf_s[x]
	<<") sign: "<<worm_sign<<" prob: "<<prob<<" extr: "<<extr<<endl;
  bool succ=(extr<=prob);
  
  //if accepted mark it down and switch state
  if(succ)
    {
      conf_s[x]=stilde;
      if(DEBUG) cout<<"  inserted a monomer (first tail) at "<<x<<": "<<stilde<<endl;
    }
  
  return succ;
}

//insert the secont tail of a monomer tail
void worm_t::insert_monomer_second_tail()
{
  bool succ=false;
  
  while(!succ)
    {
      int x1=(int)get_unif_double(V);
      int cur_s=conf_s[x1];
      int stilde=mod[cur_s*3+(2-worm_sign)];
      double extr=get_unif_double(1);
      double prob=prob_M[stilde*3+cur_s];
      if(DEBUG) cout<<" trying to insert a monomer (second tail) at "<<x1<<", (currrently: "<<cur_s
		    <<") sign: "<<worm_sign<<" prob: "<<prob<<" extr: "<<extr<<endl;
      succ=(extr<=prob);
      
      //if accepted mark it and switch to normal state
      if(succ)
	{
	  conf_s[x1]=stilde;
	  x=x1;
	  if(DEBUG) cout<<"  inserted a monomer (second tail) at "<<x<<": "<<stilde<<endl;
	}
    }
}

//insert a pair of fimer
void worm_t::insert_dimer_pair()
{
  if(insert_monomer_first_tail())
    insert_monomer_second_tail();
}

//initialize worm
void worm_t::init()
{
  //extract the starting site and sign
  x0=(int)get_unif_double(V);
  worm_sign=(int)get_unif_double(2)*2;

  if(DEBUG) cout<<"starting at "<<x0<<" sign: "<<worm_sign<<endl;
  
  //initialize worm lengthm state and position
  worm_length=0;
  x=x0;
  worm_complete=false;  
}

//perform an ordinary move
void worm_t::move(int nu)
{
  int rho=abs(nu)-1;
  int x1=((nu>0)?x:neighdw(x,rho));
  int cur_b=conf_b[x1*NDIMS+rho];
  int btilde=mod[cur_b*3+((nu>0)?worm_sign:(2-worm_sign))];
  if(DEBUG) cout<<" trying to move from "<<x<<" (currently: "<<cur_b
		<<") in direction: "<<rho<<(nu>0?'+':'-')<<" to change it to "<<btilde<<endl;
  if(get_unif_double(1)<=prob_B[btilde*3+cur_b])
    {
      conf_b[x1*NDIMS+rho]=btilde;
      x=((nu>0)?neighup(x,rho):neighdw(x,rho));
      if(DEBUG) cout<<"  moved in direction: "<<rho<<endl;
      worm_length++;
    }
}

//make a whole step
void worm_t::update()
{
  //initialize the worm
  init();
  
  //loop until the worm is closed
  while(!worm_complete)
    {
      int nu=(int)get_unif_double(2*NDIMS+1)-NDIMS;
      
      //insert the monomer first tail
      if(nu==0) insert_dimer_pair();
      else move(nu);
      
      //check closing
      if(x==x0 && worm_length!=0) worm_complete=true;
    }

  if(DEBUG) cout<<"Worm length: "<<worm_length<<", sign: "<<worm_sign<<endl;
  if(DEBUG) check_constraint();
}

int main(int narg,char **arg)
{
  init();  
  
  ofstream energy_file("/tmp/worm/energy.xmg");
  ofstream action_file("/tmp/worm/action.xmg");
  for(int iconf=0;iconf<NTRAJ;iconf++)
    {
      worm_t worm;
      for(int i=0;i<NSWEEP_PER_TRAJ;i++) worm.update();
      
      action_ingr_t act;
      energy_file<<act.E()<<endl;
      action_file<<act.U()<<endl;
    }
  
  close();
  
  return 0;
}
