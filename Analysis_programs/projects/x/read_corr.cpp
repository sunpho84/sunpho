#include "include.h"

int L=24;
int njacks=15;
  
int nind_mom;
int *deg,*ind_of_lx;
double *m2;

typedef int coords[4];
coords size={L,L,L,L};
int V=L*L*L*L;
double phase[4]={1,0.01,0.02,0.03};

//Return the index of site of coord x in a box of sides s
int lx_of_coord(coords x)
{
  int ilx=0;
  for(int mu=0;mu<4;mu++) ilx=ilx*size[mu]+x[mu];
  return ilx;
}
void coord_of_lx(coords x,int ilx)
{
  for(int mu=3;mu>=0;mu--)
    {
      x[mu]=ilx%size[mu];
      ilx/=size[mu];
    }
}

//find symmetric coord
void reduced_coord(coords r,coords c)
{
  for(int mu=0;mu<4;mu++)
    if(c[mu]<=size[mu]/2) r[mu]=c[mu];
    else r[mu]=size[mu]-c[mu];
}

//check if it is independent
bool is_ind(int imom)
{
  //get coord
  coords c,r;
  coord_of_lx(c,imom);
  reduced_coord(r,c);
  
  bool is=true;
  for(int mu=0;mu<4;mu++)
    {
      is&=(c[mu]<=size[mu]/2);
	  
      //check if reduced coordinate is smaller or equal to following (excluded time)
      if(mu>0)
	for(int nu=mu+1;nu<4;nu++)
	  is&=(r[mu]<=r[nu]);
    }
  return is;
}

//count the number of independent momenta
void count_nind_mom()
{
  //mark independent if it is
  for(int imom=0;imom<V;imom++)
    if(is_ind(imom))
      nind_mom++;
}

//symmetrize a certain direction
void symmetrize(jvec &corr,int mu)
{
  for(int imom=0;imom<V;imom++)
    {
      //take coords
      coords c;
      coord_of_lx(c,imom);
      
      //take symm
      c[mu]=(size[mu]-c[mu])%size[mu];
      int jmom=lx_of_coord(c);
      jack t=(corr[imom]+corr[jmom])/2;
      
      //overwrite
      corr[imom]=corr[jmom]=t;
    }
}
void symmetrize(jvec &corr)
{for(int mu=0;mu<4;mu++) symmetrize(corr,mu);}

//permute
void permute(jvec &corr)
{
  int mu[6]={1,2,3,3,1,2};
  int nu[6]={2,3,1,2,3,1};
  int rh[6]={3,1,2,1,2,3};
  
  for(int imom=0;imom<V;imom++)
    {
      //take coords
      coords c;
      coord_of_lx(c,imom);
      
      int jmom[6];
      jack t(njacks);
      t=0;
      for(int iperm=0;iperm<6;iperm++)
	{
	  //permute coords
	  coords d;
	  d[0]=c[0];
	  d[1]=c[mu[iperm]];
	  d[2]=c[nu[iperm]];
	  d[3]=c[rh[iperm]];
	  
	  //take note and average
	  jmom[iperm]=lx_of_coord(d);
	  t+=corr[jmom[iperm]];
	}
      
      //overwrite
      t/=6;
      for(int iperm=0;iperm<6;iperm++) corr[jmom[iperm]]=t;
    }
}

//set m2
void compute_momenta()
{
  for(int imom=0;imom<V;imom++)
    {
      //compute reduced coord
      coords c,r;
      coord_of_lx(c,imom);
      reduced_coord(r,c);
      
      //compute distances
      m2[imom]=0;
      for(int mu=0;mu<4;mu++)
	{
	  double m=(2*r[mu]+phase[mu])*M_PI/size[mu];
	  m2[imom]+=m*m;
	}
    }
}

bool a(int imom)
{
  //compute reduced coord
  coords c,r;
  coord_of_lx(c,imom);
  reduced_coord(r,c);
  
  return true;//(c[0]+c[1]+c[2]+c[3]==1);
}

int main()
{
  //count the momenta
  count_nind_mom();
  cout<<"Number of independent momenta: "<<nind_mom<<endl;
  debug_load=0;
  
  //load data
  jvec corr(V,njacks);
  corr.load("XspaceCorrClust",0);
  
  //compute distances
  m2=new double[V];
  compute_momenta();
  
  //open output
  ofstream out("plot.xmg");
  out<<"@type xydy"<<endl;
  out<<"@yaxes scale Logarithmic"<<endl;
  
  //write all
  out<<"@s0 line type 0"<<endl;
  for(int imom=0;imom<V;imom++) if(a(imom)) out<<sqrt(m2[imom])<<" "<<corr[imom]<<endl;
  out<<"&"<<endl;
  //symmetrize
  out<<"@s1 line type 0"<<endl;
  symmetrize(corr);
  for(int imom=0;imom<V;imom++) if(a(imom)) out<<sqrt(m2[imom])+0.05*1<<" "<<corr[imom]<<endl; 
  out<<"&"<<endl;
  //permute
  out<<"@s2 line type 0"<<endl;
  permute(corr);
  for(int imom=0;imom<V;imom++) if(a(imom)) out<<sqrt(m2[imom])+0.05*2<<" "<<corr[imom]<<endl;
  out<<"&"<<endl;
  //filter
  out<<"@s3 line type 0"<<endl;
  for(int imom=0;imom<V;imom++) if(a(imom)) if(is_ind(imom)) out<<sqrt(m2[imom])+0.05*3<<" "<<corr[imom]<<endl;
  out<<"&"<<endl;

  delete[] m2;
  
  return 0;
}
