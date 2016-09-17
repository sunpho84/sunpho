#include "utils.hpp"
#include "vectors.hpp"

#include <algorithm>
#include <functional>
#include <array>
#include <iostream>
#include <vector>

using namespace std;

void test_dvect()
{
  printf("Double vector number of entries: %lu\n",dvect_nel);
  
  dvect a={0,1,2,3};
  dvect b={2,3,4,5};
  dregv ra=dvect_load(a);
  dregv rb=dvect_load(b);
  dregv rc=dregv_summ(ra,rb);
  
  dvect c;
  dregv_store(c,rc);
  
  dvect_printf(c);
  
  printf("-----\n");
  dregv_store(c,dregv_shift_dw(ra,1));
  dvect_printf(a);
  dvect_printf(c);
  
  printf("-----\n");
  dregv_store(c,dregv_shift_dw(ra,2));
  dvect_printf(a);
  dvect_printf(c);
  
  printf("-----\n");
  dregv_store(c,dregv_shift_dw(ra,3));
  dvect_printf(a);
  dvect_printf(c);
  
  printf("-----\n");
  dregv_store(c,dregv_shift_up(ra,3));
  dvect_printf(a);
  dvect_printf(c);
}

void test_cdvect()
{
  //dregv d;
  //dregv_zero(d);
  
  cdvect a={{0,2,3,4},{1,0,0,0}};
  cdregv ra;
  cdvect_load(ra,a);
  printf("a=");cdvect_printf(a);
  cdvect b={{0,0,0,0},{1,1,1,1}};
  cdregv rb;
  cdvect_load(rb,b);
  printf("b=");cdvect_printf(a);
  
  cdregv rc;
  cdregv_prod(rc,ra,rb);
  cdvect c;
  cdregv_store(c,rc);
  printf("aXb=");cdvect_printf(c);
  
  ///////////////////////////////
  cdvect d={{1,1,1,1},{1,1,1,1}};
  cdregv rd;
  cdvect_load(rd,d);
  printf("d=");cdvect_printf(d);
  cdregv_summ_the_prod(rd,ra,rb);
  cdvect e;
  cdregv_store(e,rd);
  printf("d+aXb=");cdvect_printf(e);
  
  cdvect_load(rd,d);
  cdregv_subt_the_prod(rd,ra,rb);
  cdregv_store(e,rd);
  printf("d-aXb=");cdvect_printf(e);
}

const int ndim=2;
const int nneigh=2*ndim;
typedef array<int,ndim> coord_t;
typedef array<int,nneigh> neigh_t;
typedef pair<int,int> vind_t;

#define FOR_DIR(mu) for(int mu=0;mu<ndim;mu++)

//! return the coordinate of site inside the given box
coord_t coord_of_site(coord_t box_size,int site)
{
  coord_t c;
  for(int mu=ndim-1;mu>=0;mu--)
    {
      c[mu]=site%box_size[mu];
      site/=box_size[mu];
    }
  return c;
}

//! return the site of a given coordinates in the given box
int site_of_coord(coord_t box_size,coord_t c)
{
  int site=0;
  FOR_DIR(mu) site=site*box_size[mu]+c[mu];
  return site;
}

int main(int narg,char **arg)
{
  //! common size
  int L=32;
  //! size per dir
  coord_t lx_size_per_dir={L,L};
  //! total size
  int vol=L*L;
  
  //! number of vranks
  const int dvgrid_nvranks=dvect_nel;
  //! nvranks per dir
  coord_t dvgrid_nvranks_per_dir={2,2};
  
  //! volume per vrank
  int dvrank_vol=vol/dvgrid_nvranks;
  //! size per dir of a vectorized rank
  coord_t dvrank_size_per_dir;
  transform(lx_size_per_dir.begin(),lx_size_per_dir.end(),dvgrid_nvranks_per_dir.begin(),dvrank_size_per_dir.begin(),divides<double>());
  
  //! coordinate in the lx grid
  vector<coord_t> lx_coord_of_lx_site(vol);
  for(int lx=0;lx<vol;lx++) lx_coord_of_lx_site[lx]=coord_of_site(lx_size_per_dir,lx);
  
  //! index in the vector ordering of an lx site and vice-versa
  vector<vind_t> dvgrid_site_of_lx_site(vol);
  vector<array<int,dvgrid_nvranks>> lx_site_of_dvgrid_site(dvrank_vol);
  for(int lx=0;lx<vol;lx++)
    {
      coord_t dvrank_coord,dvgrid_site_coord;
      FOR_DIR(mu)
	{
	  dvrank_coord[mu]=lx_coord_of_lx_site[lx][mu]/dvrank_size_per_dir[mu];
	  dvgrid_site_coord[mu]=lx_coord_of_lx_site[lx][mu]%dvrank_size_per_dir[mu];
	}
      int dvgrid_site=site_of_coord(dvrank_size_per_dir,dvgrid_site_coord);
      int dvgrid_rank=site_of_coord(dvgrid_nvranks_per_dir,dvrank_coord);
      dvgrid_site_of_lx_site[lx]=make_pair(dvgrid_site,dvgrid_rank);
      lx_site_of_dvgrid_site[dvgrid_site][dvgrid_rank]=lx;
      
      //cout<<lx<<" "<<dvgrid_site_of_lx_site[lx].first<<" "<<dvgrid_site_of_lx_site[lx].second<<"  "<<lx_site_of_dvgrid_site[dvgrid_site][dvgrid_rank]<<endl;
    }
  
  //! border size in the vector ordering
  int dvrank_bord_vol=0;
  FOR_DIR(mu) dvrank_bord_vol+=dvrank_size_per_dir[mu]*2;
  
  //! neighbors in all directions, including border, in the dvgrid
  vector<neigh_t> dvgrid_site_neigh(dvrank_vol+dvrank_bord_vol);
  for(auto &it : dvgrid_site_neigh) it.fill(-1);
  //search neighbors
  int iv_bord=dvrank_vol;
  for(int orie=0;orie<2;orie++)
    FOR_DIR(mu)
      for(int dvgrid_site=0;dvgrid_site<dvrank_vol;dvgrid_site++)
	{
	  //! site in the lx grid
	  int lx=lx_site_of_dvgrid_site[dvgrid_site][0];
	  //! lx coordinate
	  coord_t &c=lx_coord_of_lx_site[lx];
	  //! full direction
	  int dir=orie*ndim+mu;
	  //! shift backward or forward
	  int sign_orie[2]={-1,+1};
	  coord_t n=c;
	  //shift
	  n[mu]+=sign_orie[orie];
	  
	  //check if it outside
	  if(n[mu]<0||n[mu]>=dvrank_size_per_dir[mu])
	    {
	      dvgrid_site_neigh[dvgrid_site][dir]=iv_bord;
	      dvgrid_site_neigh[iv_bord][dir]=dvgrid_site;
	      iv_bord++;
	    }
	  else dvgrid_site_neigh[dvgrid_site][dir]=site_of_coord(lx_size_per_dir,n);
	}
  
  for(int dvgrid_site=0;dvgrid_site<dvrank_vol+dvrank_bord_vol;dvgrid_site++)
    {
      cout<<dvgrid_site;
      for(int dir=0;dir<2*ndim;dir++) cout<<"\t"<<dvgrid_site_neigh[dvgrid_site][dir];
      cout<<endl;
    }
  
  //test_dvect();
  //test_cdvect();
  
  return 0;
}
