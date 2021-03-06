#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>
#include <string.h>

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace bissa
{
  //Return the index of site of coord x in the border mu
  int bordlx_of_coord(coords x,int mu)
  {
    int ilx=0;  
    for(int nu=0;nu<NDIM;nu++)
      if(nu!=mu)
	ilx=ilx*loc_size[nu]+x[nu];
    
    return ilx;
  }
  
  //Return the index of site of coord x in a box of sides s
  int lx_of_coord(coords x,coords s)
  {
    int ilx=0;
    
    for(int mu=0;mu<NDIM;mu++)
      ilx=ilx*s[mu]+x[mu];
    
    return ilx;
  }
  void coord_of_lx(coords x,int ilx,coords s)
  {
    for(int mu=3;mu>=0;mu--)
      {
	x[mu]=ilx%s[mu];
	ilx/=s[mu];
      }
  }
  
  //wrappers
  int loclx_of_coord(coords x)
  {return lx_of_coord(x,loc_size);}
  
  //wrappers
  int glblx_of_coord(coords x)
  {return lx_of_coord(x,glb_size);}
  
  //combine two points
  int glblx_of_comb(int b,int wb,int c,int wc)
  {
    coords co;
    for(int mu=0;mu<NDIM;mu++)
      {
	co[mu]=glb_coord_of_loclx[b][mu]*wb+glb_coord_of_loclx[c][mu]*wc;
	while(co[mu]<0) co[mu]+=glb_size[mu];
	co[mu]%=glb_size[mu];
      }
    
    return glblx_of_coord(co);
  }
  
  void glb_coord_of_glblx(coords x,int gx)
  {
    for(int mu=3;mu>=0;mu--)
      {
	int next=gx/glb_size[mu];
	x[mu]=gx-next*glb_size[mu];
	gx=next;
      }
  }
  
  int glblx_of_diff(int b,int c)
  {return glblx_of_comb(b,+1,c,-1);}
  
  int glblx_of_summ(int b,int c)
  {return glblx_of_comb(b,+1,c,+1);}
  
  int glblx_opp(int b)
  {return glblx_of_diff(0,b);}
  
  //Return the coordinate of the rank containing the global coord
  void rank_coord_of_site_of_coord(coords rank_coord,coords glb_coord)
  {for(int mu=0;mu<NDIM;mu++) rank_coord[mu]=glb_coord[mu]/loc_size[mu];}
  
  //Return the rank of passed coord
  int rank_of_coord(coords x)
  {return lx_of_coord(x,nrank_dir);}
  void coord_of_rank(coords c,int x)
  {coord_of_lx(c,x,nrank_dir);}
  
  //Return the rank containing the global coordinates
  int rank_hosting_site_of_coord(coords x)
  {
    coords p;
    rank_coord_of_site_of_coord(p,x);
    
    return rank_of_coord(p);
  }
  //Return the rank containing the glblx passed
  int rank_hosting_glblx(int gx)
  {
    coords c;
    glb_coord_of_glblx(c,gx);
    return rank_hosting_site_of_coord(c);
  }
  
  //Return the local site and rank containing the global coordinates
  void get_loclx_and_rank_of_coord(int *ivol,int *rank,coords g)
  {
    coords l,p;
    for(int mu=0;mu<NDIM;mu++)
      {
	p[mu]=g[mu]/loc_size[mu];
	l[mu]=g[mu]-p[mu]*loc_size[mu];
      }
    
    (*rank)=rank_of_coord(p);
    (*ivol)=loclx_of_coord(l);
  }
  
  //Return the global index of site addressed by rank and loclx
  int get_glblx_of_rank_and_loclx(int irank,int loclx)
  {
    coords p;
    coord_of_rank(p,irank);
    
    int iglblx=0;
    for(int mu=0;mu<NDIM;mu++)
      iglblx=iglblx*glb_size[mu]+loc_coord_of_loclx[loclx][mu];
    
    return iglblx;
  }
  
  //return the index of the site of passed "pseudolocal" coordinate
  //if the coordinates are local, return the index according to the function loclx_of_coord
  //if exactly one of the coordinate is just out return its index according to bordlx_of_coord, incremented of previous border and loc_vol
  int full_lx_of_coords(coords ext_x)
  {
    int ort_dir_bord[4][3]={{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
#if NDIM != 4
    #error puppa
#endif
    
    //pseudo-localize it
    coords x;
    for(int mu=0;mu<NDIM;mu++)
      {
	x[mu]=ext_x[mu];
	while(x[mu]<0) x[mu]+=glb_size[mu];
	while(x[mu]>=glb_size[mu]) x[mu]-=glb_size[mu];
      }
    
    //check locality
    int isloc=1;
    for(int mu=0;mu<NDIM;mu++)
      {
	isloc&=(x[mu]>=0);
	isloc&=(x[mu]<loc_size[mu]);
      }
    
    if(isloc) return loclx_of_coord(x);
    
    //check borderity
    int isbord[NDIM];
    for(int mu=0;mu<NDIM;mu++)
      {
	isbord[mu]=0;
	if(paral_dir[mu])
	  {
	    if(x[mu]==glb_size[mu]-1) isbord[mu]=-1;
	    if(x[mu]==loc_size[mu]) isbord[mu]=+1;
	  }
      }
    
    //check if it is in one of the NDIM forward or backward borders
#if NDIM != 4
    #error puppa
#endif
    for(int mu=0;mu<NDIM;mu++)
      if((isbord[ort_dir_bord[mu][0]]==0)&&(isbord[ort_dir_bord[mu][1]]==0)&&(isbord[ort_dir_bord[mu][2]]==0))
	{
	  if(isbord[mu]==-1) return loc_vol+bord_offset[mu]+bordlx_of_coord(x,mu);             //backward border comes first
	  if(isbord[mu]==+1) return loc_vol+bord_vol/2+bord_offset[mu]+bordlx_of_coord(x,mu);  //forward border comes after
	}
    
    return -1;
  }
  
  //return the border site adiacent at surface
  int bordlx_of_surflx(int loclx,int mu)
  {
    if(!paral_dir[mu]) return -1;
    if(loc_size[mu]<2) crash("not working if one dir is smaller than 2");
    
    if(loc_coord_of_loclx[loclx][mu]==0) return loclx_neighdw[loclx][mu]-loc_vol;
    if(loc_coord_of_loclx[loclx][mu]==loc_size[mu]-1) return loclx_neighup[loclx][mu]-loc_vol;
    
    return -1;
  }
  
  //label all the sites: bulk and border
  void label_all_sites()
  {
    coords x;
    
#if NDIM != 4
    #error puppa
#endif
    for(x[0]=-paral_dir[0];x[0]<loc_size[0]+paral_dir[0];x[0]++) 
      for(x[1]=-paral_dir[1];x[1]<loc_size[1]+paral_dir[1];x[1]++) 
	for(x[2]=-paral_dir[2];x[2]<loc_size[2]+paral_dir[2];x[2]++) 
	  for(x[3]=-paral_dir[3];x[3]<loc_size[3]+paral_dir[3];x[3]++) 
	    {
	      //check if it is defined
	      int iloc=full_lx_of_coords(x);
	      if(iloc!=-1)
		{
		  //compute global coordinates, assigning
		  for(int nu=0;nu<NDIM;nu++)
		    glb_coord_of_loclx[iloc][nu]=(x[nu]+rank_coord[nu]*loc_size[nu]+glb_size[nu])%glb_size[nu];
		  
		  //find the global index
		  int iglb=glblx_of_coord(glb_coord_of_loclx[iloc]);
		  
		  //if it is on the bulk store it
		  if(iloc<loc_vol)
		    {
		      for(int nu=0;nu<NDIM;nu++) loc_coord_of_loclx[iloc][nu]=x[nu];
		      glblx_of_loclx[iloc]=iglb;
		    }
		  
		  //if it is on the border store it
		  if(iloc>=loc_vol&&iloc<loc_vol+bord_vol)
		    {
		      int ibord=iloc-loc_vol;
		      glblx_of_bordlx[ibord]=iglb;
		      loclx_of_bordlx[ibord]=iloc;
		    }
		}
	    }
  }
  
  //find the neighbours
  void find_neighbouring_sites()
  {
    //loop over the four directions
    for(int ivol=0;ivol<loc_vol+bord_vol;ivol++)
      for(int mu=0;mu<NDIM;mu++)
	{
	  //copy the coords
	  coords n;
	  for(int nu=0;nu<NDIM;nu++) n[nu]=glb_coord_of_loclx[ivol][nu]-loc_size[nu]*rank_coord[nu];
	  
	  //move forward
	  n[mu]++;
	  int nup=full_lx_of_coords(n);
	  //move backward
	  n[mu]-=2;
	  int ndw=full_lx_of_coords(n);
	  
	  //if "local" assign it (automatically -1 otherwise)
	  loclx_neighup[ivol][mu]=nup;
	  loclx_neighdw[ivol][mu]=ndw;
	}
  }
  
  //finds how to fill the borders with opposite surface (up b->dw s)
  void find_surf_of_bord()
  {
    BISSA_LOC_VOL_LOOP(loclx)
      for(int mu=0;mu<NDIM;mu++)
	{
	  int bordlx=bordlx_of_surflx(loclx,mu);
	  if(bordlx!=-1) surflx_of_bordlx[bordlx]=loclx;
	}
  }
  
  //index all the sites on bulk
  void find_bulk_sites()
  {
    //check surfacity
    int ibulk=0,inon_fw_surf=0,inon_bw_surf=0;
    int isurf=0,ifw_surf=0,ibw_surf=0;
    BISSA_LOC_VOL_LOOP(ivol)
    {
      //find if it is on bulk or non_fw or non_bw surf
      int is_bulk=true,is_non_fw_surf=true,is_non_bw_surf=true;
      for(int mu=0;mu<NDIM;mu++)
	if(paral_dir[mu])
	  {
	    if(loc_coord_of_loclx[ivol][mu]==loc_size[mu]-1) is_bulk=is_non_fw_surf=false;
	    if(loc_coord_of_loclx[ivol][mu]==0)              is_bulk=is_non_bw_surf=false;
	  }
      
      //mark it
      if(is_bulk) loclx_of_bulklx[ibulk++]=ivol;
      else        loclx_of_surflx[isurf++]=ivol;
      if(is_non_fw_surf) loclx_of_non_fw_surflx[inon_fw_surf++]=ivol;
      else               loclx_of_fw_surflx[ifw_surf++]=ivol;
      if(is_non_bw_surf) loclx_of_non_bw_surflx[inon_bw_surf++]=ivol;
      else               loclx_of_bw_surflx[ibw_surf++]=ivol;
    }
    
    if(ibulk!=bulk_vol) crash("mismatch in bulk id");
    if(isurf!=surf_vol) crash("mismatch in surf id");
    if(inon_fw_surf!=non_fw_surf_vol) crash("mismatch in non_fw_surf id");
    if(inon_bw_surf!=non_bw_surf_vol) crash("mismatch in non_bw_surf id");
    if(ifw_surf!=fw_surf_vol) crash("mismatch in fw_surf id");
    if(ibw_surf!=bw_surf_vol) crash("mismatch in bw_surf id");
  }  
  
  //indexes run as t,x,y,z (faster:z)
  void set_lx_geometry()
  {
    if(lx_geom_inited==1) crash("cartesian geometry already intialized!");
    lx_geom_inited=1;
    
    if(grid_inited!=1) crash("grid not initialized!");
    
    //find the rank of the neighbour in the various dir
    for(int mu=0;mu<NDIM;mu++)
      MPI_Cart_shift(cart_comm,mu,1,&(rank_neighdw[mu]),&(rank_neighup[mu]));
    memcpy(rank_neigh[0],rank_neighdw,sizeof(coords));
    memcpy(rank_neigh[1],rank_neighup,sizeof(coords));
    
    loc_coord_of_loclx=bissa_malloc("loc_coord_of_loclx",loc_vol,coords);
    glb_coord_of_loclx=bissa_malloc("glb_coord_of_loclx",loc_vol+bord_vol,coords);
    loclx_neigh[0]=loclx_neighdw=bissa_malloc("loclx_neighdw",loc_vol+bord_vol,coords);
    loclx_neigh[1]=loclx_neighup=bissa_malloc("loclx_neighup",loc_vol+bord_vol,coords);  
    ignore_borders_communications_warning(loc_coord_of_loclx);
    ignore_borders_communications_warning(glb_coord_of_loclx);
    ignore_borders_communications_warning(loclx_neighup);
    ignore_borders_communications_warning(loclx_neighdw);
    
    //local to global
    glblx_of_loclx=bissa_malloc("glblx_of_loclx",loc_vol,int);
    
    //borders
    glblx_of_bordlx=bissa_malloc("glblx_of_bordlx",bord_vol,int);
    loclx_of_bordlx=bissa_malloc("loclx_of_bordlx",bord_vol,int);
    surflx_of_bordlx=bissa_malloc("surflx_of_bordlx",bord_vol,int);
    
    //bulk and surfs
    loclx_of_bulklx=bissa_malloc("loclx_of_bulklx",bulk_vol,int);
    loclx_of_surflx=bissa_malloc("loclx_of_surflx",surf_vol,int);
    loclx_of_non_bw_surflx=bissa_malloc("loclx_of_non_bw_surflx",non_bw_surf_vol,int);
    loclx_of_non_fw_surflx=bissa_malloc("loclx_of_non_fw_surflx",non_fw_surf_vol,int);
    loclx_of_bw_surflx=bissa_malloc("loclx_of_bw_surflx",bw_surf_vol,int);
    loclx_of_fw_surflx=bissa_malloc("loclx_of_fw_surflx",fw_surf_vol,int);
    
    //label the sites and neighbours
    label_all_sites();
    find_neighbouring_sites();
    
    //matches surface and opposite border
    find_surf_of_bord();
    
    //find bulk sites
    find_bulk_sites();
    
    //init sender and receiver points for borders
    for(int mu=0;mu<NDIM;mu++)
      if(paral_dir[mu]!=0)
	{
	  coords zero;
	  for(int nu=0;nu<NDIM;nu++) zero[nu]=0;
	  start_lx_bord_send_up[mu]=loclx_of_coord(zero);
	  start_lx_bord_rece_up[mu]=(loc_vol+bord_offset[mu]+bord_vol/2);
	  coords x;
	  for(int nu=0;nu<NDIM;nu++)
	    if(nu==mu) x[nu]=loc_size[mu]-1;
	    else x[nu]=0;
	  start_lx_bord_send_dw[mu]=loclx_of_coord(x);
	  start_lx_bord_rece_dw[mu]=loc_vol+bord_offset[mu];
	}
    
    //allocate a buffer large enough to allow communications of double lx border
    recv_buf_size=std::max(recv_buf_size,bord_vol*sizeof(double));
    send_buf_size=std::max(send_buf_size,bord_vol*sizeof(double));
    
    master_printf("Cartesian geometry intialized\n");
  }
  
  //global movements
  int glblx_neighup(int gx,int mu)
  {
    coords c;
    glb_coord_of_glblx(c,gx);
    c[mu]=(c[mu]+1)%glb_size[mu];
    
    return glblx_of_coord(c);
  }
  int glblx_neighdw(int gx,int mu)
  {
    coords c;
    glb_coord_of_glblx(c,gx);
    c[mu]=(c[mu]+glb_size[mu]-1)%glb_size[mu];
    
    return glblx_of_coord(c);
  }
  
  //wrapper for a previous defined function
  void get_loclx_and_rank_of_glblx(int *lx,int *rx,int gx)
  {
    coords c;
    glb_coord_of_glblx(c,gx);
    get_loclx_and_rank_of_coord(lx,rx,c);
  }
  
  //unset cartesian geometry
  void unset_lx_geometry()
  {
    if(lx_geom_inited!=1) crash("cartesian geometry not initialized!");
    
    master_printf("Unsetting cartesian geometry\n");
    lx_geom_inited=0;
    
#if defined BGQ && defined SPI
    free(recv_buf);
    free(send_buf);
#else
    bissa_free(recv_buf);
    bissa_free(send_buf);
#endif
    
    bissa_free(loc_coord_of_loclx);
    bissa_free(glb_coord_of_loclx);
    bissa_free(loclx_neighup);
    bissa_free(loclx_neighdw);
    
    bissa_free(glblx_of_loclx);
    bissa_free(glblx_of_bordlx);
    bissa_free(loclx_of_bordlx);
    bissa_free(surflx_of_bordlx);
    
    bissa_free(loclx_of_bulklx);
    bissa_free(loclx_of_surflx);
    bissa_free(loclx_of_non_fw_surflx);
    bissa_free(loclx_of_fw_surflx);
    bissa_free(loclx_of_non_bw_surflx);
    bissa_free(loclx_of_bw_surflx);
  }
  
  //definitions of lexical ordered sender for borders
  void initialize_lx_bord_senders_of_kind(MPI_Datatype *MPI_BORD_SEND,MPI_Datatype *base)
  {
#if NDIM != 4
    #error puppa
#endif
    //Various type useful for sub-borders
    MPI_Datatype MPI_3_SLICE;
    MPI_Datatype MPI_23_SLICE;
    MPI_Type_contiguous(loc_size[3],*base,&MPI_3_SLICE);
    MPI_Type_contiguous(loc_size[2]*loc_size[3],*base,&MPI_23_SLICE);
    
    ///////////define the sender for the NDIM kinds of borders////////////
    MPI_Type_contiguous(loc_size[1]*loc_size[2]*loc_size[3],*base,&(MPI_BORD_SEND[0]));
    MPI_Type_vector(loc_size[0],1,loc_size[1],MPI_23_SLICE,&(MPI_BORD_SEND[1]));
    MPI_Type_vector(loc_size[0]*loc_size[1],1,loc_size[2],MPI_3_SLICE,&(MPI_BORD_SEND[2]));
    MPI_Type_vector(loc_size[0]*loc_size[1]*loc_size[2],1,loc_size[3],*base,&(MPI_BORD_SEND[3]));
    //Commit
    for(int ibord=0;ibord<NDIM;ibord++) MPI_Type_commit(&(MPI_BORD_SEND[ibord]));
  }
  
  //definitions of lexical ordered receivers for borders
  void initialize_lx_bord_receivers_of_kind(MPI_Datatype *MPI_BORD_RECE,MPI_Datatype *base)
  {
#if NDIM != 4
    #error puppa
#endif
    //define the NDIM dir borders receivers, which are contiguous in memory
    MPI_Type_contiguous(loc_size[1]*loc_size[2]*loc_size[3],*base,&(MPI_BORD_RECE[0]));
    MPI_Type_contiguous(loc_size[0]*loc_size[2]*loc_size[3],*base,&(MPI_BORD_RECE[1]));
    MPI_Type_contiguous(loc_size[0]*loc_size[1]*loc_size[3],*base,&(MPI_BORD_RECE[2]));
    MPI_Type_contiguous(loc_size[0]*loc_size[1]*loc_size[2],*base,&(MPI_BORD_RECE[3]));
    for(int ibord=0;ibord<NDIM;ibord++) MPI_Type_commit(&(MPI_BORD_RECE[ibord]));
  }
  
  //initalize senders and receivers for borders of lexically ordered vectors
  void set_lx_bord_senders_and_receivers(MPI_Datatype *MPI_BORD_SEND,MPI_Datatype *MPI_BORD_RECE,MPI_Datatype *base)
  {
    initialize_lx_bord_senders_of_kind(MPI_BORD_SEND,base);
    initialize_lx_bord_receivers_of_kind(MPI_BORD_RECE,base);
  }
  
  //define all the local lattice momenta
  void define_local_momenta(momentum_t *k,double *k2,momentum_t *ktilde,double *ktilde2,momentum_t bc)
  {
    if(!lx_geom_inited) set_lx_geometry();
    
    //first of all, defines the local momenta for the various directions
    BISSA_LOC_VOL_LOOP(imom)
    {
      k2[imom]=ktilde2[imom]=0;
      for(int mu=0;mu<NDIM;mu++)
	{
	  k[imom][mu]=M_PI*(2*glb_coord_of_loclx[imom][mu]+bc[mu])/glb_size[mu];
	  ktilde[imom][mu]=sin(k[imom][mu]);
	  
	  k2[imom]+=k[imom][mu]*k[imom][mu];
	  ktilde2[imom]+=ktilde[imom][mu]*ktilde[imom][mu];
	}
    }
  }
}
