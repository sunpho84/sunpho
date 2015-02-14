#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/global_variables.hpp"
#include "base/debug.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/new_types_definitions.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace bissa
{
  //compute the parity of a global site
  int glb_coord_parity(coords c)
  {
    int par=0;
    for(int mu=0;mu<4;mu++) par+=c[mu];
    par%=2;
    
    return par;
  }
  int glblx_parity(int glx)
  {
    coords c;
    glb_coord_of_glblx(c,glx);
    
    return glb_coord_parity(c);
  }
  
  //set the eo geometry
  void set_eo_geometry()
  {
    if(!use_eo_geom) crash("E/O Geometry was not to be used!");
    if(eo_geom_inited) crash("E/O Geometry already initialized!");
    
    //check that all local sizes are multiples of 2
    int ok=1;
    for(int mu=0;mu<4;mu++) ok&=(loc_size[mu]%2==0);
    if(!ok) crash("local lattice size odd!");
    
    //set half the vol and bord
    glb_volh=glb_vol/2;
    loc_volh=loc_vol/2;
    
    //set the various time-slice types
    loclx_parity=bissa_malloc("loclx_parity",loc_vol+bord_vol,int);
    ignore_borders_communications_warning(loclx_parity);
    
    loceo_of_loclx=bissa_malloc("loceo_of_loclx",loc_vol+bord_vol,int);
    ignore_borders_communications_warning(loceo_of_loclx);
    
    for(int par=0;par<2;par++) loclx_of_loceo[par]=bissa_malloc("loclx_of_loceo",loc_volh+bord_volh,int);
    for(int par=0;par<2;par++) loceo_neighup[par]=bissa_malloc("loceo_neighup",loc_volh+bord_volh,coords);
    for(int par=0;par<2;par++) loceo_neighdw[par]=bissa_malloc("loceo_neighdw",loc_volh+bord_volh,coords);
    for(int par=0;par<2;par++) surfeo_of_bordeo[par]=bissa_malloc("surfeo_of_bordeo",bord_volh,int);
    for(int par=0;par<2;par++) ignore_borders_communications_warning(loclx_of_loceo[par]);
    for(int par=0;par<2;par++) ignore_borders_communications_warning(loceo_neighup[par]);
    for(int par=0;par<2;par++) ignore_borders_communications_warning(loceo_neighdw[par]);
    
    //Label the sites
    int iloc_eo[2]={0,0};
    for(int loclx=0;loclx<loc_vol+bord_vol;loclx++)
      {
	//fix parity of local index
	int par=loclx_parity[loclx]=glb_coord_parity(glb_coord_of_loclx[loclx]);
	
	//associate the e/o index to lx sites and vice-versa
	loceo_of_loclx[loclx]=iloc_eo[par];
	loclx_of_loceo[par][iloc_eo[par]]=loclx;
	iloc_eo[par]++;
      }
    
    //Fix the movements among e/o ordered sites
    for(int loclx=0;loclx<loc_vol+bord_vol;loclx++)
      for(int mu=0;mu<4;mu++)
	{
	  //take parity and e/o corresponding site
	  int par=loclx_parity[loclx];
	  int loceo=loceo_of_loclx[loclx];
	  
	  //up movements
	  int loclx_up=loclx_neighup[loclx][mu];
	  if(loclx_up>=0 && loclx_up<loc_vol+bord_vol)
	    loceo_neighup[par][loceo][mu]=loceo_of_loclx[loclx_up];
	  
	  //dw movements
	  int loclx_dw=loclx_neighdw[loclx][mu];
	  if(loclx_dw>=0 && loclx_dw<loc_vol+bord_vol)
	    loceo_neighdw[par][loceo][mu]=loceo_of_loclx[loclx_dw];
	}
    
    //finds how to fill the borders with surface
    for(int bordlx=0;bordlx<bord_vol;bordlx++)
      {
	int surflx=surflx_of_bordlx[bordlx];
	surfeo_of_bordeo[loclx_parity[surflx]][loceo_of_loclx[bordlx+loc_vol]-loc_volh]=loceo_of_loclx[surflx];
      }
    
    //init sender and receiver points for borders
    for(int mu=0;mu<4;mu++)
      if(paral_dir[mu]!=0)
	{
	  start_eo_bord_send_up[mu]=loceo_of_loclx[start_lx_bord_send_up[mu]];
	  start_eo_bord_rece_up[mu]=loceo_of_loclx[start_lx_bord_rece_up[mu]];
	  start_eo_bord_send_dw[mu]=loceo_of_loclx[start_lx_bord_send_dw[mu]];
	  start_eo_bord_rece_dw[mu]=loceo_of_loclx[start_lx_bord_rece_dw[mu]];
	}
    
    master_printf("E/O Geometry intialized\n");
    
    eo_geom_inited=1;
  }
  
  //definitions of e/o split sender for borders
  void initialize_eo_bord_senders_of_kind(MPI_Datatype *MPI_EO_BORD_SEND_TXY,MPI_Datatype *MPI_EV_BORD_SEND_Z,MPI_Datatype *MPI_OD_BORD_SEND_Z,MPI_Datatype *base)
  {
    //Various type useful for sub-borders
    MPI_Datatype MPI_EO_3_SLICE;
    MPI_Datatype MPI_EO_23_SLICE;
    MPI_Type_contiguous(loc_size[3]/2,*base,&MPI_EO_3_SLICE);
    MPI_Type_contiguous(loc_size[2]*loc_size[3]/2,*base,&MPI_EO_23_SLICE);
    
    ///////////define the sender for the 4 kinds of borders////////////
    MPI_Type_contiguous(loc_size[1]*loc_size[2]*loc_size[3]/2,*base,&(MPI_EO_BORD_SEND_TXY[0]));
    MPI_Type_vector(loc_size[0],1,loc_size[1],MPI_EO_23_SLICE,&(MPI_EO_BORD_SEND_TXY[1]));
    MPI_Type_vector(loc_size[0]*loc_size[1],1,loc_size[2],MPI_EO_3_SLICE,&(MPI_EO_BORD_SEND_TXY[2]));
    //Commit
    for(int ibord=0;ibord<3;ibord++) MPI_Type_commit(&(MPI_EO_BORD_SEND_TXY[ibord]));
    
    //the z sending border is a mess
    int eo_bord_z_size=loc_volh/loc_size[3];
    int *ev_bord_z_pos_disp_dw=bissa_malloc("ev_bord_z_disp_dw",eo_bord_z_size,int);
    int *ev_bord_z_pos_disp_up=bissa_malloc("ev_bord_z_disp_up",eo_bord_z_size,int);
    int *od_bord_z_pos_disp_dw=bissa_malloc("od_bord_z_disp_dw",eo_bord_z_size,int);
    int *od_bord_z_pos_disp_up=bissa_malloc("od_bord_z_disp_up",eo_bord_z_size,int);
    int *single=bissa_malloc("single",eo_bord_z_size,int);
    int ev_izdw=0,ev_izup=0;
    int od_izdw=0,od_izup=0;
    BISSA_LOC_VOLH_LOOP(ieo)
    {
      int ev_ilx=loclx_of_loceo[0][ieo];
      int od_ilx=loclx_of_loceo[1][ieo];
      int ev_x3=loc_coord_of_loclx[ev_ilx][3];
      int od_x3=loc_coord_of_loclx[od_ilx][3];
      if(ev_x3==0) ev_bord_z_pos_disp_dw[ev_izdw++]=ieo;
      if(ev_x3==loc_size[3]-1) ev_bord_z_pos_disp_up[ev_izup++]=ieo;
      if(od_x3==0) od_bord_z_pos_disp_dw[od_izdw++]=ieo;
      if(od_x3==loc_size[3]-1) od_bord_z_pos_disp_up[od_izup++]=ieo;
    }
    for(int ibord_z=0;ibord_z<eo_bord_z_size;ibord_z++)
      single[ibord_z]=1;
    
    MPI_Type_indexed(eo_bord_z_size,single,ev_bord_z_pos_disp_dw,*base,&(MPI_EV_BORD_SEND_Z[0]));
    MPI_Type_indexed(eo_bord_z_size,single,ev_bord_z_pos_disp_up,*base,&(MPI_EV_BORD_SEND_Z[1]));
    
    MPI_Type_indexed(eo_bord_z_size,single,od_bord_z_pos_disp_dw,*base,&(MPI_OD_BORD_SEND_Z[0]));
    MPI_Type_indexed(eo_bord_z_size,single,od_bord_z_pos_disp_up,*base,&(MPI_OD_BORD_SEND_Z[1]));
    
    //commit the mess
    MPI_Type_commit(&(MPI_EV_BORD_SEND_Z[0]));
    MPI_Type_commit(&(MPI_EV_BORD_SEND_Z[1]));
    
    MPI_Type_commit(&(MPI_OD_BORD_SEND_Z[0]));
    MPI_Type_commit(&(MPI_OD_BORD_SEND_Z[1]));
    
    bissa_free(single);
    bissa_free(ev_bord_z_pos_disp_dw);
    bissa_free(ev_bord_z_pos_disp_up);
    bissa_free(od_bord_z_pos_disp_dw);
    bissa_free(od_bord_z_pos_disp_up);
  }
  
  //definitions of e/o split receivers for borders
  void initialize_eo_bord_receivers_of_kind(MPI_Datatype *MPI_EO_BORD_RECE,MPI_Datatype *base)
  {
    //define the 4 dir borders receivers, which are contiguous in memory
    MPI_Type_contiguous(loc_size[1]*loc_size[2]*loc_size[3]/2,*base,&(MPI_EO_BORD_RECE[0]));
    MPI_Type_contiguous(loc_size[0]*loc_size[2]*loc_size[3]/2,*base,&(MPI_EO_BORD_RECE[1]));
    MPI_Type_contiguous(loc_size[0]*loc_size[1]*loc_size[3]/2,*base,&(MPI_EO_BORD_RECE[2]));
    MPI_Type_contiguous(loc_size[0]*loc_size[1]*loc_size[2]/2,*base,&(MPI_EO_BORD_RECE[3]));
    for(int ibord=0;ibord<4;ibord++) MPI_Type_commit(&(MPI_EO_BORD_RECE[ibord]));
  }
  
  //initalize senders and receivers for borders of e/o split ordered vectors
  void set_eo_bord_senders_and_receivers(MPI_Datatype *MPI_EO_BORD_SEND_TXY,MPI_Datatype *MPI_EV_BORD_SEND_Z,MPI_Datatype *MPI_OD_BORD_SEND_Z,MPI_Datatype *MPI_EO_BORD_RECE,MPI_Datatype *base)
  {
    initialize_eo_bord_senders_of_kind(MPI_EO_BORD_SEND_TXY,MPI_EV_BORD_SEND_Z,MPI_OD_BORD_SEND_Z,base);
    initialize_eo_bord_receivers_of_kind(MPI_EO_BORD_RECE,base);
  }
  
  //unset the eo geometry
  void unset_eo_geometry()
  {
    if(!eo_geom_inited)
      crash("asking to unset never initialized E/O Geometry!");
    
    master_printf("Unsetting E/O Geometry\n");
    
    for(int par=0;par<2;par++)
      {
	bissa_free(loclx_of_loceo[par]);
	bissa_free(surfeo_of_bordeo[par]);
	bissa_free(loceo_neighup[par]);
	bissa_free(loceo_neighdw[par]);
      }
    bissa_free(loclx_parity);
    bissa_free(loceo_of_loclx);
    
    eo_geom_inited=0;
  }
}
