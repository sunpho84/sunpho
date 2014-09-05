#include "include.h"

const int njacks=16;
int deb=0;
int T,TH,L,t0,tsep,ibeta;
int nlevls,nlevls_sto;
int reorder;
char infile[100];
double ZV_mom[4]={0.5816,0.6103,0.6451,0.686};
double lat_med[4]={0.486508,0.422773,0.335339,0.268402};
double th[4]={1.10,1,1.04,0.98};
int *op_map;
int tfit_impr_min,tfit_impr_max;
int tfit_3pts_fexc_min,tfit_3pts_fexc_max;
int tfit_PP_min,tfit_PP_max;
int tfit_VV_min,tfit_VV_max;
int tfit_op;
jvec VKVK_st;
jvec VKVK_mv;
gevp_pars_t *g;

double momentum(double th)
{return M_PI/TH*th;}

//relativistic lattice energy: lattice...
jack latt_en(jack M,double th)
{return 2*asinh(sqrt(3*sqr(sin(momentum(th)/2))+sqr(sinh(M/2))));}

void scan_input(const char *path)
{
  FILE *fin=open_file(path,"r");
  
  //size
  read_formatted_from_file_expecting((char*)&T,fin,"%d","T");
  L=TH=T/2;
  read_formatted_from_file_expecting((char*)&tsep,fin,"%d","tsep");
  
  //ibeta
  read_formatted_from_file_expecting((char*)&ibeta,fin,"%d","ibeta");
  
  //file path
  read_formatted_from_file_expecting(infile,fin,"%s","infile");
  
  //read the stored number of nlevls
  read_formatted_from_file_expecting((char*)&nlevls_sto,fin,"%d","nlevls_sto");

  //allocate nlevls-depending stuff
  read_formatted_from_file_expecting((char*)&nlevls,fin,"%d","nlevls");
  op_map=(int*)malloc(nlevls*sizeof(int));
  for(int ilev=0;ilev<nlevls;ilev++) read_formatted_from_file((char*)&(op_map[ilev]),fin,"%d","lev");
  
  //timeslice for normalization
  read_formatted_from_file_expecting((char*)&t0,fin,"%d","t0");
  read_formatted_from_file_expecting((char*)&tfit_op,fin,"%d","tfit_op");
  
  //read wheter reorder or not, and debug
  read_formatted_from_file_expecting((char*)&reorder,fin,"%d","reorder");
  printf("Reordering: %d\n",reorder);
  read_formatted_from_file_expecting((char*)&deb,fin,"%d","debug");
  printf("Debug: %d\n",deb);
  
  //interval to fit improved operators
  read_formatted_from_file_expecting((char*)&tfit_impr_min,fin,"%d","tfit_impr");
  read_formatted_from_file((char*)&tfit_impr_max,fin,"%d","tfit_impr");

  //interval to fit non improved one
  read_formatted_from_file_expecting((char*)&tfit_PP_min,fin,"%d","tfit_PP");
  read_formatted_from_file((char*)&tfit_PP_max,fin,"%d","tfit_PP");
  
  //interval to fit V
  read_formatted_from_file_expecting((char*)&tfit_VV_min,fin,"%d","tfit_VV");
  read_formatted_from_file((char*)&tfit_VV_max,fin,"%d","tfit_VV");
  
  //interval to fit first excited matrix element
  read_formatted_from_file_expecting((char*)&tfit_3pts_fexc_min,fin,"%d","tfit_3pts_fexc");
  read_formatted_from_file((char*)&tfit_3pts_fexc_max,fin,"%d","tfit_3pts_fexc");
  
  fclose(fin);
  
  /////////////////////////////////// read data/////////////////////////////////////

  //load standing V
  VKVK_st=jvec_load(infile,T,njacks,0).simmetrized(1);
  VKVK_mv=jvec_load(infile,T,njacks,1).simmetrized(1);
  
  //load data for gevp
  g=new gevp_pars_t(nlevls,njacks,TH,t0);
  g->load_raw_data("raw_data.xmg",infile,op_map,nlevls_sto,2);
}
