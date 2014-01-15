#include "include.h"

const int njacks=16;
int deb=0;
int T,TH,L,t0,tfit_op;
int parity[100];
int nlevls,nlevls_sto;
int reorder;
int ibeta;
double mc;
char infile[100];

jack F(jack M,jack Z)
{
  return fabs(Z*2*mc/sinh(M)/M);
}

int main(int narg,char **arg)
{
  if(narg<2) crash("use %s input" ,arg[0]);
  FILE *fin=open_file(arg[1],"r");
  
  //size
  read_formatted_from_file_expecting((char*)&T,fin,"%d","T");
  L=TH=T/2;
  
  //beta and charm mass
  read_formatted_from_file_expecting((char*)&ibeta,fin,"%d","ibeta");
  read_formatted_from_file_expecting((char*)&mc,fin,"%lg","mc");
  
  //file path
  read_formatted_from_file_expecting(infile,fin,"%s","infile");
  
  //read the stored number of nlevls
  read_formatted_from_file_expecting((char*)&nlevls_sto,fin,"%d","nlevls_sto");

  //allocate nlevls-depending stuff
  read_formatted_from_file_expecting((char*)&nlevls,fin,"%d","nlevls");
  int map[nlevls];
  for(int ilev=0;ilev<nlevls;ilev++) read_formatted_from_file((char*)&(map[ilev]),fin,"%d","lev");
  
  //timeslice for normalization
  read_formatted_from_file_expecting((char*)&t0,fin,"%d","t0");
  read_formatted_from_file_expecting((char*)&tfit_op,fin,"%d","tfit_op");
  
  //read wheter reorder or not, and debug
  read_formatted_from_file_expecting((char*)&reorder,fin,"%d","reorder");
  printf("Reordering: %d\n",reorder);
  read_formatted_from_file_expecting((char*)&deb,fin,"%d","debug");
  printf("Debug: %d\n",deb);
  
  //interval to fit ground state
  int tfit_ground_min,tfit_ground_max;
  read_formatted_from_file_expecting((char*)&tfit_ground_min,fin,"%d","tfit_ground");
  read_formatted_from_file((char*)&tfit_ground_max,fin,"%d","tfit_ground");
  
  //interval to fit first state
  int tfit_first_min,tfit_first_max;
  read_formatted_from_file_expecting((char*)&tfit_first_min,fin,"%d","tfit_first");
  read_formatted_from_file((char*)&tfit_first_max,fin,"%d","tfit_first");
  
  //interval to fit first state
  int tfit_second_min,tfit_second_max;
  read_formatted_from_file_expecting((char*)&tfit_second_min,fin,"%d","tfit_second");
  read_formatted_from_file((char*)&tfit_second_max,fin,"%d","tfit_second");
  
  fclose(fin);
  
  //load data
  gevp_pars_t g(nlevls,njacks,TH,t0);
  g.load_raw_data("raw_data.xmg",infile,map,nlevls_sto,0);
  
  ////////////////////////////// finished reading input ///////////////////////////
  
  //resolve gevp
  g.gevp();
  g.check_orthogonality();
  
  //reorder
  if(reorder) g.reorder_eig();
  g.convert_to_full_eig_ve();
  
  for(int t=0;t<TH;t++)
    {
      for(int ilev=0;ilev<nlevls;ilev++) cout<<smart_print(g.eig_va[ilev][t])<<" ";
      cout<<endl;
    }
  cout<<endl;
  
  {
    ofstream out("gevp.xmg");
    out<<"@type xydy"<<endl;
    for(int ilev=0;ilev<nlevls;ilev++)
      {
	jvec m=effective_mass(g.eig_va[ilev]);
	out<<m<<endl;
	out<<"&"<<endl;
      }
  }
  
  jvec M_P(nlevls,njacks),Z_P(nlevls,njacks);
  {
    ofstream fits("fits.xmg");
    jack Z2(njacks);
    
    two_pts_fit(M_P[0],Z2,g.eig_va[0],tfit_ground_min,tfit_ground_max);
    Z_P[0]=sqrt(Z2);
    fits<<write_constant_fit_plot(effective_mass(g.eig_va[0]),M_P[0],tfit_ground_min,tfit_ground_max,0);
    cout<<"M0: "<<smart_print(M_P[0])<<endl;

    if(nlevls>1)
      {
	two_pts_fit(M_P[1],Z2,g.eig_va[1],tfit_first_min,tfit_first_max);
	Z_P[1]=sqrt(Z2);
	fits<<write_constant_fit_plot(effective_mass(g.eig_va[1]),M_P[1],tfit_first_min,tfit_first_max,3);
	
	cout<<"M1/M0: "<<smart_print(M_P[1]/M_P[0])<<endl;
	if(nlevls>2)
	  {
	    two_pts_fit(M_P[2],Z2,g.eig_va[2],tfit_second_min,tfit_second_max);
	    Z_P[2]=sqrt(Z2);
	    fits<<write_constant_fit_plot(effective_mass(g.eig_va[2]),M_P[2],tfit_second_min,tfit_second_max,6);
	    
	    cout<<"M2/M0: "<<smart_print(M_P[2]/M_P[0])<<endl;
	  }
      }
  }
  
  //define a new map or copy old
  int *map2;
  map2=new int[nlevls+1];
  map2[0]=0;
  for(int ilev=0;ilev<nlevls;ilev++) map2[ilev+1]=map[ilev];
  
  //load also lev 0
  gevp_pars_t g2(nlevls+1,njacks,TH,t0);
  g2.load_raw_data(NULL,infile,map2,nlevls_sto,0);
  
  //build sink opt corr
  jvec impr_sink_corr[nlevls];
  for(int ilev=0;ilev<nlevls;ilev++)
    {
      impr_sink_corr[ilev]=jvec(TH+1,njacks);
      impr_sink_corr[ilev]*=0;
      
      for(int iop=0;iop<nlevls;iop++)
	{
	  jack w=g.full_eig_ve[iop*nlevls+ilev][tfit_op];
	  impr_sink_corr[ilev]+=g2.data[iop+1]*w;
	}
    }
  
  //plot coeff
  for(int ilev=0;ilev<nlevls;ilev++)
    {
      ofstream out_coeff_plot(combine("out_coeff_%d_plot.xmg",ilev).c_str());
      out_coeff_plot<<"@type xydy"<<endl;
      for(int iop=0;iop<nlevls;iop++)
	{
	  for(int t=1;t<=TH;t++) out_coeff_plot<<g.eig_ve[iop*nlevls+ilev][t]<<endl;
	  out_coeff_plot<<"&"<<endl;
	}
      out_coeff_plot.close();
    }
  
  //plot fit
  ofstream out_sink_plot("out_sink_plot.xmg");
  out_sink_plot<<"@type xydy"<<endl;
  
  jack ZLZOPT(njacks);
  jvec ZL(nlevls,njacks),M_P2(nlevls,njacks);

  two_pts_fit(M_P2[0],ZLZOPT,impr_sink_corr[0],tfit_ground_min,tfit_ground_max);
  ZL[0]=ZLZOPT/Z_P[0];
  double lat_med[4]={0.486508,0.422773,0.335339,0.268402};
  cout<<"M0: "<<smart_print(M_P2[0])<<endl;
  jack FETA=F(M_P2[0],ZL[0]);
  FETA.write_to_binfile("FETA");
  cout<<"F[0]: "<<smart_print(FETA/lat_med[ibeta])<<endl;
  out_sink_plot<<write_constant_fit_plot(effective_mass(impr_sink_corr[0]),M_P2[0],tfit_ground_min,tfit_ground_max,0);

  if(nlevls>1)
    {
      two_pts_fit(M_P2[1],ZLZOPT,impr_sink_corr[1],tfit_first_min,tfit_first_max);
      ZL[1]=ZLZOPT/Z_P[1];
      jack FETA_P=F(M_P2[1],ZL[1]);
      jack FETA_P_FETA=FETA_P/FETA;
      FETA_P_FETA.write_to_binfile("FETA_P_FETA");
      cout<<"F[1]: "<<smart_print(FETA_P/lat_med[ibeta])<<endl;
      cout<<"F[1]/F[0]: "<<smart_print(FETA_P_FETA)<<endl;
      out_sink_plot<<write_constant_fit_plot(effective_mass(impr_sink_corr[1]),M_P2[1],tfit_first_min,
					     tfit_first_max,3);
      
      cout<<"M1/M0: "<<smart_print(M_P2[1]/M_P2[0])<<endl;
      if(nlevls>2)
	{
	  two_pts_fit(M_P2[2],ZLZOPT,impr_sink_corr[2],tfit_second_min,tfit_second_max);
	  ZL[2]=ZLZOPT/Z_P[2];
	  jack FETA_PP=F(M_P2[2],ZL[2]);
	  jack FETA_PP_FETA=FETA_PP/FETA;
	  FETA_PP_FETA.write_to_binfile("FETA_PP_FETA");
	  cout<<"F[2]: "<<smart_print(FETA_PP/lat_med[ibeta])<<endl;
	  cout<<"F[2]/F[0]: "<<smart_print(FETA_PP_FETA)<<endl;
	  out_sink_plot<<write_constant_fit_plot(effective_mass(impr_sink_corr[2]),M_P2[2],
						 tfit_second_min,tfit_second_max,6);
	  
	  cout<<"M2/M0: "<<smart_print(M_P2[2]/M_P2[0])<<endl;
	}
    }
  
  return 0;
}
