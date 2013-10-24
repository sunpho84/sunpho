#include "include.h"

const int njack=16;
int ibeta,nmass,nth,ncombo;
int T,TH,tsep,tmin_2pts,tmax_2pts,tmin_3pts,tmax_3pts;
double *theta,*mass,*q2;

///////////////////////////////////////////////////////////////////////////////////

int icombo_2pts(int r,int im1,int im2,int ith)
{return im1+nmass*(r+2*(im2+nmass*ith));}

int icombo_3pts(int im_S1,int ith_S0)
{return ith_S0+nth*im_S1;}

jvec load_2pts(const char *path,int im1,int im2,int ith)
{return (jvec_load(path,T,njack,icombo_2pts(0,im1,im2,ith))+jvec_load(path,T,njack,icombo_2pts(1,im1,im2,ith)))/2;}

jvec load_3pts(const char *name,int im_S1,int ith_S0,int parity)
{
  jvec a=-jvec_load(combine("../CORRELATORS/3pts_%s_sp0_30_30",name).c_str(),T,njack,icombo_3pts(im_S1,ith_S0));
  
  return a.simmetrized(parity);
}

void read_set_pars(const char *path)
{
  FILE *set_pars_file=open_file(path,"r");
  
  read_formatted_from_file_expecting((char*)&T,set_pars_file,"%d","T");
  TH=T/2;
  read_formatted_from_file_expecting((char*)&tsep,set_pars_file,"%d","TSep");
  
  //beta
  read_formatted_from_file_expecting((char*)&ibeta,set_pars_file,"%d","ibeta");
  
  //read the masses
  read_formatted_from_file_expecting((char*)&nmass,set_pars_file,"%d","nmass");
  expect_string_from_file(set_pars_file,"mass_list");
  mass=(double*)malloc(sizeof(double)*nmass);
  for(int imass=0;imass<nmass;imass++) read_formatted_from_file((char*)&(mass[imass]),set_pars_file,"%lg","mass");
  
  //read theta
  read_formatted_from_file_expecting((char*)&nth,set_pars_file,"%d","ntheta");
  expect_string_from_file(set_pars_file,"theta_list");
  theta=(double*)malloc(sizeof(double)*nth);
  q2=(double*)malloc(sizeof(double)*nth);
  for(int ith=0;ith<nth;ith++)
    {
      read_formatted_from_file((char*)&(theta[ith]),set_pars_file,"%lg","theta");
      q2[ith]=3*sqr(M_PI/TH*theta[ith]);
    }
  
  ncombo=nth*nmass;
  
  fclose(set_pars_file);
}

void read_input(const char *path)
{
  FILE *input_file=open_file(path,"r");
  
  read_formatted_from_file_expecting((char*)&tmin_2pts,input_file,"%d","tint_2pts");
  read_formatted_from_file((char*)&tmax_2pts,input_file,"%d","tint_2pts");

  read_formatted_from_file_expecting((char*)&tmin_3pts,input_file,"%d","tint_3pts");
  read_formatted_from_file((char*)&tmax_3pts,input_file,"%d","tint_3pts");
  
  fclose(input_file);
}



int main()
{
  read_set_pars("../data_pars");
  read_input("analysis_pars");
  
  ////////////////////////// P5 /////////////////////
  
  jvec E(ncombo,njack),ZS(ncombo,njack),ZL(ncombo,njack);
  jvec P5_corr_SL[ncombo],P5_corr_SS[ncombo];
  ofstream raw_corr_plot("plots/raw_corr_P5.xmg");
  ofstream fit_corr_plot_SL("plots/fit_corr_P5_SL.xmg");
  ofstream fit_corr_plot_SS("plots/fit_corr_P5_SS.xmg");
  raw_corr_plot<<"@type xydy"<<endl;
  fit_corr_plot_SL<<"@type xydy"<<endl;
  fit_corr_plot_SS<<"@type xydy"<<endl;
  for(int ith=0;ith<nth;ith++)
    for(int im_S0=0;im_S0<nmass;im_S0++)
      {
	int icombo=ith*nmass+im_S0;
	
	//load P5
	jvec SL=P5_corr_SL[icombo]=load_2pts("../CORRELATORS/2pts_P5P5_30_00",im_S0,0,ith).simmetrized(1);
	jvec SS=P5_corr_SS[icombo]=load_2pts("../CORRELATORS/2pts_P5P5_30_30",im_S0,0,ith).simmetrized(1);
	
	//fit P5 Kaon
	two_pts_SL_fit(E[icombo],ZL[icombo],ZS[icombo],SL,SS,tmin_2pts,tmax_2pts,tmin_2pts,tmax_2pts);
	
	raw_corr_plot<<SL<<"&"<<endl<<SS<<"&"<<endl;
	fit_corr_plot_SL<<write_constant_fit_plot(effective_mass(SL),E[icombo],tmin_2pts,tmax_2pts,3*icombo);
	fit_corr_plot_SS<<write_constant_fit_plot(effective_mass(SS),E[icombo],tmin_2pts,tmax_2pts,3*icombo);
      }
  
  //plot dispertion relations
  ofstream E_plot("plots/E.xmg");
  ofstream ZL_plot("plots/ZL.xmg");
  ofstream ZS_plot("plots/ZS.xmg");
  E_plot<<"@type xydy"<<endl;
  ZL_plot<<"@type xydy"<<endl;
  ZS_plot<<"@type xydy"<<endl;
  for(int im_S0=0;im_S0<nmass;im_S0++)
    {
      for(int ith=0;ith<nth;ith++)
	{
	  int icombo=ith*nmass+im_S0;
	  E_plot<<q2[ith]<<" "<<E[icombo]<<endl;
	  ZL_plot<<q2[ith]<<" "<<ZL[icombo]<<endl;
	  ZS_plot<<q2[ith]<<" "<<ZS[icombo]<<endl;
	}
      E_plot<<"&"<<endl;
      ZL_plot<<"&"<<endl;
      ZS_plot<<"&"<<endl;
    }
    
  //////////////////////////////////// P5 {V0,VK,TK} P5 ////////////////////////////////
  
  //compute matrix elements
  const char ME_flag[3][10]={"V0P5","VKP5","TKP5"};
  const int parity_ME[3]={-1,+1,-1};
  jvec ME(nth,njack);
  for(int im_S1=0;im_S1<nmass;im_S1++)
    {
      //open the output
      ofstream raw_corr_3pts_plot(combine("plots/raw_3pts_%d.xmg",im_S1).c_str());
      ofstream ME_plot[3];
      ME_plot[0].open(combine("plots/V0_corr_%d.xmg",im_S1).c_str());
      ME_plot[1].open(combine("plots/VK_corr_%d.xmg",im_S1).c_str());
      ME_plot[2].open(combine("plots/TK_corr_%d.xmg",im_S1).c_str());
      raw_corr_3pts_plot<<"@type xydy"<<endl;
      for(int i_ME=0;i_ME<3;i_ME++) ME_plot[i_ME]<<"@type xydy"<<endl;
	
      //loop over theta
      for(int ith_S0=0;ith_S0<nth;ith_S0++)
	{
	  //define the temporal dependance
	  jack E_eta=/*latt_en*/(E[ith_S0*nmass]/*,theta[ith_S0]*/),ZS_eta=ZS[ith_S0*nmass];
	  jack M_Bc=E[im_S1],ZS_Bc=ZS[im_S1];
	  jvec dT(tsep+1,njack);
	  for(int t=0;t<=tsep;t++)
	    dT[t]=(ZS_eta*ZS_Bc)/(2*E_eta*2*M_Bc)*exp(-E_eta*t)*exp(-M_Bc*(tsep-t));
	  
	  //load the corr and divide it by dT
	  for(int i_ME=0;i_ME<3;i_ME++)
	    {
	      jvec raw_corr_3pts=load_3pts(ME_flag[i_ME],im_S1,ith_S0,parity_ME[i_ME]);
	      raw_corr_3pts_plot<<raw_corr_3pts<<"&"<<endl;
	      
	      //divide
	      jvec ME_corr=raw_corr_3pts/dT;

	      ME_plot[i_ME]<<ME_corr<<"&"<<endl;
	    }
	}
    }
      
  return 0;
}
