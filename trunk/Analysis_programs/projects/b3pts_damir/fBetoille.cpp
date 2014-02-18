#include "common.cpp"

int tmin_K,tmax_K;
int tmin_H,tmax_H;
int tmin_V0,tmax_V0;
int tmin_VK,tmax_VK;
int tmin_TK,tmax_TK;

int im_spec,im_S0,im_S1;

///////////////////////////////////////////////////////////////////////////////////

int icombo_2pts(REIM ri,int r,S_2pts_mass im1,S_2pts_mass im2)
{return ri+2*(im1+nm*(r+2*im2));}

jvec load_2pts(const char *path,REIM ri,S_2pts_mass im1,S_2pts_mass im2)
{return (jvec_load(path,T,njack,icombo_2pts(ri,0,im1,im2))+jvec_load(path,T,njack,icombo_2pts(ri,1,im1,im2)))/2;}

int main(int narg,char **arg)
{
  int im_spec=0;

  FILE *input_file=open_file("analysis_pars","r");

  int tminP,tmaxP,tminV,tmaxV;
  read_formatted_from_file_expecting((char*)&tminP,input_file,"%d","tintP");
  read_formatted_from_file((char*)&tmaxP,input_file,"%d","tintP");
  read_formatted_from_file_expecting((char*)&tminV,input_file,"%d","tintV");
  read_formatted_from_file((char*)&tmaxV,input_file,"%d","tintV");
  
  fclose(input_file);
  
  ///////
  
  read_set_pars("../data_pars");

  ofstream P5_SL_out("M_P5_fit_SL.xmg");
  ofstream P5_SS_out("M_P5_fit_SS.xmg");
  ofstream VK_SL_out("M_VK_fit_SL.xmg");
  ofstream VK_SS_out("M_VK_fit_SS.xmg");
  
  //results
  jvec M_P5(10,njack),ZL_P5(10,njack);
  jvec M_VK(10,njack),ZL_VK(10,njack);

  ofstream table("table_for_Damir.txt");
  table.precision(16);
  for(int im_S0=0;im_S0<10;im_S0++)
    {
      //temp
      jack ZS_P5,ZS_VK;
      cout<<"Considering "<<prod_meson[im_spec][min(im_S0,1)]<<endl;
      
      //load P5 and fit
      jvec P5_corr_SL=load_2pts("../CORRELATORS/2pts_P5P5_30_00",RE,H_S_2PTS[im_S0],L_S_2PTS[im_spec]).simmetrized(1);
      jvec P5_corr_SS=load_2pts("../CORRELATORS/2pts_P5P5_30_30",RE,H_S_2PTS[im_S0],L_S_2PTS[im_spec]).simmetrized(1);
      two_pts_SL_fit(M_P5[im_S0],ZL_P5[im_S0],ZS_P5,P5_corr_SL,P5_corr_SS,tminP,tmaxP,tminP,tmaxP);
      P5_SL_out<<write_constant_fit_plot(effective_mass(P5_corr_SL),M_P5[im_S0],tminP,tmaxP,3*im_S0)<<endl;
      P5_SS_out<<write_constant_fit_plot(effective_mass(P5_corr_SS),M_P5[im_S0],tminP,tmaxP,3*im_S0)<<endl;
      
      //load VK and fit
      jvec VK_corr_SL=load_2pts("../CORRELATORS/2pts_VKVK_30_00",RE,H_S_2PTS[im_S0],L_S_2PTS[im_spec]).simmetrized(1);
      jvec VK_corr_SS=load_2pts("../CORRELATORS/2pts_VKVK_30_30",RE,H_S_2PTS[im_S0],L_S_2PTS[im_spec]).simmetrized(1);
      two_pts_SL_fit(M_VK[im_S0],ZL_VK[im_S0],ZS_VK,VK_corr_SL,VK_corr_SS,tminV,tmaxV,tminV,tmaxV);
      VK_SL_out<<write_constant_fit_plot(effective_mass(VK_corr_SL),M_VK[im_S0],tminV,tmaxV,3*im_S0)<<endl;
      VK_SS_out<<write_constant_fit_plot(effective_mass(VK_corr_SS),M_VK[im_S0],tminV,tmaxV,3*im_S0)<<endl;
    }

  for(int im_S0=0;im_S0<10;im_S0++){for(int ijack=0;ijack<=njack;ijack++) table<<M_P5[im_S0][ijack]<<" ";table<<endl;}
  for(int im_S0=0;im_S0<10;im_S0++){for(int ijack=0;ijack<=njack;ijack++) table<<ZL_P5[im_S0][ijack]<<" ";table<<endl;}
  for(int im_S0=0;im_S0<10;im_S0++){for(int ijack=0;ijack<=njack;ijack++) table<<M_VK[im_S0][ijack]<<" ";table<<endl;}
  for(int im_S0=0;im_S0<10;im_S0++){for(int ijack=0;ijack<=njack;ijack++) table<<ZL_VK[im_S0][ijack]<<" ";table<<endl;}
  
  return 0;
}
