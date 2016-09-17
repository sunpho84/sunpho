#include "common.cpp"

int icombo_2pts(REIM ri,int r,S_2pts_mass im1,S_2pts_mass im2)
{return ri+2*(im1+nm*(r+2*im2));}

int icombo_3pts(REIM ri,S0_3pts_mass im_S0,S1_3pts_mass im_S1,int ith_S0)
{return ri+2*(im_S0+nm_S0*(ith_S0+nth_S0*im_S1));}

jvec load_2pts(const char *path,REIM ri,S_2pts_mass im1,S_2pts_mass im2)
{return (jvec_load(path,T,njack,icombo_2pts(ri,0,im1,im2))+jvec_load(path,T,njack,icombo_2pts(ri,1,im1,im2)))/2;}

jvec load_3pts(info_3pts &info,combo_3pts combo)
{
  jvec a=-jvec_load(combine("../CORRELATORS/3pts_sp%d_%s_30_30",combo.im_spec,info.name).c_str(),T,njack,icombo_3pts(info.ri,combo.im_S0,combo.im_S1,combo.ith_S0));
  
  //correlators do not need to be changed, even if we are looking at mirror process
  //we checked performed explicitely another computation  
  if(info.pa!=UNC && tsep==TH) return a.simmetrized(map_pa[info.pa]);
  else return a.subset(0,tsep+1);
}

int main(int narg,char **arg)
{
  if(narg!=2) crash("Use: %s path",arg[0]);
  read_set_pars("../data_pars");
  
  const char tag_isp[3]="ls";
  for(int isp=0;isp<2;isp++)
    {
      //load 2pts for D
      jvec Detoille_corr_SL=load_2pts("../CORRELATORS/2pts_VKVK_30_00",RE,H0_S_2PTS,L_S_2PTS[isp]).simmetrized(1);
      jvec Detoille_corr_SS=load_2pts("../CORRELATORS/2pts_VKVK_30_30",RE,H0_S_2PTS,L_S_2PTS[isp]).simmetrized(1);
      stamp(combine("%s/D%c_etoille_SL",arg[1],tag_isp[isp]),Detoille_corr_SL);
      stamp(combine("%s/D%c_etoille_SS",arg[1],tag_isp[isp]),Detoille_corr_SS);
      
      for(int ih=0;ih<10;ih++)
	{
	  //load 2pts for B
	  jvec B_corr_SL=load_2pts("../CORRELATORS/2pts_P5P5_30_00",RE,H_S_2PTS[ih],L_S_2PTS[isp]).simmetrized(1);
	  jvec B_corr_SS=load_2pts("../CORRELATORS/2pts_P5P5_30_30",RE,H_S_2PTS[ih],L_S_2PTS[isp]).simmetrized(1);
	  stamp(combine("%s/B%d%c_SL",arg[1],ih,tag_isp[isp]),B_corr_SL);
	  stamp(combine("%s/B%d%c_SS",arg[1],ih,tag_isp[isp]),B_corr_SS);
	  
	  //load 3pts
	  combo_3pts combo(SPEC[isp],CHARM_S0_3PTS,H_S1_3PTS[ih],0);
	  jvec BD_AK=load_3pts(AKVK,combo);
	  stamp(combine("%s/B%d%c_Ax_D%c",arg[1],ih,tag_isp[isp],tag_isp[isp]),BD_AK);
	}
    }
  
  return 0;
}
