#include "include.h"

int T=48,TH=T/2,tsep=24;
int nmass=10;
int njacks=16;

int icombo_3pts(int itheta,int imass,int ispec,int ri)
{return ri+2*(ispec+2*(imass+nmass*itheta));}

int icombo_2pts(int itheta,int imass,int ispec)
{return /*always real*/0+2*(/*pion*/1+ispec+2*(imass+nmass*itheta));}

int main()
{
  //load Pi corr
  jvec Pi_sl=jvec_load("DATA/2pts_sl_corr_P5P5",T,njacks,0).simmetrized(1);
  jvec Pi_ss=jvec_load("DATA/2pts_ss_corr_P5P5",T,njacks,0).simmetrized(1);
  
  //fit Pi mass and z
  jack M_Pi(njacks),ZL_Pi(njacks),ZS_Pi(njacks);
  two_pts_SL_fit(M_Pi,ZL_Pi,ZS_Pi,Pi_sl,Pi_ss,10,TH,10,TH,"PLOTS/Pi.xmg");

  int ih=4;
  int iD_st_lspec=icombo_2pts(1,ih,0);
  int iD_M_cspec=icombo_2pts(0,ih,1);
  int iD_P_cspec=icombo_2pts(2,ih,1);
  
  //---------
  
  //load D corr
  jvec D_sl=jvec_load("DATA/2pts_sl_corr_P5P5",T,njacks,iD_st_lspec).simmetrized(1);
  jvec D_ss=jvec_load("DATA/2pts_ss_corr_P5P5",T,njacks,iD_st_lspec).simmetrized(1);
  
  //fit D mass and z
  jack M_D(njacks),ZL_D(njacks),ZS_D(njacks);
  two_pts_SL_fit(M_D,ZL_D,ZS_D,D_sl,D_ss,10,TH,10,TH,"PLOTS/D.xmg");

  //load D corr in motion
  jvec D_cspec_sl=(jvec_load("DATA/2pts_sl_corr_P5P5",T,njacks,iD_M_cspec).simmetrized(1)+
		   jvec_load("DATA/2pts_sl_corr_P5P5",T,njacks,iD_P_cspec).simmetrized(1))/2;
  jvec D_cspec_ss=(jvec_load("DATA/2pts_ss_corr_P5P5",T,njacks,iD_M_cspec).simmetrized(1)+
		   jvec_load("DATA/2pts_ss_corr_P5P5",T,njacks,iD_P_cspec).simmetrized(1))/2;
  
  //fit D energy and z
  jack E_D(njacks),ZL_D_cspec(njacks),ZS_D_cspec(njacks);
  two_pts_SL_fit(E_D,ZL_D_cspec,ZS_D_cspec,D_cspec_sl,D_cspec_ss,10,TH,10,TH,"PLOTS/D_cspec.xmg");

  cout<<"Variation of ZD: "<<smart_print(ZS_D_cspec/ZS_D*100-100)<<"%"<<endl;
  
  //---------
  
  //load D* corr
  jvec Dst_sl=jvec_load("DATA/2pts_sl_corr_VKVK",T,njacks,iD_st_lspec).simmetrized(1);
  jvec Dst_ss=jvec_load("DATA/2pts_ss_corr_VKVK",T,njacks,iD_st_lspec).simmetrized(1);
  
  //fit D* mass and z
  jack M_Dst(njacks),ZL_Dst(njacks),ZS_Dst(njacks);
  two_pts_SL_fit(M_Dst,ZL_Dst,ZS_Dst,Dst_sl,Dst_ss,8,20,8,20,"PLOTS/Dst.xmg");
  
  //---------
  
  int iDDst_rest=icombo_3pts(1/*theta*/,ih/*mass*/,0/*spec*/,0/*rim*/);
  int iDDst_lspec_M=icombo_3pts(0/*theta*/,ih/*mass*/,0/*spec*/,0/*rim*/);
  int iDDst_lspec_P=icombo_3pts(2/*theta*/,ih/*mass*/,0/*spec*/,0/*rim*/);
  int iDDst_cspec_M=icombo_3pts(0/*theta*/,ih/*mass*/,1/*spec*/,0/*rim*/);
  int iDDst_cspec_P=icombo_3pts(2/*theta*/,ih/*mass*/,1/*spec*/,0/*rim*/);
  
  jvec AKVK_rest_corr=jvec_load("DATA/3pts_corr_AKVK",T,njacks,iDDst_rest);
  jvec AKVK_lspec_M_corr=jvec_load("DATA/3pts_corr_AKVK",T,njacks,iDDst_lspec_M);
  jvec AKVK_lspec_P_corr=jvec_load("DATA/3pts_corr_AKVK",T,njacks,iDDst_lspec_P);
  jvec AKVK_cspec_M_corr=jvec_load("DATA/3pts_corr_AKVK",T,njacks,iDDst_cspec_M);
  jvec AKVK_cspec_P_corr=jvec_load("DATA/3pts_corr_AKVK",T,njacks,iDDst_cspec_P);
  jvec AKVK_lspec_corr=(AKVK_lspec_M_corr+AKVK_lspec_P_corr)/2;
  jvec AKVK_cspec_corr=(AKVK_cspec_M_corr+AKVK_cspec_P_corr)/2;
  AKVK_rest_corr.print_to_file("PLOTS/AKVK_rest.xmg");  
  AKVK_lspec_M_corr.print_to_file("PLOTS/AKVK_lspec_M.xmg");  
  AKVK_lspec_P_corr.print_to_file("PLOTS/AKVK_lspec_P.xmg");  
  AKVK_cspec_M_corr.print_to_file("PLOTS/AKVK_cspec_M.xmg");  
  AKVK_cspec_P_corr.print_to_file("PLOTS/AKVK_cspec_P.xmg");  
  AKVK_lspec_corr.print_to_file("PLOTS/AKVK_lspec.xmg");
  AKVK_cspec_corr.print_to_file("PLOTS/AKVK_cspec.xmg");

  //symmetrize or take only the relevant subset
  if(tsep==TH)
    {
      AKVK_cspec_corr=AKVK_cspec_corr.simmetrized(1);
      AKVK_lspec_corr=AKVK_lspec_corr.simmetrized(1);
    }
  else
    {
      AKVK_cspec_corr=AKVK_cspec_corr.subset(0,tsep+1);
      AKVK_lspec_corr=AKVK_lspec_corr.subset(0,tsep+1);
    }
  
  ofstream AK_eff_mass_rest_out("PLOTS/AK_eff_mass_rest.xmg");
  AK_eff_mass_rest_out<<"@type xydy"<<endl;
  AK_eff_mass_rest_out<<aperiodic_effective_mass(AKVK_rest_corr)<<endl;
  AK_eff_mass_rest_out<<"&"<<endl;
  AK_eff_mass_rest_out<<write_line_with_error(M_Dst-M_D,M_D*0,0,TH,2)<<endl;
  ofstream AK_eff_mass_mot_out("PLOTS/AK_eff_mass_mot.xmg");
  AK_eff_mass_mot_out<<"@type xydy"<<endl;
  AK_eff_mass_mot_out<<aperiodic_effective_mass(AKVK_cspec_corr)<<endl;
  AK_eff_mass_mot_out<<"&"<<endl;
  AK_eff_mass_mot_out<<aperiodic_effective_mass(AKVK_lspec_corr)<<endl;
  AK_eff_mass_mot_out<<"&"<<endl;
  
  (-0.5*D_ss[tsep]/jvec_load("DATA/3pts_corr_V0P5",T,njacks,iDDst_rest)).subset(1,TH-1).print_to_file("PLOTS/Zv.xmg");

  //compute dt
  jvec dt_cspec(tsep+1,njacks);
  for(int t=0;t<=tsep;t++)
    dt_cspec[t]=(ZS_D_cspec*ZS_Dst)*
      exp((-M_Dst*t)+(-E_D*(tsep-t)))/
       (2*E_D*2*M_Dst);

  jack AKVK_cspec_mel=constant_fit(AKVK_cspec_corr/dt_cspec,4,8,"PLOTS/AKVK_ratio_an.xmg");

  jack A1=AKVK_cspec_mel/(M_Dst+M_D);
  jack fpi_gDvDPi_pt1=AKVK_cspec_mel;
  double Za=0.746;
  jack gc_pt1=fpi_gDvDPi_pt1*Za/(2*sqrt(M_D*M_Dst));
  cout<<"gc: "<<smart_print(gc_pt1)<<endl;
  return 0;
}
