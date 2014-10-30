#include "../b3pts_damir/common.cpp"

const int i2_P5P5=0;
//const int i2_VKVK=1;

struct info_3pts_new
{
  int id;
  PARITY pa;
public:
  info_3pts_new(int id,PARITY pa) : id(id),pa(pa) {}
private:
  info_3pts_new();
};

info_3pts_new iS0P5(0,EVN);
info_3pts_new iS0VK(1,ODD);
info_3pts_new iVKP5(2,EVN);
info_3pts_new iV0P5(3,ODD);
info_3pts_new iTKP5(4,ODD);
info_3pts_new iVJVK(5,ODD);
info_3pts_new iVKVJ(6,ODD);
info_3pts_new iP5VK(7,EVN);
info_3pts_new iAKVK(8,EVN);
info_3pts_new iAJVK(9,EVN); //Horr
info_3pts_new iA0VK(10,ODD);
info_3pts_new iAKV0(11,ODD);
info_3pts_new iA0V0(12,ODD); //Horr
info_3pts_new iTJVK(13,ODD); //Horr
info_3pts_new iTKVJ(14,ODD); //Horr
info_3pts_new iBKVK(15,ODD);
info_3pts_new iBJVK(16,ODD);
info_3pts_new iBKVJ(17,ODD);
info_3pts_new iBKV0(18,EVN);

const int njacks=16;

int icombo_2pts(int iel,int im1,int ith2,int im2,int r=0)
{return iel+2*(im2+12*(r+2*(im1+12*ith2)));}

int icombo_3pts(int iel,int ith1,int im1,int ith2,int im2)
{return iel+19*(im1+3*(ith1+9*(im2+10*ith2)));}

/*
jvec load_2pts_VKVK(int sme,int im1,int im2,int ith=4)
{
  jvec a=jvec_load(combine("../DATA/2pts_30_%02d",sme).c_str(),T,njacks,icombo_2pts(i2_VKVK,im1,ith,im2,0));
  jvec b=jvec_load(combine("../DATA/2pts_30_%02d",sme).c_str(),T,njacks,icombo_2pts(i2_VKVK,im1,ith,im2,1));
  jvec c=jvec_load(combine("../DATA/2pts_30_%02d",sme).c_str(),T,njacks,icombo_2pts(i2_VKVK,im1,8-ith,im2,0));
  jvec d=jvec_load(combine("../DATA/2pts_30_%02d",sme).c_str(),T,njacks,icombo_2pts(i2_VKVK,im1,8-ith,im2,1));
  
  return (a+b+c+d).simmetrized(1)/4;
}
*/

jvec load_2pts_P5P5(int sme,int im1,int im2,int ith=4)
{
  jvec a=jvec_load(combine("../DATA/2pts_30_%02d",sme).c_str(),T,njacks,icombo_2pts(i2_P5P5,im1,ith,im2,0));
  jvec b=jvec_load(combine("../DATA/2pts_30_%02d",sme).c_str(),T,njacks,icombo_2pts(i2_P5P5,im1,ith,im2,1));
  jvec c=jvec_load(combine("../DATA/2pts_30_%02d",sme).c_str(),T,njacks,icombo_2pts(i2_P5P5,im1,8-ith,im2,0));
  jvec d=jvec_load(combine("../DATA/2pts_30_%02d",sme).c_str(),T,njacks,icombo_2pts(i2_P5P5,im1,8-ith,im2,1));
  
  return (a+b+c+d).simmetrized(1)/4;
}

jvec load_3pts(info_3pts_new &info,int im_spec,int ith1,int im1,int ith2,int im2)
{
  jvec a=jvec_load(combine("../DATA/3pts_sp%d_30_30",im_spec).c_str(),T,njacks,icombo_3pts(info.id,ith1,im1,ith2,im2));
  int pa_map[2]={1,-1};
  if(tsep==T/2 && info.pa!=UNC)return a.simmetrized(pa_map[info.pa]);
  else return a.subset(0,tsep+1);
}

int main(int narg,char **arg)
{
  debug_load=debug_fit=0;
  
  if(narg!=2) crash("Use: %s path",arg[0]);
  
  read_set_pars("../data_pars");
  
  const char tag_isp[3]="ls";

  jvec Petoille_corr_SL[2][5];
  jvec Petoille_corr_SS[2][5];
  
  jvec H_corr_SL[2][10];
  jvec H_corr_SS[2][10];
  
  //load 2pts for Pi, K, B and Bs
  for(int isp=0;isp<2;isp++)
    {
      for(int ith=4;ith<9;ith++)
	{      
	  Petoille_corr_SL[isp][ith-4]=load_2pts_P5P5( 0,isp,0,ith);
	  Petoille_corr_SS[isp][ith-4]=load_2pts_P5P5(30,isp,0,ith);
	  char tag_product[2][3]={"Pi","K"};
	  stamp(combine("%s/%s_etoille_SL_th%d",arg[1],tag_product[isp],ith-4),Petoille_corr_SL[isp][ith-4]);
	  stamp(combine("%s/%s_etoille_SS_th%d",arg[1],tag_product[isp],ith-4),Petoille_corr_SS[isp][ith-4]);
	}

      for(int ih=0;ih<10;ih++)
	{
	  //load 2pts for H
	  H_corr_SL[isp][ih]=load_2pts_P5P5( 0,isp,2+ih,4);
	  H_corr_SS[isp][ih]=load_2pts_P5P5(30,isp,2+ih,4);
	  stamp(combine("%s/B%c_m%d_SL",arg[1],tag_isp[isp],ih),H_corr_SL[isp][ih]);
	  stamp(combine("%s/B%c_m%d_SS",arg[1],tag_isp[isp],ih),H_corr_SS[isp][ih]);
	}
    }
  
  //int iS0_2p_l[3]={0,1,1}; //B->Pi, B->K, Bs->K
  int iS0_3p_l[3]={0,1,0}; //B->Pi, B->K, Bs->K
  int isp_l[3]   ={0,0,1}; //B->Pi, B->K, Bs->K
  const char tag_3pts[3][3]={"Pi","K","K"};
  for(int icombo=0;icombo<3;icombo++)
    for(int ih=0;ih<10;ih++)
      for(int ith=4;ith<9;ith++)
	{
	  //int iS0_2p=iS0_2p_l[icombo];
	  int iS0_3p=iS0_3p_l[icombo];
	  int isp=isp_l[icombo];
	  
	  //jvec dt_nu=Petoille_corr_SL[iS0_2p][ith-4].subset(0,tsep+1)*H_corr_SL[isp][ih].subset(0,tsep+1).inverted();
	  
	  jvec HD_V0=(load_3pts(iV0P5,isp,ith,iS0_3p,4,ih)+load_3pts(iV0P5,isp,8-ith,iS0_3p,4,ih))/2;
	  stamp(combine("%s/B%c_m%d_V0_%s_th%d",arg[1],tag_isp[isp],ih,tag_3pts[icombo],ith-4),HD_V0);
	  //(HD_V0/dt_nu).print_to_file(combine("%s/B%c_m%d_V0_%s_th%d.xmg",arg[1],tag_isp[isp],ih,tag_3pts[icombo],ith-4).c_str());

	  jvec HD_VK=(load_3pts(iVKP5,isp,ith,iS0_3p,4,ih)-load_3pts(iVKP5,isp,8-ith,iS0_3p,4,ih))/2;
	  stamp(combine("%s/B%c_m%d_VK_%s_th%d",arg[1],tag_isp[isp],ih,tag_3pts[icombo],ith-4),HD_VK);
	  //(HD_VK/dt_nu).print_to_file(combine("%s/B%c_m%d_VK_%s_th%d.xmg",arg[1],tag_isp[isp],ih,tag_3pts[icombo],ith-4).c_str());
	  
	  jvec HD_TK=(load_3pts(iTKP5,isp,ith,iS0_3p,4,ih)-load_3pts(iTKP5,isp,8-ith,iS0_3p,4,ih))/2;
	  stamp(combine("%s/B%c_m%d_TK_%s_th%d",arg[1],tag_isp[isp],ih,tag_3pts[icombo],ith-4),HD_TK);
	  //(HD_TK/dt_nu).print_to_file(combine("%s/B%c_m%d_TK_%s_th%d.xmg",arg[1],tag_isp[isp],ih,tag_3pts[icombo],ith-4).c_str());
	}
  
  return 0;
}
