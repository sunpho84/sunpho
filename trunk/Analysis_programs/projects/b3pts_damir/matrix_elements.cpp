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

void read_input(const char *path)
{
  FILE *input_file=open_file(path,"r");
  
  read_formatted_from_file_expecting((char*)&im_spec,input_file,"%d","im_spec");
  read_formatted_from_file_expecting((char*)&im_S0,input_file,"%d","im_S0");
  read_formatted_from_file_expecting((char*)&im_S1,input_file,"%d","im_S1");
  
  cout<<"Considering "<<ori_meson[im_spec][im_S0]<<" -> "<<prod_meson[im_spec][min(im_S1,1)]<<" transition"<<endl;
  
  read_formatted_from_file_expecting((char*)&tmin_K,input_file,"%d","tint_K");
  read_formatted_from_file((char*)&tmax_K,input_file,"%d","tint_K");

  read_formatted_from_file_expecting((char*)&tmin_H,input_file,"%d","tint_H");
  read_formatted_from_file((char*)&tmax_H,input_file,"%d","tint_K");
  
  read_formatted_from_file_expecting((char*)&tmin_V0,input_file,"%d","tint_V0");
  read_formatted_from_file((char*)&tmax_V0,input_file,"%d","tint_V0");
  
  read_formatted_from_file_expecting((char*)&tmin_VK,input_file,"%d","tint_VK");
  read_formatted_from_file((char*)&tmax_VK,input_file,"%d","tint_VK");
  
  read_formatted_from_file_expecting((char*)&tmin_TK,input_file,"%d","tint_TK");
  read_formatted_from_file((char*)&tmax_TK,input_file,"%d","tint_TK");
  
  fclose(input_file);
}

int main()
{
  read_set_pars("../data_pars");
  read_input("analysis_pars");
  
  
  ////////////////////////// P5 Kaon /////////////////////
  
  jack M_K_P5,ZL_K_P5,ZS_K_P5;
  
  {
    //load P5 Kaon
    jvec K_P5_corr_SL=load_2pts("../CORRELATORS/2pts_P5P5_30_00",RE,L_S_2PTS[im_S0],L_S_2PTS[im_spec]).simmetrized(1);
    jvec K_P5_corr_SS=load_2pts("../CORRELATORS/2pts_P5P5_30_30",RE,L_S_2PTS[im_S0],L_S_2PTS[im_spec]).simmetrized(1);
    
    //fit P5 Kaon
    two_pts_SL_fit(M_K_P5,ZL_K_P5,ZS_K_P5,K_P5_corr_SL,K_P5_corr_SS,tmin_K,tmax_K,tmin_K,tmax_K,"masses_Z/plots/M_K_P5_fit.xmg");
    
    //write the mass and Z
    {
      ofstream raw("masses_Z/plots/raw_K.xmg");
      raw<<"@type xydy"<<endl;
      raw<<K_P5_corr_SL<<"&"<<endl<<K_P5_corr_SS<<endl;
    }
    M_K_P5.write_to_binfile("masses_Z/M_K_P5");
    ZL_K_P5.write_to_binfile("masses_Z/Z_K_P5");

    cout<<"MK "<<M_K_P5<<endl;
    cout<<"ZL2K "<<sqr(ZL_K_P5)<<", ZLK "<<ZL_K_P5<<endl;
    cout<<"ZS2K "<<sqr(ZS_K_P5)<<", ZSK "<<ZS_K_P5<<endl;
  }
    
  /////////////////////////// VK Kaon /////////////////////
  
  /*
    
  jack M_K_VK,ZL_K_VK,ZS_K_VK;
  
  {
    //load VK Kaon
    jvec K_VK_corr_SL=load_2pts("../CORRELATORS/2pts_VKVK_30_00",RE,L_S_2PTS[im_S0],L_S_2PTS[im_spec]).simmetrized(1);
    jvec K_VK_corr_SS=load_2pts("../CORRELATORS/2pts_VKVK_30_30",RE,L_S_2PTS[im_S0],L_S_2PTS[im_spec]).simmetrized(1);
    
    //fit VK Kaon
    two_pts_SL_fit(M_K_VK,ZL_K_VK,ZS_K_VK,K_VK_corr_SL,K_VK_corr_SS,tmin_K,tmax_K,tmin_K,tmax_K,"masses_Z/plots/M_K_VK_fit.xmg");
    
    //write the mass and Z
    M_K_VK.write_to_binfile("masses_Z/M_K_VK");
    ZL_K_VK.write_to_binfile("masses_Z/Z_K_VK");
  }
  */
  
  ////////////////////////// P5 Heavy //////////////////////
  
  //determine M and Z for H meson
  jack M_H_P5,ZL_H_P5,ZS_H_P5;
  
  {
    //load H
    jvec H_P5_corr_SL=load_2pts("../CORRELATORS/2pts_P5P5_30_00",RE,H_S_2PTS[im_S1],L_S_2PTS[im_spec]).simmetrized(1);
    jvec H_P5_corr_SS=load_2pts("../CORRELATORS/2pts_P5P5_30_30",RE,H_S_2PTS[im_S1],L_S_2PTS[im_spec]).simmetrized(1);
    
    //fit H
    two_pts_SL_fit(M_H_P5,ZL_H_P5,ZS_H_P5,H_P5_corr_SL,H_P5_corr_SS,tmin_H,tmax_H,tmin_H,tmax_H,
		   "masses_Z/plots/M_H_P5_fit_SL.xmg","masses_Z/plots/M_H_P5_fit_SS.xmg");
    
    //write the mass and Z
    {
      ofstream raw("masses_Z/plots/raw_H.xmg");
      raw<<"@type xydy"<<endl;
      raw<<H_P5_corr_SL<<"&"<<endl<<H_P5_corr_SS<<endl;
    }
    M_H_P5.write_to_binfile("masses_Z/M_H_P5");
    ZL_H_P5.write_to_binfile("masses_Z/Z_H_P5");
    
    cout<<"MH "<<M_H_P5<<endl;
    cout<<"ZL2H "<<sqr(ZL_H_P5)<<endl;
    cout<<"ZS2H "<<sqr(ZS_H_P5)<<endl;
  }
  
  //////////////////////////////////// P5 {AK,BK,P5,VJ,TJ} VK ////////////////////////////////
  
  //compute matrix elements
  for(int i_ME=0;i_ME<n_ME;i_ME++)
    {
      info_3pts info=info_ME[i_ME];
      
      //create the name with "_"
      char spaced_name[12];
      sprintf(spaced_name,"P5_%c%c_%s",*info.name,*(info.name+1),info.name+2);
      
      //open the output
      ofstream out_ME(combine("matrix_elements/plots/%s.xmg",spaced_name).c_str());
      
      //stored TK for null theta
      jvec null_TKP5,null_VKP5;
      
      //loop over theta
      jvec ME(nth_S0,njack);
      for(int ith=0;ith<nth_S0;ith++)
	{
	  //find the combo
	  combo_3pts combo(SPEC[im_spec],L_S0_3PTS[im_S0],H_S1_3PTS[im_S1],ith);
	  
	  //define the temporal dependance
	  //jack E_K_PV=latt_en((info.name[2]=='P')?M_K_P5:M_K_VK,th_S0[ith]);
	  jack E_K_PV=latt_en(M_K_P5,th_S0[ith]);
	  jvec dT(tsep+1,njack);
	  for(int t=0;t<=tsep;t++)
	    dT[t]=(ZS_K_P5*ZS_H_P5)/(2*E_K_PV*2*M_H_P5)*exp(-E_K_PV*t)*exp(-M_H_P5*(tsep-t));
	  
	  //load the corr and divide it by dT
	  jvec H_P5_ME_K_PV_corr=load_3pts(info,combo);
	  
	  //store "0" for TK
	  if(i_ME==3)
	    {
	      if(ith==0) null_TKP5=H_P5_ME_K_PV_corr;
	      else H_P5_ME_K_PV_corr-=null_TKP5;
	    }
	  
	  //print raw for V0
	  if(i_ME==2) H_P5_ME_K_PV_corr.print_to_file(combine("matrix_elements/plots/raw_V0%d.xmg",ith).c_str());

	  //same for VK
	  if(i_ME==1)
	    {
	      H_P5_ME_K_PV_corr.print_to_file(combine("matrix_elements/plots/raw_VK%d.xmg",ith).c_str());
	      
	      if(ith==0) null_VKP5=H_P5_ME_K_PV_corr;
	      else H_P5_ME_K_PV_corr-=null_VKP5;
	    }
	  
	  H_P5_ME_K_PV_corr/=dT;

	  int tmin_ME=tmin_VK,tmax_ME=tmax_VK;
	  if(i_ME==1){tmin_ME=tmin_VK;tmax_ME=tmax_VK;}
	  if(i_ME==2){tmin_ME=tmin_V0;tmax_ME=tmax_V0;}
	  if(i_ME==3){tmin_ME=tmin_TK;tmax_ME=tmax_TK;}
	  ME[ith]=constant_fit(H_P5_ME_K_PV_corr,tmin_ME,tmax_ME);

	  //output
	  if(ith<5) out_ME<<"@type xydy\n"<<H_P5_ME_K_PV_corr<<"&\n@type xy\n"<<write_constant_with_error(ME[ith],tmin_ME,tmax_ME)<<"&\n";
	}
      
      //write 
      ME.write_to_binfile(combine("matrix_elements/%s",spaced_name).c_str());
    }
  
  //divita ratios
  {
    ofstream out_1("matrix_elements/plots/divita_ratio1.xmg");
    ofstream out_2("matrix_elements/plots/divita_ratio2.xmg");
    ofstream out_3("matrix_elements/plots/divita_ratio3.xmg");
    out_1<<"@type xydy"<<endl;
    out_2<<"@type xydy"<<endl;
    out_3<<"@type xydy"<<endl;
    jvec rat_1(nth_S0,njack);
    jvec rat_2(nth_S0,njack);
    jvec rat_3(nth_S0,njack);
    for(int ith=0;ith<nth_S0;ith++)
	{
	  //find the combo
	  int iVK=1,iV0=2,iTK=3;
	  combo_3pts combo_th(SPEC[im_spec],L_S0_3PTS[im_S0],H_S1_3PTS[im_S1],ith);
	  combo_3pts combo_0(SPEC[im_spec],L_S0_3PTS[im_S0],H_S1_3PTS[im_S1],0);
	  
	  //compute deltaE
	  jack E_K_P5=latt_en(M_K_P5,th_S0[ith]);
	  jack D_K_P5=E_K_P5-M_K_P5;

	  //load and take note
	  jvec rat_1_corr=load_3pts(info_ME[iV0],combo_th)/load_3pts(info_ME[iV0],combo_0);
	  for(int t=0;t<=tsep;t++) rat_1_corr[t]*=exp(D_K_P5*t)*E_K_P5/M_K_P5;
	  rat_1[ith]=constant_fit(rat_1_corr,tmin_VK,tmax_VK);
	  
	  //load and take note
	  jvec rat_2_corr=load_3pts(info_ME[iVK],combo_th)/load_3pts(info_ME[iV0],combo_th)-
	    load_3pts(info_ME[iVK],combo_0)/load_3pts(info_ME[iV0],combo_0);
	  rat_2[ith]=constant_fit(rat_2_corr,tmin_VK,tmax_VK);

	  jvec rat_3_corr=(load_3pts(info_ME[iTK],combo_th)-load_3pts(info_ME[iTK],combo_0))/
	    (load_3pts(info_ME[iVK],combo_th)-load_3pts(info_ME[iVK],combo_0));
	  rat_3[ith]=constant_fit(rat_3_corr,tmin_VK,tmax_VK);


          out_1<<"@type xydy\n"<<rat_1_corr<<"&\n@type xy\n"<<write_constant_with_error(rat_1[ith],tmin_VK,tmax_VK)<<"&\n";
          out_2<<"@type xydy\n"<<rat_2_corr<<"&\n@type xy\n"<<write_constant_with_error(rat_2[ith],tmin_VK,tmax_VK)<<"&\n";
          if(ith) out_3<<"@type xydy\n"<<rat_3_corr<<"&\n@type xy\n"<<write_constant_with_error(rat_3[ith],tmin_VK,tmax_VK)<<"&\n";
	}
    
    rat_1.write_to_binfile("matrix_elements/divita_ratio_1");
    rat_2.write_to_binfile("matrix_elements/divita_ratio_2");
    rat_3.write_to_binfile("matrix_elements/divita_ratio_3");
  }
  
  return 0;
}
