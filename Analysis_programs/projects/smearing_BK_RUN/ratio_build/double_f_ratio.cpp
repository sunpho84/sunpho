#include "../../nf2/common.cpp"
#include "../../nf2/interpolate/interpolate_lib.cpp"

int ico(int il,int ih,int nm)
{return il*nm+ih;}

int main()
{
  //load ensemble list path and data path
  FILE *an_input_file=open_file("analysis_pars","r");
  int mode;
  char ens_list_path[1024],base_MZ_path[1024],obs_name[1024],meson_name[1024];
  read_formatted_from_file_expecting(ens_list_path,an_input_file,"%s","ens_list_path");
  read_formatted_from_file_expecting(base_MZ_path,an_input_file,"%s","base_MZ_path");
  read_formatted_from_file_expecting((char*)&mode,an_input_file,"%d","mode");
  read_formatted_from_file_expecting(obs_name,an_input_file,"%s","obs_name");
  read_formatted_from_file_expecting(meson_name,an_input_file,"%s","meson_name");
  fclose(an_input_file);
  
  //load ensembles list and parameters
  char **base_corrs_path,**ens_name;
  int nens,*T,*ibeta,*iml_un,*nlights,*nmass;
  double **mass;
  load_ensembles_list(base_corrs_path,ens_name,nens,T,ibeta,nmass,mass,iml_un,nlights,ens_list_path);
  
  //load all ensembles data
  bvec *aM,*Z;
  load_all_ensembles_MZ(aM,Z,nens,T,ibeta,nlights,nmass,base_MZ_path,obs_name,ens_name,base_corrs_path,mode);
  init_latpars();
  
  //prepare the list of mass
  int nh=nmass[2]-nlights[2];
  bvec mh(nh,nboot,njack);
  for(int ih=0;ih<nh;ih++)
    {
      mh[ih]=mass[1][ih+nlights[1]]/lat[1]/Zp[1];
      //cout<<ih<<" "<<mh[ih]<<endl;
    }
  
  bvec M[nens],f[nens],ratio(nens,nboot,njack);

  cout<<mass[2][0]<<endl;
  cout<<mass[2][2]<<endl;
  cout<<mass[2][nlights[2]+1]<<endl;
  
  for(int iens=0;iens<nens;iens++)
    {
      int b=ibeta[iens];
      M[iens]=aM[iens]/lat[b];
      f[iens]=sqrt(Z[iens])/(aM[iens]*sinh(aM[iens]))/lat[b];
      for(int il=0;il<nlights[iens];il++)
	for(int ih=0;ih<nmass[iens];ih++)
	  {
	    int ic=il*nmass[iens]+ih;
	    f[iens].data[ic]*=mass[iens][il]+mass[iens][ih];
	  }
      boot fPi=f[iens].data[ico(0,0,nmass[iens])];
      boot fK=f[iens].data[ico(2,0,nmass[iens])];
      f[iens]*=sqrt(M[iens]);
      
      boot fD=f[iens].data[ico(0,nlights[iens]+1,nmass[iens])]; 
      boot fDs=f[iens].data[ico(2,nlights[iens]+1,nmass[iens])]; 
      
      boot ratio_D=fDs/fD/sqrt(19685.0/18696);
      boot ratio_L=fK/fPi;
      
      //ratio[iens]=ratio_D/ratio_L;
      ratio[iens]=ratio_L;
      
      cout<<iens<<" "<<ratio[iens]<<endl;
    } 
  
  ratio.write_to_binfile("ratio");
  
  return 0;
}
