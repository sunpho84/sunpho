#include "common.cpp"
#include "interpolate/interpolate_lib.cpp"

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
  
  cout.precision(7);
  
  char beta[4][10]={"3.80","3.90","4.05","4.20"};
      cout<<" Beta\tL\tmsea\tm1\tm2\tM\t\tZ2"<<endl;
  for(int iens=0;iens<nens;iens++)
    {
      for(int im1=max(nlights[iens]-3,iml_un[iens]+1);im1<nlights[iens];im1++)
	cout<<beta[ibeta[iens]]<<"\t"<<T[iens]/2<<"\t"<<mass[iens][iml_un[iens]]<<"\t"<<mass[iens][iml_un[iens]]<<"\t"<<mass[iens][im1]<<"\t"<<smart_print(aM[iens][icombo(iml_un[iens],im1,nmass[iens],nlights[iens],mode)])<<"\t"<<smart_print(Z[iens][icombo(iml_un[iens],im1,nmass[iens],nlights[iens],mode)])<<endl;
    }
  
  return 0;
}
