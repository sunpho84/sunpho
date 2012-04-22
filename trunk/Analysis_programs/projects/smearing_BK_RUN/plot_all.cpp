#include <include.h>
#include <iostream>
#include <sstream>

#include "../nf2/common.cpp"
#include "../nf2/interpolate/interpolate_lib.cpp"

using namespace std;

int main()
{
  init_latpars();

  //load ensemble list path and data path
  FILE *an_input_file=open_file("analysis_pars","r");
  char ens_list_path[1024],base_MZ_path[1024],obs_name[1024],meson_name[1024];
  int mode;
  read_formatted_from_file_expecting(ens_list_path,an_input_file,"%s","ens_list_path");
  read_formatted_from_file_expecting(base_MZ_path,an_input_file,"%s","base_MZ_path");
  read_formatted_from_file_expecting((char*)&mode,an_input_file,"%s","mode");
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
  load_all_ensembles_MZ(aM,Z,nens,T,ibeta,nlights,nmass,base_MZ_path,obs_name,ens_name,base_corrs_path,1);
  
  //compute M
  bvec M[nens];
  for(int iens=0;iens<nens;iens++)
    {
      int b=ibeta[iens];
      M[iens]=aM[iens]/lat[b];
    }
  
  ofstream out("all_heavy.xmg");
  out<<"@type xydy"<<endl;
  for(int iens=0;iens<nens;iens++)
    {
      for(int ic=nlights[iens];ic<nmass[iens];ic++)
	out<<(mass[iens][ic]/lat[ibeta[iens]]/Zp[ibeta[iens]]).med()<<" "<<M[iens][ic]<<endl;
      out<<"&"<<endl;
    }
  
  return 0;
}
