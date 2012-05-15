#include "../../nf2/common.cpp"
#include "../../nf2/interpolate/interpolate_lib.cpp"

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
  
  //output
  bvec Mout(nref_hmass*nens,nboot,njack);
  
  for(int iens=0;iens<nens;iens++)
    {
      ///here we will store stranged interpolated
      bvec Mt(nmass[iens],nboot,njack);
      //now interpolate to the strange each heavy
      for(int im=0;im<nmass[iens];im++)
	{
	  //sort out light at fixed heavy
	  bvec ytemp(nlights[iens],nboot,njack);
	  for(int iml=0;iml<nlights[iens];iml++)
	    ytemp[iml]=aM[iens][iml*nmass[iens]+im]/lat[ibeta[iens]];
	  
	  //interpolate at fixed heavy
	  Mt[im]=interpolate_single(mass[iens],ytemp,ams_phys[ibeta[iens]],combine("interpolating_strange_ens_%02d_heavy_%02d.xmg",iens,im).c_str());
	}
      
      //produce heavy ref bare mass
      bvec mass_out(nref_hmass,nboot,njack);
      for(int iref_hmass=0;iref_hmass<nref_hmass;iref_hmass++)
        mass_out[iref_hmass]=ref_hmass[iref_hmass]*lat[ibeta[iens]]*Zp[ibeta[iens]];
      
      //interpolate
      bvec temp_Mout=interpolate_multi(mass[iens],Mt,mass_out,combine("interpolating_heavy_ens_%02d.xmg",iens).c_str());
      
      for(int iref_hmass=0;iref_hmass<nref_hmass;iref_hmass++)
        {
          Mout[iens+iref_hmass*nens]=temp_Mout[iref_hmass];
          cout<<iref_hmass<<" iens="<<iens<<" x="<<mass_out[iref_hmass]<<" y="<<temp_Mout[iref_hmass]<<endl;
        }
    }
  
  Mout.write_to_binfile("intM");

  return 0;
}
