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
  
  //prepare the list of mass
  int nh=nmass[2]-nlights[2];
  bvec mh(nh,nboot,njack);
  for(int ih=0;ih<nh;ih++)
    {
      mh[ih]=mass[1][ih+nlights[1]]/lat[1]/Zp[1];
      cout<<ih<<" "<<mh[ih]<<endl;
    }
  
  //compute phi
  bvec M[nens],f[nens],phi[nens];
  
  if(string(obs_name)==string("P5P5"))
    for(int iens=0;iens<nens;iens++)
    {
      int b=ibeta[iens];
      M[iens]=aM[iens]/lat[b];
      f[iens]=sqrt(Z[iens])/(aM[iens]*sinh(aM[iens]))/lat[b];
      phi[iens]=f[iens];
      for(int il=0;il<nlights[iens];il++)
	for(int ih=0;ih<nmass[iens];ih++)
	  {
	    int ic=il*nmass[iens]+ih;
	    f[iens].data[ic]*=mass[iens][il]+mass[iens][ih];
	    phi[iens].data[ic]=f[iens][ic]*sqrt(M[iens][ic]);
	    
	    cout<<iens<<" ens, ic="<<ic<<", il="<<mass[iens][il]<<", ih="<<mass[iens][ih]<<", Z="<<Z[iens][ic]<<", M="<<M[iens][ic]<<endl;
	  }
    }

  if(string(obs_name)==string("S0S0"))
    for(int iens=0;iens<nens;iens++)
      {
	int b=ibeta[iens];
	M[iens]=aM[iens]/lat[b];
	f[iens]=(1/Zp_fr_Zs[b])*sqrt(Z[iens])/(sinh(aM[iens])*aM[iens])/lat[b];
	phi[iens]=f[iens];
	for(int il=0;il<nlights[iens];il++)
	  for(int ih=0;ih<nmass[iens];ih++)
	    {
	      int ic=il*nmass[iens]+ih;
	      f[iens].data[ic]*=-mass[iens][il]+mass[iens][ih];
	      phi[iens].data[ic]=f[iens][ic]*sqrt(M[iens][ic]);
	      
	      cout<<iens<<" ens, ic="<<ic<<", il="<<mass[iens][il]<<", ih="<<mass[iens][ih]<<", Z="<<Z[iens][ic]<<", M="<<M[iens][ic]<<endl;
	    }
      }
  
  bvec intphi[nh];
  //prepare Mhl
  for(int ih=0;ih<nh;ih++)
    {
      intphi[ih]=bvec(nens,nboot,njack);
      for(int iens=0;iens<nens;iens++)
	{
	  int ib=ibeta[iens];
	  int nl=nlights[iens];
	  int nm=nmass[iens];
	  int ni=nm-nl;
	  
	  //prepare input
	  double mi[ni];
	  bvec phii(ni,nboot,njack);
	  for(int im=0;im<ni;im++)
	    {
	      mi[im]=mass[iens][im+nl];
	      phii[im]=phi[iens][im+nl];
	    }
	  
	  //interpolate
	  boot x=mh[ih]*lat[ib]*Zp[ib];
	  intphi[ih][iens]=interpolate_single(phii,mi,x);
	  
	  cout<<ih<<" iens="<<iens<<" x="<<x.med()<<" y="<<intphi[ih][iens]<<endl;
	}
      if(ih==0) intphi[ih].write_to_binfile("intphi");
      else      intphi[ih].append_to_binfile("intphi");
    }
  
  return 0;
}
