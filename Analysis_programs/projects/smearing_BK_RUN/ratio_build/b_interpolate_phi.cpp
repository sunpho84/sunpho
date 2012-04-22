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
  
  for(int iens=0;iens<nens;iens++)
    {
      int ib=ibeta[iens];
      int nl=nlights[iens];
      int nm=nmass[iens];
      int ni=nm-nl;
      
      ofstream out_rat(combine("rat_%d.out",iens).c_str());
      ofstream out_phi(combine("phi_%d.out",iens).c_str());
      
      out_rat<<"@type xydy"<<endl;
      out_phi<<"@type xydy"<<endl;
      out_phi<<1/mass[iens][nl]<<" "<<phi[iens][nl]<<endl;

      double x[nh];
      bvec ratio(nh,nboot,njack);
      for(int ih1=1;ih1<nh;ih1++)
	{
	  int ih0=ih1-1;
	  
	  x[ih0]=1/mass[iens][ih1+nl];
	  ratio[ih0]=phi[iens][ih1+nl]/phi[iens][ih0+nl];
	  
	  out_phi<<1/mass[iens][ih1+nl]<<" "<<phi[iens][ih1+nl]<<endl;
	  out_rat<<x[ih0]<<" "<<ratio[ih0]<<endl;
	}
      x[nh-1]=0;
      ratio[nh-1]=pow(ratio[0],1/100.0);
      
      bvec d=poly_fit(x,ratio,2);
      
      out_rat<<"&"<<endl;
      for(double t=0;t<x[0];t+=0.1)
	out_rat<<t<<" "<<d[0]+d[1]*t+d[2]*t*t<<endl;
      
      out_phi<<"&"<<endl;
      boot run=phi[iens][nl];
      out_phi<<1/mass[iens][nl]<<" "<<run<<endl;
      for(int ih1=1;ih1<nh;ih1++)
	{
	  double p=1/mass[iens][ih1+nl];
	  run*=(d[0]+d[1]*p+d[2]*p*p);
	  out_phi<<p<<" "<<run<<endl;
	}
      
      out_rat.close();
      out_phi.close();
      
      cout<<mass[iens][iml_un[iens]]/lat[ib].med()/Zp[ib].med()<<" "<<" "<<run/pow(lat[ib],1.5)<<endl;
    }
  
  return 0;
}
