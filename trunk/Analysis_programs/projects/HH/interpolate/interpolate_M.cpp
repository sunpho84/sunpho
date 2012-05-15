#include "../HH_common.cpp"

int main()
{
  int frac=0;
  
  //load ensemble list path and data path
  FILE *an_input_file=open_file("analysis_pars","r");
  char ens_list_path[1024],base_M_path[1024],obs_name[1024];
  read_formatted_from_file_expecting(ens_list_path,an_input_file,"%s","ens_list_path");
  read_formatted_from_file_expecting(base_M_path,an_input_file,"%s","base_M_path");
  read_formatted_from_file_expecting(obs_name,an_input_file,"%s","obs_name");
  fclose(an_input_file);
  
  //load ensembles list and parameters
  char **base_corrs_path,**ens_name;
  int nens,*T,*ibeta,*nmass;
  double **mass,*sea_mass;
  load_ensembles_list(base_corrs_path,ens_name,nens,T,ibeta,nmass,mass,sea_mass,ens_list_path);
  
  //load all ensembles data
  bvec *aM;
  if(strlen(obs_name)==4)
    load_all_ensembles_M(aM,nens,T,ibeta,nmass,base_M_path,obs_name,ens_name,base_corrs_path,"M");
  else
    {
      char sep[2]="_";
      char *last;
      char *first=strtok_r(obs_name,sep,&last);
      char *opera=strtok_r(NULL,sep,&last);
      char *secon=strtok_r(NULL,sep,&last);
      cout<<first<<" "<<opera<<" "<<secon<<endl;
      bvec *aM2;
      load_all_ensembles_M(aM,nens,T,ibeta,nmass,base_M_path,first,ens_name,base_corrs_path,"M");
      load_all_ensembles_M(aM2,nens,T,ibeta,nmass,base_M_path,secon,ens_name,base_corrs_path,"M");
      for(int iens=0;iens<nens;iens++)
	{
	  if(string(opera)=="minus")
	    aM[iens]-=aM2[iens];
	  if(string(opera)=="frac")
	    {
	      aM[iens]/=aM2[iens];
	      frac=1;
	    }	
	}
    }
  
  init_latpars();
  
  bvec Mout(nref_hmass*nens,nboot,njack);
  //prepare Mhl
  for(int iens=0;iens<nens;iens++)
    {
      int ib=ibeta[iens];
      
      //prepare input
      bvec mass_out(nref_hmass,nboot,njack);
      bvec Min(nmass[iens],nboot,njack);
      for(int im=0;im<nmass[iens];im++)
	if(!frac)
	  Min[im]=aM[iens][icombo(im,im,nmass[iens])]/lat[ib];
	else
	  Min[im]=aM[iens][icombo(im,im,nmass[iens])]-1;
      for(int iref_hmass=0;iref_hmass<nref_hmass;iref_hmass++)
	mass_out[iref_hmass]=ref_hmass[iref_hmass]*lat[ib]*Zp[ib];
      
      //interpolate
      bvec temp_Mout=interpolate_multi(mass[iens],Min,mass_out);
      
      ofstream fout(combine("interpolating_ens%02d.xmg",iens).c_str());
      fout<<"@type xydy"<<endl;
      for(int iref_hmass=0;iref_hmass<nref_hmass;iref_hmass++)
	{
	  Mout[iens+iref_hmass*nens]=temp_Mout[iref_hmass];
	  cout<<iref_hmass<<" iens="<<iens<<" x="<<mass_out[iref_hmass]<<" y="<<temp_Mout[iref_hmass]<<endl;
	  fout<<ref_hmass[iref_hmass].med()<<" "<<temp_Mout[iref_hmass]<<endl;
	}
      fout<<"&"<<endl;
    }
  
  int n390=0;
  ofstream out380("380charm.xmg");
  ofstream out405("405charm.xmg");
  ofstream out420("420charm.xmg");
  out380<<"@type xydy"<<endl;
  out405<<"@type xydy"<<endl;
  out420<<"@type xydy"<<endl;
  for(int iens=0;iens<nens;iens++)
    {
      if(ibeta[iens]==0) out380<<sea_mass[iens]<<" "<<aM[iens][icombo(1,1,nmass[iens])]<<endl;
      if(ibeta[iens]==1) n390++;
      if(ibeta[iens]==2) out405<<sea_mass[iens]<<" "<<aM[iens][icombo(1,1,nmass[iens])]<<endl;
      if(ibeta[iens]==3) out420<<sea_mass[iens]<<" "<<aM[iens][icombo(1,1,nmass[iens])]<<endl;
    }
  
  int i390=0;
  double m390[n390];
  bvec M390(n390,nboot,njack);
  for(int iens=0;iens<nens;iens++)
    if(ibeta[iens]==1)
      {
	m390[i390]=sea_mass[iens];
	M390[i390]=aM[iens][icombo(1,1,nmass[iens])];
	i390++;
      }
  ofstream out390("390charm.xmg");
  bvec parfit=poly_fit(m390,M390,1,0,0.0120);
  out390<<"@    title \"Chiral extrapolation at fixed lattice spacing\""<<endl;
  out390<<"@    yaxis  label \"aM\""<<endl;
  out390<<"@    xaxis  label \"am\\ssea\""<<endl;
  out390<<"@    yaxis  label char size 1.390000"<<endl<<
    "@    xaxis  label char size 1.390000"<<endl;
    
  out390<<"@ s0 line color 9"<<endl;
  out390<<"@ s0 fill color 9"<<endl;
  out390<<"@ s0 fill type 1"<<endl;

  out390<<"@ s1 line color 4"<<endl;
  out390<<"@ s1 line linewidth 2"<<endl;

  out390<<"@ s2 line type 0"<<endl;
  out390<<"@ s2 symbol 2"<<endl;
  out390<<"@ s2 errorbar color 4"<<endl;
  out390<<"@ s2 symbol color 4"<<endl;
  out390<<"@ s2 symbol linewidth 2"<<endl;
  out390<<"@ s2 errorbar linewidth 2"<<endl;
  out390<<"@ s2 errorbar riser linewidth 2"<<endl;
  out390<<write_poly_with_error(parfit,0,0.0120);
  out390<<"@type xydy"<<endl;
  for(int iens=0;iens<nens;iens++) if(ibeta[iens]==1) out390<<sea_mass[iens]<<" "<<aM[iens][icombo(1,1,nmass[iens])]<<endl;


  Mout.write_to_binfile("intM");
  
  return 0;
}
