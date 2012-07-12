#pragma once

#include "include.h"

void read_ensemble_pars(char *base_path,int &T,int &TSEP,int &ibeta,int &nmass,double *&mass,int &iml_un,int &nlight,int &ntheta,double *&theta,int &ist_th,int &njack,const char *data_list_file)
{
  FILE *input=open_file(data_list_file,"r");
  
  read_formatted_from_file_expecting(base_path,input,"%s","base_path");
  read_formatted_from_file_expecting((char*)&T,input,"%d","T");
  read_formatted_from_file_expecting((char*)&TSEP,input,"%d","TSEP");
  read_formatted_from_file_expecting((char*)&ibeta,input,"%d","Beta");
  
  read_formatted_from_file_expecting((char*)&nmass,input,"%d","nmass");
  expect_string_from_file(input,"mass_list");
  mass=(double*)malloc(sizeof(double)*nmass);
  for(int imass=0;imass<nmass;imass++) read_formatted_from_file((char*)&(mass[imass]),input,"%lg","mass");
  read_formatted_from_file_expecting((char*)&iml_un,input,"%d","iml_un");
  read_formatted_from_file_expecting((char*)&nlight,input,"%d","nlight");
  
  read_formatted_from_file_expecting((char*)&ntheta,input,"%d","ntheta");
  expect_string_from_file(input,"theta_list");
  theta=(double*)malloc(sizeof(double)*ntheta);
  for(int itheta=0;itheta<ntheta;itheta++) read_formatted_from_file((char*)&(theta[itheta]),input,"%lg","theta");
  
  read_formatted_from_file_expecting((char*)&ist_th,input,"%d","ist_th");
  
  read_formatted_from_file_expecting((char*)&njack,input,"%d","njack");
  
  fclose(input);
}

//load all data
void load_ensembles_list(char **&base_corrs_path,char **&ens_name,int &nens,int *&T,int *&TSEP,int *&ibeta,int *&nmass,double **&mass,int *&iml_un,int *&nlight,int *&ntheta,double **&theta,int *&ist_th,int *&njack,const char *ens_list_path)
{
  FILE *ens_list_file=open_file(ens_list_path,"r");
  
  //load the number of ensembles
  read_formatted_from_file_expecting((char*)(&nens),ens_list_file,"%d","nens");
  
  //base address of all correlations
  char base_ens_path[1024];
  read_formatted_from_file_expecting((char*)(&base_ens_path),ens_list_file,"%s","base_ens_path");
  
  //ensemble parameters
  base_corrs_path=(char**)malloc(nens*sizeof(char*));
  ens_name=(char**)malloc(nens*sizeof(char*));
  T=(int*)malloc(nens*sizeof(int));
  TSEP=(int*)malloc(nens*sizeof(int));
  ibeta=(int*)malloc(nens*sizeof(int));
  nmass=(int*)malloc(nens*sizeof(int));
  mass=(double**)malloc(nens*sizeof(double*));
  iml_un=(int*)malloc(nens*sizeof(int));
  ntheta=(int*)malloc(nens*sizeof(int));
  theta=(double**)malloc(nens*sizeof(double*));
  ist_th=(int*)malloc(nens*sizeof(int));
  nlight=(int*)malloc(nens*sizeof(int));
  njack=(int*)malloc(nens*sizeof(int));
  
  //Loop over all ensemble. Ensemble specification is supposed to be stored in a file named [base_ens_path]/[ens_path]/data_list
  for(int iens=0;iens<nens;iens++)
    {
      //allocate ensemble name and correlations path
      base_corrs_path[iens]=(char*)malloc(1024);
      ens_name[iens]=(char*)malloc(1024);
      
      //load ensemble name and parameters
      read_formatted_from_file(ens_name[iens],ens_list_file,"%s","ens_name");
      read_ensemble_pars(base_corrs_path[iens],T[iens],TSEP[iens],ibeta[iens],nmass[iens],mass[iens],iml_un[iens],nlight[iens],ntheta[iens],theta[iens],ist_th[iens],njack[iens],combine("%s/%s/data_list",base_ens_path,ens_name[iens]).c_str());
    }
  
  fclose(ens_list_file);
}
