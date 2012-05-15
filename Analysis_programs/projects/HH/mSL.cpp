#include <include.h>
#include <iostream>
#include <sstream>

using namespace std;

int nmass;
int njack;
double *mass;
int T,TH;
int sour_nsmear;
char base_path[1024];
char corr_name[1024],out_file[1024];
int tmin[2],tmax[2];

int icombo(int im1,int im2,int r1,int r2,int ri)
{
  int ic=ri+2*(r1+2*(im1+nmass*(r2+2*im2)));

  return ic;
}

jvec load_corr(char *base_path,int sm_lev,int im1,int im2,int r1,int r2,const char *obs_name)
{
  char path[1024];
  int tag[2]={0,sour_nsmear};
  sprintf(path,"%s/2pts_%s_%02d_%02d",base_path,obs_name,sour_nsmear,tag[sm_lev]);

  jvec l=(jvec_load(path,T,njack,icombo(im1,im2,!r1,!r2,0))+
	  jvec_load(path,T,njack,icombo(im1,im2, r1, r2,0)))/2;
  
  return l;
}

void read_ensemble_pars(const char *data_list_file)
{
  FILE *input=open_file(data_list_file,"r");
  
  read_formatted_from_file_expecting(base_path,input,"%s","base_path");
  read_formatted_from_file_expecting((char*)&T,input,"%d","T");
  read_formatted_from_file_expecting((char*)&njack,input,"%d","njack");
  int ibeta;
  read_formatted_from_file_expecting((char*)&ibeta,input,"%d","Beta");
  read_formatted_from_file_expecting((char*)&sour_nsmear,input,"%d","sour_nsmear");
  
  read_formatted_from_file_expecting((char*)&nmass,input,"%d","nmass");
  expect_string_from_file(input,"mass_list");
  mass=(double*)malloc(sizeof(double)*nmass);
  for(int imass=0;imass<nmass;imass++) read_formatted_from_file((char*)&(mass[imass]),input,"%lg","mass");
  
  fclose(input);
}

void read_pars(const char *input)
{
  FILE *fin=open_file(input,"r");
  char data_list_file[256];
  read_formatted_from_file_expecting(data_list_file,fin,"%s","data_list_file");

  read_formatted_from_file_expecting(corr_name,fin,"%s","corr_name");
  
  read_formatted_from_file_expecting((char*)(&tmin[0]),fin,"%d","tint_low");
  read_formatted_from_file((char*)(&tmax[0]),fin,"%d","tint_low");
  
  read_formatted_from_file_expecting((char*)(&tmin[1]),fin,"%d","tint_high");
  read_formatted_from_file((char*)(&tmax[1]),fin,"%d","tint_high");
  
  read_formatted_from_file_expecting(out_file,fin,"%s","out_file");
  
  fclose(fin);
  
  read_ensemble_pars(data_list_file);
  
  TH=T/2;
}

int main(int narg,char **arg)
{
  read_pars("input");
  jvec c;
  //loop over heavier mass
  int ncombo=nmass*nmass;
  jvec aM(ncombo,njack);
  for(int im1=0;im1<nmass;im1++)
    for(int im2=0;im2<nmass;im2++)
      {
	int ic=im1*nmass+im2;
	
	//load the corrs
	c=load_corr(base_path,0,im1,im2,0,0,corr_name); //TM
	
	double f=sqrt((im1*im1+im2*im2)/2)/(nmass-1);
	int tmi=int((tmin[1]-tmin[0])*f+tmin[0]+0.5);
	int tma=int((tmax[1]-tmax[0])*f+tmax[0]+0.5);
	cout<<f<<endl;
	aM.data[ic]=constant_fit(effective_mass(c.simmetrized(1)),tmi,tma,combine("eff_plot_%02d_%02d.xmg",im1,im2).c_str());
	
	ofstream out(combine("plot_%02d_%02d.xmg",im1,im2).c_str());
	out<<c<<endl;
	cout<<"im1:"<<im1<<" im2:"<<im2<<" ic:"<<ic<<"/"<<ncombo<<" aM: "<<aM[ic]<<
	  " c[0]/c["<<TH<<"]="<<c[0]/c[TH]<<
	  endl;
      }
  
  aM.write_to_binfile("M");

  return 0;
}
