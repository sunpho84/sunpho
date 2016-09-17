#include "include.h"

int njacks;
int T,term_conf;
char data_path[200];
int conf_min,conf_max;
jvec *data;
int tmin_fit[3],tmax_fit[3];

void read_input(const char *input_path)
{
  FILE *fin=open_file(input_path,"r");
  
  read_formatted_from_file_expecting((char*)&njacks,fin,"%d","njacks");
  read_formatted_from_file_expecting(data_path,fin,"%s","data_path");
  read_formatted_from_file_expecting((char*)&term_conf,fin,"%d","term_conf");
  read_formatted_from_file_expecting((char*)&T,fin,"%d","T");
  for(int icombo=0;icombo<3;icombo++)
    {
      read_formatted_from_file_expecting((char*)&(tmin_fit[icombo]),fin,"%d",combine("tfit_int%d",icombo).c_str());
      read_formatted_from_file((char*)&(tmax_fit[icombo]),fin,"%d","tfit_int");
    }
  
  fclose(fin);
}

void read_data()
{
  FILE *fin=open_file(data_path,"r");
      
  int rc;
  int ncorr=0;
  do
    {      
      int iconf,read_t;
      double m1,m2,fuf;
      rc=fscanf(fin," # iconf %d , m1 = %lg , m2 = %lg\n",&iconf,&m1,&m2);
      if(rc>0)
	{
	  if(ncorr==0||iconf<conf_min) conf_min=iconf;
          if(ncorr==0||iconf>conf_max) conf_max=iconf;
	  
	  for(int t=0;t<T;t++)
	    if(fscanf(fin,"%d %lg",&read_t,&fuf)!=2||read_t!=t) crash("incomplete correlation");
	  ncorr++;
	}
    }
  while(rc>0);

  int ncombo=ncorr/(conf_max-conf_min);
  
  printf("ncorr: %d, conf_min: %d, conf_max: %d, ncombo: %d\n",ncorr,conf_min,conf_max,ncombo);
  
  //get back
  fseek(fin,0,SEEK_SET);
  
  conf_min=max(term_conf,conf_min);
  conf_max=max(term_conf,conf_max);
  int nconf_ava=conf_max-conf_min+1;
  int clust_size=nconf_ava/njacks;
  int nconf=clust_size*njacks;
  conf_max=conf_min+nconf;
  
  printf("conf used: %d %d\n",conf_min,conf_max);
  
  //allocate data
  data=(jvec*)malloc(ncombo*sizeof(jvec));
  for(int icombo=0;icombo<ncombo;icombo++)
    {
      data[icombo]=jvec(T,njacks);
      data[icombo]*=0;
    }
  
  //clusterize
  ncorr=0;
  do
    {      
      int iconf,read_t;
      double m1,m2,fuf;
      rc=fscanf(fin," # iconf %d , m1 = %lg , m2 = %lg\n",&iconf,&m1,&m2);
      if(rc>0)
	{
	  int icombo=ncorr%3;
	  for(int t=0;t<T;t++)
	    {
	      if(fscanf(fin,"%d %lg",&read_t,&fuf)!=2||read_t!=t) crash("incomplete correlation");
	      if(iconf>=conf_min && iconf<conf_max)
		{
		  int iclust=(iconf-conf_min)/clust_size;
		  data[icombo].data[t].data[iclust]+=fuf;
		}
	    }
	  ncorr++;
	}
    }
  while(rc>0);
  
  fclose(fin);
  
  //jack
  for(int icombo=0;icombo<ncombo;icombo++)
    for(int t=0;t<T;t++)
      {
	for(int iclust=0;iclust<njacks;iclust++)
          data[icombo].data[t].data[njacks]+=data[icombo].data[t].data[iclust];
        
        for(int iclust=0;iclust<njacks;iclust++)
          data[icombo].data[t].data[iclust]=(data[icombo].data[t].data[njacks]-data[icombo].data[t].data[iclust])/(nconf-clust_size);
        
        data[icombo].data[t].data[njacks]/=nconf;
      }
}


int main(int narg,char **arg)
{
  if(narg<2) crash("use %s input",arg[0]);
  read_input(arg[1]);
  
  read_data();
  
  ofstream out("pion.xmg");
  jvec M(3,njacks);
  for(int i=0;i<3;i++)
    {
      jvec cm=effective_mass(data[i].simmetrized(1));
      M[i]=constant_fit(cm,tmin_fit[i],tmax_fit[i]);
      
      //write the plot
      out<<"@type xydy"<<endl;
      out<<cm;
      out<<"&"<<endl;
      out<<write_constant_with_error(M[i],tmin_fit[i],tmax_fit[i]);
      out<<"&"<<endl;
    }
  
  cout<<smart_print(M[1]/M[0])<<endl;
  cout<<494.0/135<<endl;
  
  return 0;
}
