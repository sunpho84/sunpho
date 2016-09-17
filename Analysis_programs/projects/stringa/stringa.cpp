#include "include.h"

char path[100];
int DMin,DMax,TMin,TMax,nconfs,nsme,njacks;
int dD,dT;

void read_input(char *path_input)
{
  FILE *fin=open_file(path_input,"r");
  
  read_formatted_from_file_expecting(path,fin,"%s","measures_path");
  read_formatted_from_file_expecting((char*)&DMin,fin,"%d","DMin");
  read_formatted_from_file_expecting((char*)&DMax,fin,"%d","DMax");
  read_formatted_from_file_expecting((char*)&TMin,fin,"%d","TMin");
  read_formatted_from_file_expecting((char*)&TMax,fin,"%d","TMax");
  read_formatted_from_file_expecting((char*)&nconfs,fin,"%d","NConfs");
  read_formatted_from_file_expecting((char*)&nsme,fin,"%d","NSme");
  
  dD=DMax-DMin+1;
  dT=TMax-TMin+1;
  njacks=nconfs;

  fclose(fin);
}

int main(int narg,char **arg)
{
  if(narg<2) crash("Use %s file",arg[0]);
  
  read_input(arg[1]);
  
  ifstream fin(path);
  double *data_raw=new double[nconfs*nsme*dT*dD*3];
  
  int nr=0;
  for(int iconf=0;iconf<nconfs;iconf++)
    for(int isme=0;isme<nsme;isme++)
      for(int t=0;t<dT;t++)
	for(int d=0;d<dD;d++)
	  {
	    int jconf,jsme,u,e;
	    if(!(fin>>jconf>>jsme>>u>>e)) crash("reading nr=%d iconf %d isme %d t %d d %d",nr,iconf,isme,t,d);
	    
	    //if(jconf!=iconf) crash("nr=%d jconf=%d!=iconf=%d",nr,jconf,iconf);
	    if(jsme!=isme) crash("nr=%d jsme=%d!=isme=%d",nr,jsme,isme);
	    if(u!=t+TMin) crash("nr=%d u=%d!=t=%d",nr,u,t+TMin);
	    if(e!=d+DMin) crash("nr=%d e=%d!=d=%d",nr,e,d+DMin);
	    
	    for(int i=0;i<3;i++) //x,y,z
	      if(!(fin>>data_raw[i+3*(d+dD*(t+dT*(isme+nsme*iconf)))])) crash("reading xyz=%d",i);
	    
	    nr++;	      
	  }
  
  double dum;
  if(fin>>dum) crash("not reached EOF");
  
  //for(int isme=0;isme<nsme;isme++)
  int isme=4;
    {
      jvec eff[3]={jvec(dT-1,njacks),jvec(dT-1,njacks),jvec(dT-1,njacks)};
	  
      for(int d=0;d<dD;d++)
	{
	  for(int i=0;i<3;i++)
	    {
	      //create all clusters
	      jvec corr(dT,njacks);
	      for(int iconf=0;iconf<nconfs;iconf++)
		for(int t=0;t<dT;t++)
		  corr[t].data[iconf]=data_raw[i+3*(d+dD*(t+dT*(isme+nsme*iconf)))];
	      
	      //create jacknives
	      for(int t=0;t<dT;t++)
		{
		  corr[t][nconfs]=0;
		  for(int iconf=0;iconf<nconfs;iconf++) corr[t][nconfs]+=corr[t][iconf];
		  for(int iconf=0;iconf<nconfs;iconf++) corr[t][iconf]=(corr[t][nconfs]-corr[t][iconf])/(njacks-1);
		  corr[t][nconfs]/=nconfs;
		}
	      
	      eff[i]=aperiodic_effective_mass(corr);

	      ofstream fout(combine("dir%d_d%02d",i,d+DMin).c_str());
	      fout<<"@type xydy"<<endl;
	      for(int t=0;t<dT-1;t++)
		fout<<t+TMin<<" "<<eff[i][t]<<endl;
	      fout.close();
	    }
	  
	  //plot diff
	  ofstream fout(combine("d%02d",d).c_str());
	  fout<<"@type xydy"<<endl;
	  for(int t=0;t<dT-1;t++)
	    fout<<t+TMin<<" "<<(eff[0][t]+eff[1][t])/2-eff[2][t]<<endl;
	  
	  fout.close();
	}
    }
  
  delete[] data_raw;

  return 0;
}
