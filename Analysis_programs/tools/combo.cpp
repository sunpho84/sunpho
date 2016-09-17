#include <iostream>
#include <stdlib.h>

using namespace std;

int main(int narg,char **arg)
{
  if(narg<4 ||(narg-4)%2)
    {
      cerr<<"Use "<<arg[0]<<" out file1 coef1 [file2 coef2 ...]"<<endl;
      exit(1);
    }
  
  FILE *fout=fopen(arg[1],"w");
  if(fout==NULL)
    {
      cerr<<"Error opening "<<arg[1]<<endl;
      exit(1);
    }
  
  int nconts=(narg-2)/2;
  FILE *fin[nconts];
  double coef[nconts];
  
  for(int icont=0;icont<nconts;icont++)
    {
      fin[icont]=fopen(arg[2+icont*2],"r");
      if(fin[icont]==NULL)
	{
	  cerr<<"Error opening "<<arg[2+icont*2]<<endl;
	  exit(1);
	}
      
      if(sscanf(arg[3+icont*2],"%lg",&(coef[icont]))!=1)
	{
	  cerr<<"Error scanning "<<arg[3+icont*2]<<endl;
          exit(1);
	}
      cout<<"Coef"<<icont<<" = "<<coef[icont]<<endl;
    }
  
  int nr;
  do
    {
      double to=0;
      
      nr=0;
      for(int icont=0;icont<nconts;icont++)
	{
	  double t;
	  nr+=fread((void*)&t,sizeof(double),1,fin[icont]);
	  to+=t*coef[icont];
	}
      cout<<to<<endl;
      if(nr!=0 && nr!=nconts)
	{
	  cerr<<"Error reading, found "<<nr<<endl;
	  exit(1);
	}
      
      if(nr==nconts)
	{
	  if(fwrite((void*)&to,sizeof(double),1,fout)!=1)
	    {
	      cerr<<"Error writing"<<endl;
	      exit(1);
	    }
	}
      else
	cout<<"Reached end of file"<<endl;
    }
  while(nr!=0);
  
  for(int icont=0;icont<nconts;icont++)
    fclose(fin[icont]);
  
  fclose(fout);
  
  return 0;
}
