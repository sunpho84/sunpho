#include "include.h"
#include <cmath>

void scan(std::vector<double> &data,int narg,char **arg)
{
  FILE *fin;
  if(narg<2||(strcasecmp(arg[1],"-")==0)) fin=stdin;
  else fin=open_file(arg[1],"r");
  
  //////////////////////////////////////////
  
  double temp;
  while(fscanf(fin,"%lg",&temp)==1) data.push_back(temp);
      
  //////////////////////////////////////////
  
  if(narg<2) fclose(fin);
}

void err_block(std::vector<double> &data,int block_size,double &err,double &err_err)
{
  double s2=0,s=0;
  int n=0;
  int data_size=data.size();
  for(int i=0;i<=data_size-block_size;i+=block_size)
    {
      double p=0;
      for(int j=i;j<i+block_size;j++) p+=data[j];
      p/=block_size;
      s2+=p*p;
      s+=p;
      n++;
    }
  
  s/=n;
  s2/=n;
  
  s2-=s*s;
    
  err=sqrt(s2/(n-1));
  err_err=err*sqrt(1-2*exp(2*(lgamma(n/2.0)-lgamma((n-1.0)/2)))/(n-1));
}

int main(int narg,char **arg)
{
  std::vector<double> data;
  scan(data,narg,arg);
  int data_size=data.size();
  cerr<<"data size: "<<data_size<<endl;
  
  int ib=0;
  for(int block_size=1;block_size<data_size/4;block_size*=2)
    {
      double err,err_err;
      err_block(data,block_size,err,err_err);
      cout<<ib++<<" "<<err<<" "<<err_err<<endl;
    }
  
  return 0;
}
