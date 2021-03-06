#include "include.h"

int nens;
double *aml;
int *T;
int njack=15;

int main()
{
  FILE *fin=open_file("input","r");
  
  read_formatted_from_file_expecting((char*)&nens,fin,"%d","nensemble");
  aml=new double[nens];
  T=new int[nens];
  jvec aM(nens,njack);
  jvec Z2(nens,njack);
  jvec af(nens,njack);
  jvec mist(nens,njack);
  
  for(int iens=0;iens<nens;iens++)
    {
      read_formatted_from_file((char*)&(aml[iens]),fin,"%lg","aml");
      read_formatted_from_file((char*)&(T[iens]),fin,"%d","T");
      
      char path[1024];
      read_formatted_from_file(path,fin,"%s","path");
      
      aM.data[iens].load(path,0);
      Z2.data[iens].load(path,1);
      
      af[iens]=2*aml[iens]*sqrt(Z2[iens])/sqr(aM[iens]);
      mist[iens]=sqr(aM[iens])*pow(af[iens],1.5)/aml[iens];
    }
  
  jack m,q;
  linear_fit(m,q,aml,mist,0,0.01,"mist.xmg");
  cout<<"chir: "<<q<<endl;
  
  for(int iens=0;iens<nens;iens++) cout<<aml[iens]<<" "<<mist[iens]<<endl;
  
  {
    ofstream fout("M2pi_fr_aml.xmg");
    fout<<"@type xydy"<<endl;
    for(int iens=0;iens<nens;iens++) fout<<aml[iens]<<" "<<aM[iens]*aM[iens]/aml[iens]<<endl;
  }
  
  fclose(fin);
  
  return 0;
}
