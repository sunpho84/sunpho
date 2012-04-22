#include "include.h"

int njack=16;

int main()
{
  int nsea=6;
  int L[6] ={32,32,24,24,24,24};
  int am_sea[6]={30,40,40,64,85,100};
  double am[15]={0,0.0159,0.0177,0.0195,0.1828,0.2150,0.2529,0.2974,0.3498,0.4114,0.4839,0.5691,0.6694,0.7873,0.9260};
  int nl=4;
  int nm=15;
  
  jvec tempM[6];
  jvec tempZ[6];
  for(int il=0;il<6;il++)
    {
      string a=combine("/home/francesco/QCD/LAVORI/SMEARING_BK_RUNS/FITTED_MASS_Z/P5P5/3.90/%2d/0.0%03d/results",L[il],am_sea[il]);

      tempM[il].create(nl*nm,njack);
      tempZ[il].create(nl*nm,njack);
      
      tempM[il].load(a.c_str(),0);
      tempZ[il].load(a.c_str(),1);
    }
  
  for(int ic=nl;ic<nm;ic++)
    {
      ofstream out(combine("Phi%02d.xmg",ic).c_str());
      out<<"@type xydy"<<endl;
      
      double a=0.086/0.197;
      
      for(int isea=0;isea<nsea;isea++)
	{
	  jack aM=tempM[isea][ic];
	  jack Z=tempZ[isea][ic];
	  
	  am[0]=am_sea[isea]/10000.0;
	  
	  jack aMs=tempM[isea][ic+nm*3];
	  jack Zs=tempZ[isea][ic+nm*3];
	  
	  jack af=(am[0]+am[ic])*sqrt(Z)/(sinh(aM)*aM);
	  jack aphi=af*sqrt(aM);
	  jack afs=(am[3]+am[ic])*sqrt(Zs)/(sinh(aMs)*aMs);
	  jack aphis=afs*sqrt(aMs);
	  cout<<aM/a<<" "<<af/a<<" "<<endl;
	  cout<<am[0]<<" "<<aphi/pow(a,1.5)<<endl;
	  out<<am[0]<<" "<<aphis/aphi<<endl;
	} 
      cout<<endl;
      out.close();
    }
  
  return 0;
}
