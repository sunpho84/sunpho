#include "include.h"

int nboot=100;
int njack=16;
int nm=11;
double mass[11]={0.984873,1.15836,1.36255,1.6023,1.88462,2.2165,2.60711,3.06615,3.60653,4.24174,4.98902};
double inv_mass[11];
bvec phih(nm,nboot,njack);

int main()
{
  FILE *an_input_file=open_file("analysis_pars","r");
  char chiral_data[1024];
  read_formatted_from_file_expecting(chiral_data,an_input_file,"%s","chiral_data");
  fclose(an_input_file);
  
  phih.load(chiral_data,0);
  
  for(int imass=0;imass<nm;imass++) inv_mass[imass]=1/mass[imass];
  
  int last=0;
  while(!isnan(phih[last].err()) && last<nm)
    last++;
  
  bvec par(3,nboot,njack);
  for(int iboot=0;iboot<nboot+1;iboot++)
    {
      double y[3];
      double *x=inv_mass+last-3;
      y[0]=phih[last-3][iboot];
      y[1]=phih[last-2][iboot];
      y[2]=phih[last-1][iboot];
      parabolic_spline(par.data[0].data[iboot],par.data[1].data[iboot],par.data[2].data[iboot],x,y);
    }
  
  grace out("phi_vs_m.xmg");
  out.fout<<"@type xydy"<<endl;
  for(int im=0;im<nm;im++) out.fout<<1/mass[im]<<" "<<phih[im]<<endl;

  cout<<"pusc"<<endl;
  
  return 0;
}
