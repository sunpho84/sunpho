#include "include.h"
#include <TMatrixD.h>
#include <TVectorD.h>

const int deb=0;

const int njacks=16;
int T,L,Tsep;
int nlevls;
int nthetaS0;
int nthetaS1;

jack relativistic_energy(jack M,double th0)
{
  double qi=M_PI/L*th0;
  //return sqrt(M*M+3*qi*qi);
  return 2*asinh(sqrt(3*sqr(sin(qi/2))+sqr(sinh(M/2))));
}

jack compute_Q2(jack E1,double th1,jack E2,double th2)
{
  cout<<"E1="<<E1<<", E2="<<E2<<endl; 
  cout<<"Th1="<<th1<<", th2="<<th2<<endl; 
  double qi1=M_PI/L*th1;
  double qi2=M_PI/L*th2;
  
  jack Q0=E1-E2;
  double Qi=qi1-qi2;
  
  return Q0*Q0-3*Qi*Qi;
}

template<class to> to two(to Z1,to Z2,to M,int t)
{return Z1*Z2*exp(-M*L)*cosh(M*(L-t))/M;}

//load 2pts
jvec load_two_points(const char *path,int ri,int ithetaS0,int ext_ism_so_lv,int ext_ism_si_lv)
{
  int ism_si_lv=min(ext_ism_si_lv,ext_ism_so_lv);
  int ism_so_lv=max(ext_ism_si_lv,ext_ism_so_lv);
  
  jvec a(T,njacks);
  jvec b(T,njacks);
  
  a.load(path,ri+2*(0+4*(ithetaS0+nthetaS0*(ism_si_lv+nlevls*ism_so_lv))));
  b.load(path,ri+2*(3+4*(ithetaS0+nthetaS0*(ism_si_lv+nlevls*ism_so_lv))));
  
  return (a+b)/2;
}

//load the passed combination of 2pts
jvec load_two_points(const char *path,int ri,int ithetaS0,double *poly)
{
  jvec out(T,njacks);
  out=0;

  for(int ism_so=0;ism_so<nlevls;ism_so++)
    {
      jvec temp(T,njacks);
      temp=0;
      for(int ism_si=0;ism_si<nlevls;ism_si++)
	temp+=poly[ism_si]*load_two_points(path,ri,ithetaS0,ism_so,ism_si);
      out+=poly[ism_so]*temp;
    }
  
  return out;
}

//load 3pts
jvec load_three_points(const char *path,int ri,int ism_so_lv,int ism_si_lv,int ithetaS0,int ithetaS1)
{
  jvec out(T,njacks);
  
  out.load(path,ri+2*(ithetaS0+nthetaS0*(ithetaS1+nthetaS1*(ism_si_lv+nlevls*ism_so_lv))));
  
  return out;
}

//load the passed combination of 3pts
jvec load_three_points(const char *path,int ri,int ithetaS0,int ithetaS1,double *poly_source,double *poly_seq)
{
  jvec out(T,njacks);
  out=0;

  for(int ism_so=0;ism_so<nlevls;ism_so++)
    {
      jvec temp(T,njacks);
      temp=0;
      for(int ism_se=0;ism_se<nlevls;ism_se++)
	temp+=poly_seq[ism_se]*load_three_points(path,ri,ism_so,ism_se,ithetaS0,ithetaS1);
      out+=poly_source[ism_so]*temp;
    }
  
  return out;
}

int main()
{
  FILE *fin=open_file("input","r");
  
  //read T,tsep
  read_formatted_from_file_expecting((char*)&T,fin,"%d","T");
  L=T/2;
  read_formatted_from_file_expecting((char*)&Tsep,fin,"%d","Tsep");
  //read nlevls
  read_formatted_from_file_expecting((char*)&nlevls,fin,"%d","nlevls");
  
  //read poly
  double poly_P5P5[nlevls];
  double poly_VKVK[nlevls];
  read_formatted_from_file_expecting((char*)&poly_P5P5[0],fin,"%lg","poly_P5P5");
  for(int ilev=1;ilev<nlevls;ilev++) read_formatted_from_file((char*)&poly_P5P5[ilev],fin,"%lg","poly_P5P5");
  read_formatted_from_file_expecting((char*)&poly_VKVK[0],fin,"%lg","poly_VKVK");
  for(int ilev=1;ilev<nlevls;ilev++) read_formatted_from_file((char*)&poly_VKVK[ilev],fin,"%lg","poly_VKVK");
  
  //read theta
  read_formatted_from_file_expecting((char*)&nthetaS0,fin,"%d","nthetaS0");
  double thetaS0[nthetaS0];
  for(int itheta=0;itheta<nthetaS0;itheta++) read_formatted_from_file((char*)&thetaS0[itheta],fin,"%lg","thetaS0");
  read_formatted_from_file_expecting((char*)&nthetaS1,fin,"%d","nthetaS1");
  double thetaS1[nthetaS1];
  for(int itheta=0;itheta<nthetaS1;itheta++) read_formatted_from_file((char*)&thetaS1[itheta],fin,"%lg","thetaS1");
  
  //read the combo of theta to be used
  int ithetaS0,ithetaS1;
  read_formatted_from_file_expecting((char*)&ithetaS0,fin,"%d","ithetaS0");
  read_formatted_from_file_expecting((char*)&ithetaS1,fin,"%d","ithetaS1");
  
  //read tmin and tmax
  int tmin_P5P5,tmax_P5P5;
  read_formatted_from_file_expecting((char*)&tmin_P5P5,fin,"%d","tfit_P5P5");
  read_formatted_from_file((char*)&tmax_P5P5,fin,"%d","tfit_P5P5");
  int tmin_VKVK,tmax_VKVK;
  read_formatted_from_file_expecting((char*)&tmin_VKVK,fin,"%d","tfit_VKVK");
  read_formatted_from_file((char*)&tmax_VKVK,fin,"%d","tfit_VKVK");
  int tmin_R,tmax_R;
  read_formatted_from_file_expecting((char*)&tmin_R,fin,"%d","tfit_R");
  read_formatted_from_file((char*)&tmax_R,fin,"%d","tfit_R");
  
  fclose(fin);
  
  ////////////////////////////// read standing two points ///////////////////////////
  
  //load basic 2 pts
  const int RE=0,STAND=0;
  jvec corr_P5P5=load_two_points("2pts_P5P5",RE,STAND,poly_P5P5).simmetrized(1);
  jvec corr_VKVK=load_two_points("2pts_VKVK",RE,STAND,poly_VKVK).simmetrized(1);
  
  //plot
  {
    ofstream out("effective_masses.xmg");
    out<<"@type xydy"<<endl;
    out<<effective_mass(corr_P5P5)<<"&"<<endl;
    out<<effective_mass(corr_VKVK)<<"&"<<endl;
  }
  
  //fit masses
  jack M_P5,M_VK,Z_P5,Z_VK,Z2_P5,Z2_VK;
  two_pts_fit(M_P5,Z2_P5,corr_P5P5,tmin_P5P5,tmax_P5P5,"effective_mass_P5P5.xmg");
  two_pts_fit(M_VK,Z2_VK,corr_VKVK,tmin_VKVK,tmax_VKVK,"effective_mass_VKVK.xmg");
  Z_P5=sqrt(Z2_P5);
  Z_VK=sqrt(Z2_VK);
  
  //output
  cout<<"M_P5: "<<M_P5<<endl;
  cout<<"M_VK: "<<M_VK<<endl;
  
  //compute Energies according to the relativistic lattice formula
  jack Eth_P5=relativistic_energy(M_P5,thetaS1[ithetaS1]);
  jack Eth_VK=relativistic_energy(M_VK,thetaS0[ithetaS0]);
  cout<<"E_P5: "<<Eth_P5<<endl;
  cout<<"E_VK: "<<Eth_VK<<endl;
  
  //Build the time dependance
  jvec Eta_Psip_td(Tsep+1,njacks);
  for(int t=0;t<=Tsep;t++)
    Eta_Psip_td[t]=(Z_P5*Z_VK)*exp(-Eth_VK*t-Eth_P5*(Tsep-t))/(2*Eth_P5*2*Eth_VK);
  
  {
    ofstream out("Eta_Psip_td.xmg");
    out<<"@type xydy"<<endl;
    out<<Eta_Psip_td<<"&"<<endl;
  }
  
  //compute Q2
  jack Q2=compute_Q2(Eth_P5,thetaS1[ithetaS1],Eth_VK,thetaS0[ithetaS0]);
  cout<<"Q2: "<<Q2<<endl;
  
  //////////////////////////////////// read 3 points ////////////////////////////////
  
  //load 3 pts
  const int IM=1;
  jvec corr_VJVK_pt1=load_three_points("3pts_sp0_VJVK_pt1",IM,ithetaS0,ithetaS1,poly_VKVK,poly_P5P5).simmetrized(-1);
  jvec corr_VJVK_pt2=-load_three_points("3pts_sp0_VJVK_pt2",IM,ithetaS0,ithetaS1,poly_VKVK,poly_P5P5).simmetrized(-1);
  jvec corr_VJVK=(corr_VJVK_pt1+corr_VJVK_pt2)/2;

  {
    ofstream out("VJVKP5_pts.xmg");
    out<<"@type xydy"<<endl;
    out<<corr_VJVK_pt1<<"&"<<endl;
    out<<corr_VJVK_pt2<<"&"<<endl;
    out<<corr_VJVK<<"&"<<endl;
  }

  //divide by the time dependance and fit
  jvec Eta_V_Psip=corr_VJVK/Eta_Psip_td;
  jack R=constant_fit(Eta_V_Psip,tmin_R,tmax_R,"R.xmg");
  
  return 0;
}
