#include "include.h"
#include "corr.cpp"
//used to find pi and K masses

int njacks;
double aml,ams;
int tmin;

int T,TH,L;
int ncombos=6;

int main()
{
  debug_load=debug_fit=false;

  FILE *input_file=open_file("input","r");
  read_formatted_from_file_expecting((char*)(&T),input_file,"%d","T");
  read_formatted_from_file_expecting((char*)(&L),input_file,"%d","L");
  read_formatted_from_file_expecting((char*)(&njacks),input_file,"%d","njacks");
  read_formatted_from_file_expecting((char*)(&aml),input_file,"%lg","aml");
  read_formatted_from_file_expecting((char*)(&ams),input_file,"%lg","ams");
  read_formatted_from_file_expecting((char*)(&tmin),input_file,"%d","tmin");
  TH=T/2;
  
  cout<<"Njacks: "<<njacks<<endl;
  
  jvec pion=jvec_load("data",T,njacks,0).simmetrized(1);
  jvec kaon=jvec_load("data",T,njacks,3).simmetrized(1);

  jack aM_Pi(njacks),Z2_Pi(njacks);
  two_pts_fit(aM_Pi,Z2_Pi,pion,tmin,TH-1,"pion.xmg");
  jack aM_K(njacks),Z2_K(njacks);
  two_pts_fit(aM_K,Z2_K,kaon,tmin,TH-1,"kaon.xmg");
      
  jack af_Pi=sqrt(Z2_Pi)*2*aml/(2*(aM_Pi)*aM_Pi);
  jack af_K=sqrt(Z2_K)*(aml+ams)/(2*(aM_K)*aM_K);

  cout<<"Kaon/Pion mass: "<<smart_print(aM_K/aM_Pi)<<" [phys: 3.662]"<<endl;

  jack a=aM_Pi/0.135;
  cout<<"Pion mass in lattice spacing units: "<<smart_print(aM_Pi)<<endl;
  cout<<"Pion decc in lattice spacing units: "<<smart_print(af_Pi)<<endl;
  cout<<"Lattice spacing=(aM_Pi)/M_Pi: "<<smart_print(a)<<" GeV^-1, "
      <<smart_print(a*0.197)<<" fm"<<endl;
  cout<<"M_Pi*L: "<<smart_print(aM_Pi*L)<<endl;
  
  a/=1000;

  jack M_Pi=aM_Pi/a;
  jack M_K=aM_K/a;
  jack f_Pi=af_Pi/a;
  jack f_K=af_K/a;

  cout<<"M_Pi: "<<smart_print(M_Pi)<<" MeV [taken as input to fix latt.]"<<endl;
  cout<<"M_K: "<<smart_print(M_K)<<" MeV [phys: 495 MeV]"<<endl;
  
  cout<<"F_Pi: "<<smart_print(f_Pi)<<" MeV [phys: 130.2(1.4)]"<<endl;
  cout<<"F_K: "<<smart_print(f_K)<<" MeV [phys: 155.5]"<<endl;

  cerr<<smart_print(a*1000*0.197)<<" "<<smart_print(M_K/M_Pi)<<" "<<smart_print(f_K/M_Pi)<<endl;

  //compute correction
  double vm,vf;
  corr_vf(vm,vf,aM_Pi.med(),af_Pi.med(),L);
  cout<<"mass corr: "<<(1-vm)*100<<"%, f corr: "<<(1-vf)*100<<"%"<<endl;
  
  cout<<"F_Pi/M_Pi: "<<smart_print(f_Pi/M_Pi)<<", corr: "<<smart_print(f_Pi*vm/vf/M_Pi)<<", phys: "<<0.965<<endl;
  
  return 0;
}
