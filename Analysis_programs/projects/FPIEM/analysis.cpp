#include "../../src/include.h"

int T,L;
double th;
int njacks;
double Zv;

void read_input()
{
  FILE *fin=open_file("input","r");
  
  read_formatted_from_file_expecting((char*)&T,fin,"%d","T");
  L=T/2;
  read_formatted_from_file_expecting((char*)&njacks,fin,"%d","NJacks");
  read_formatted_from_file_expecting((char*)&th,fin,"%lg","Theta");
  read_formatted_from_file_expecting((char*)&Zv,fin,"%lg","Zv");
  
  fclose(fin);
}

int main()
{
  read_input();
  
  debug_fit=0;
  
  jvec plus=jvec_load(combine("mes_contr_P5P5_Pi_PLUS_%lg",th),T,njacks,0).simmetrized(1);
  jack E(njacks),Z2(njacks);
  two_pts_fit(E,Z2,plus,7,40,"plus.xmg","/dev/null");
  cout<<"Energy: "<<smart_print(E)<<endl;
  double pi=M_PI*th/L;
  
  jvec rest=jvec_load("mes_contr_P5P5_Pi_REST",T,njacks,0).simmetrized(1);
  jack M(njacks);
  two_pts_fit(M,Z2,rest,7,L,"rest.xmg","/dev/null");
  cout<<"Mass: "<<smart_print(M)<<endl;
  E=latt_e(M,pi);
  cout<<"Energy reco: "<<smart_print(E)<<endl;
  
  jvec mel=-jvec_load(combine("mes_contr_V0P5_Pi_MEL_Pi_%lg",th),T,njacks,0).simmetrized(-1);
  jvec m=(mel+mel.inverted())/2;
  mel.print_to_file("mel.xmg");
  m.print_to_file("m.xmg");
  
  TH_two_pts_fit=L;
  jvec rat_an=Zv*m/(fun_two_pts_migrad_fit(Z2,E,L)/2);
  jvec rat_nu=Zv*m/(plus[L]/2);
  rat_an.print_to_file("rat_an.xmg");
  rat_nu.print_to_file("rat_nu.xmg");
  if(Zv==1)
    {
      cout<<"Zv: "<<smart_print(1/rat_an[T/4])<<endl;
      cout<<"Zv: "<<smart_print(1/rat_nu[T/4])<<endl;
    }
  else
    {
      cout<<"f: "<<smart_print(rat_an[T/4])<<endl;
      cout<<"f: "<<smart_print(rat_nu[T/4])<<endl;
    }
  
  return 0;
}
