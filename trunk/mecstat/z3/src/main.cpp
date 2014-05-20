#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "action.hpp"
#include "close.hpp"
#include "data.hpp"
#include "debug.hpp"
#include "geometry.hpp"
#include "init.hpp"
#include "macros.hpp"
#include "meta.hpp"
#include "random.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

//update a particular site
void update_heat(int s,int itraj)
{
  //get current state
  int e=phi[s];
  
  //compute the case
  int ntypes[3]={0,0,0};
  for(int mu=0;mu<NDIMS;mu++)
    {
      ntypes[phi[neighup(s,mu)]]++;
      ntypes[phi[neighdw(s,mu)]]++;
    }
  int icase=ntypes[0]*(2*NDIMS+1)+ntypes[1];
  
  //compute transition probabilities
  double act[3],ave_act=0;;
  for(z3_t d=0;d<3;d++)
    {
      double base_act=act_contr_tab[icase][d];
      double meta_act=compute_meta_potential(glb_N+contr_N[d]-contr_N[e],glb_N0+contr_N0[d]-contr_N0[e],itraj);
      
      act[d]=base_act+meta_act;
      ave_act+=act[d];
    }
  ave_act/=3;
  double tr_prob[3];
  for(z3_t d=0;d<3;d++) tr_prob[d]=exp(-act[d]+ave_act);
  double tr_extr=get_unif_double(tr_prob[0]+tr_prob[1]+tr_prob[2]);
  
  //select
  int f;
  if(tr_extr<tr_prob[0]) f=0;
  else
    if(tr_extr<tr_prob[0]+tr_prob[1]) f=1;
    else f=2;
  
  //remove contribution
  glb_N-=contr_N[e];
  glb_N0-=contr_N0[e];
  glb_ntypes[e]--;
  glb_nequals-=ntypes[e];
  
  //add contribution
  phi[s]=f;
  glb_N+=contr_N[f];
  glb_N0+=contr_N0[f];
  glb_ntypes[f]++;
  glb_nequals+=ntypes[f];
}
void update_heat(int itraj)
{for(int s=0;s<V;s++) update_heat(s,itraj);}

//update a particular site
void update_metro(int s,int itraj)
{
  z3_t a=phi[s];
  double base_act_in=compute_action(phi);
  double meta_act_in=compute_meta_potential(glb_N,glb_N0,itraj);
  double act_in=base_act_in+meta_act_in;
  
  int b=phi[s]=get_random_z3();
  
  double base_act_fin=compute_action(phi);
  double meta_act_fin=compute_meta_potential(glb_N+contr_N[b]-contr_N[a],glb_N0+contr_N0[b]-contr_N0[a],itraj);
  double act_fin=base_act_fin+meta_act_fin;
  
  double delta_act=act_fin-act_in;
  double pacc=exp(-delta_act);
  double pext=get_unif_double(1);
  if(pext>=pacc) phi[s]=a;

  /*else
    {
      //remove contribution
      glb_N-=contr_N[a];
      glb_N0-=contr_N0[a];
      glb_ntypes[a]--;
      glb_nequals-=ntypes[a];
      
      //add contribution
      glb_N+=contr_N[b];
      glb_N0+=contr_N0[b];
      glb_ntypes[b]++;
      glb_nequals+=ntypes[b];
    }
  */
}
void update_metro(int itraj)
{
  for(int s=0;s<V;s++) update_metro(s,itraj);
  compute_action(phi);
  glb_N=compute_N(phi);
  glb_N0=compute_N0(phi);
}

int main(int narg,char **arg)
{
  init(HOT,100);
  
  //init occupation histograms
  vector<int> occ_N(2*V+1),occ_N0(V+1);
  for(int i=0;i<=2*V;i++) occ_N[i]=0;
  for(int i=0;i<=V;i++) occ_N0[i]=0;
  
  for(int iterm=0;iterm<nterm;iterm++) update_heat(0);
  
  //update
  ofstream data_file("/tmp/heat/data.xmg");
  ofstream energy_file("/tmp/heat/energy.xmg");
  for(int itraj=0;itraj<ntraj;itraj++)
    {
      update_heat(itraj);
      
      //mark down occupations
      data_N[itraj].N=glb_N;
      data_N[itraj].N0=glb_N0;
      
      //mark them in history
      occ_N[glb_N+V]++;
      occ_N0[glb_N0]++;
      //increase_meta_potentials(glb_N,glb_N0);
      
      //write data_N
      data_file<<data_N[itraj].N<<" "<<data_N[itraj].N0<<endl;
      energy_file<<compute_energy_internal()/V<<" "<<endl;
    }
  
  //write the histogram on N
  ofstream histo_file_N("/tmp/heat/histo_N.xmg");
  for(int i=0;i<=2*V;i++) if(occ_N[i]) histo_file_N<<i-V<<" "<<(double)occ_N[i]/ntraj<<endl;
  
  //write the histogram on N0
  ofstream histo_file_N0("/tmp/heat/histo_N0.xmg");
  for(int i=0;i<=V;i++) if(occ_N0[i]) histo_file_N0<<i<<" "<<(double)occ_N0[i]/ntraj<<endl;

  //write the meta potential on N
  //ofstream meta_file_N("/tmp/heat/meta_pot_N.xmg");
  //double off_N=*(std::max_element(meta_pot_N.begin(),meta_pot_N.end()));
  //for(int i=0;i<=2*V;i++)
  //{
  //meta_pot_N[i]-=off_N;
  //if(occ_N[i]) meta_file_N<<i-V<<" "<<meta_pot_N[i]<<endl;
  //}
  
  //write the meta potential on N0
  //ofstream meta_file_N0("/tmp/heat/meta_pot_N0.xmg");
  //double off_N0=*(std::max_element(meta_pot_N0.begin(),meta_pot_N0.end()));
  //for(int i=0;i<=V;i++)
  //{
  //meta_pot_N0[i]-=off_N0;
  //if(occ_N0[i]) meta_file_N0<<i-V<<" "<<meta_pot_N0[i]<<endl;
  //}
  
  //write Prob of occupying N
  //ofstream P_occ_N_file("/tmp/heat/P_occ_N.xmg");
  //std::vector<double> P_occ_N(2*V+1); 
  //double Z_N=0;
  //for(int i=0;i<=2*V;i++)
  //{
  //P_occ_N[i]=exp(meta_pot_N[i])-exp(-off_N);
  //Z_N+=P_occ_N[i];
  //}
  //for(int i=0;i<=2*V;i++)
  //{
  //P_occ_N[i]/=Z_N;
  //if(P_occ_N[i]) P_occ_N_file<<i-V<<" "<<P_occ_N[i]<<endl;
  //}
  
  /*
  ofstream P_occ_N_file_re("/tmp/heat/P_occ_N_re.xmg");
  std::vector<double> P_occ_N_re(2*V+1); 
  double mu=1.6;
  double pref_N=sqrt(3)*kappa*sinh(mu);
  double Z_re=0;
  for(int i=0;i<=2*V;i++)
    {
      P_occ_re[i]=(exp(meta_pot_N[i])-exp(-off_N))*cos(pref*(i-V));
      Z_re+=P_occ_re[i];
    }
  for(int i=0;i<=2*V;i++)
    {
      P_occ_re[i]/=Z_re;
      if(P_occ_re[i]) N_file_re<<i-V<<" "<<P_occ_re[i]<<endl;
    }
  */
  
  close();
  
  return 0;
}
