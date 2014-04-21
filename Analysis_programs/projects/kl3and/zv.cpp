#include "include.h"

const int T=64,njacks=38;
const int atw=1;

class results_t
{
public:
  jack Zv_pion;
  jack Zv_kaon_a;
  jack Zv_kaon_b;
};

int icombo(int tsep,int iv)
{return iv+3*tsep;}

jvec simmetric(jvec in)
{
  jvec out(T,njacks);
  for(int t=0;t<T;t++) out[(T-t)%T]=in[t];

  return out;
}

jvec load(const char *path,int icombo)
{
  jvec out(T,njacks);
  out.load(combine("data/total_%s",path).c_str(),icombo);
  
  return out;
}

void simple_plot(const char *path,jvec in,int tsep=-1,ios::openmode mode=ios::out)
{
  ofstream out(combine("plots/%s.xmg",path).c_str(),mode);
  out<<"@type xydy"<<endl;
  if(tsep!=-1) out<<"#tsep"<<tsep<<endl;
  for(int i=0;i<in.nel;i++)
    {
      int t;
      if(tsep==-1) t=i;
      else t=i-tsep/2;
      out<<t<<" "<<in[i]<<endl;
    }
  out<<"&"<<endl;
}

results_t study_single_tsep(int tsep,ios::openmode mode=ios::out)
{
  //fit pion mass
  jvec pion_rest=load("pion-00WW",0).simmetrized(1),pion_rest_alt=load("pion-00WP",0).simmetrized(1);
  jack M_Pi=constant_fit(effective_mass(pion_rest_alt),11,T/2-1,"plots/pion_rest_alt_effmass.xmg");
  cout<<"Pion mass: "<<smart_print(M_Pi)<<endl;
  
  //fit kaon mass
  jvec kaon_rest=load("kaon-00WW",0).simmetrized(1),kaon_rest_alt=load("kaon-00WP",0).simmetrized(1);
  jack M_K=constant_fit(effective_mass(kaon_rest_alt),14,T/2-1,"plots/kaon_rest_alt_effmass.xmg");  
  cout<<"Kaon mass: "<<smart_print(M_K)<<endl;
  
  //plot the effective mass of less preicse corrs
  simple_plot("pion_rest_effmass",effective_mass(pion_rest));
  simple_plot("kaon_rest_effmass",effective_mass(kaon_rest));
  
  ////////////////////////////////// renormalization //////////////////////////////////
  
  //load insertion of V0 between pions or kaons at rest
  jvec pp_V0_rest=(load("zpa-00",icombo(tsep,1))-atw*simmetric(load("zpa-00",icombo(T-tsep,1))))/(atw+1);
  jvec kk_V0_rest_a=(load("zka-00",icombo(tsep,1))-atw*simmetric(load("zka-00",icombo(T-tsep,1))))/(atw+1);
  jvec kk_V0_rest_b=(load("zkb-00",icombo(tsep,1))-atw*simmetric(load("zkb-00",icombo(T-tsep,1))))/(atw+1);
  
  simple_plot(combine("pp_V0_rest%02d",tsep).c_str(),aperiodic_effective_mass(pp_V0_rest));
  //define tilded
  jvec pion_rest_tilde(T,njacks);
  jvec kaon_rest_tilde(T,njacks);
  for(int t=0;t<=T/2;t++)
    {
      pion_rest_tilde[t]=pion_rest[t]-pion_rest[T/2]*exp(-M_Pi*(T/2-t))/2;
      kaon_rest_tilde[t]=kaon_rest[t]-kaon_rest[T/2]*exp(-M_K*(T/2-t))/2;
    }
  for(int t=T/2;t<T;t++)
    {
      //pion_rest_tilde[t]=pion_rest[t]-pion_rest[T/2]*exp(-M_Pi*(T/2-t))/2;
      //kaon_rest_tilde[t]=kaon_rest[t]-kaon_rest[T/2]*exp(-M_K*(T/2-t))/2;
      pion_rest_tilde[t]=pion_rest[T/2]*exp(-M_Pi*(t-T/2))/2;
      kaon_rest_tilde[t]=kaon_rest[T/2]*exp(-M_K*(t-T/2))/2;
    }
  
  //plot the tilded
  simple_plot("pion_rest_tilde",pion_rest_tilde);
  simple_plot("kaon_rest_tilde",kaon_rest_tilde);

  //Zv
  jvec Zv_pion_rest_corr=-pion_rest_tilde[tsep]/pp_V0_rest;
  cout<<"1-Corr["<<tsep<<"]: "<<1+cov(pion_rest_tilde[tsep],pp_V0_rest[tsep/2])/sqrt(var(pion_rest_tilde[tsep])*var(pp_V0_rest[tsep/2]))<<endl;
  simple_plot("Zv_pion_rest",Zv_pion_rest_corr.subset(0,tsep)+0.005*tsep/4,tsep,mode);
  jvec Zv_kaon_rest_a_corr=-kaon_rest_tilde[tsep]/kk_V0_rest_a;
  simple_plot("Zv_kaon_rest_a",Zv_kaon_rest_a_corr.subset(0,tsep),tsep,mode);
  jvec Zv_kaon_rest_b_corr=-kaon_rest_tilde[tsep]/kk_V0_rest_b;
  simple_plot("Zv_kaon_rest_b",Zv_kaon_rest_b_corr.subset(0,tsep),tsep,mode);
  
  //extract Zv from pion and kaon_a at rest
  jack Zv_pion=constant_fit(Zv_pion_rest_corr.subset(0,tsep),tsep/2-tsep/4,tsep/2+tsep/4,combine("plots/Zv_pion%02d.xmg",tsep).c_str());
  jack Zv_kaon_a=constant_fit(Zv_kaon_rest_a_corr,tsep/2-tsep/4,tsep/2+tsep/4);
  jack Zv_kaon_b=constant_fit(Zv_kaon_rest_b_corr,tsep/2-tsep/4,tsep/2+tsep/4);
  cerr<<
    " tsep: "<<tsep<<
    " pi_rest: "<<smart_print(Zv_pion)<<
    " k_rest_a: "<<smart_print(Zv_kaon_a)<<
    " k_rest_b: "<<smart_print(Zv_kaon_b)<<endl;
  
  //save results
  results_t results;
  results.Zv_pion=Zv_pion;
  results.Zv_kaon_a=Zv_kaon_a;
  results.Zv_kaon_b=Zv_kaon_b;
  
  return results;
}

int main()
{
  std::map<int,results_t> results;
  for(int tsep=4;tsep<=40;tsep+=4) results[tsep]=study_single_tsep(tsep,ios::out|ios::app);
  
  ofstream results_out("plots/results.xmg");
  results_out<<"@type xydy"<<endl;
  for(auto it=results.begin();it!=results.end();it++)
    results_out<<it->first<<" "<</*smart_print*/(it->second.Zv_pion)<<endl;
  results_out<<"&"<<endl;
  for(auto it=results.begin();it!=results.end();it++)
    results_out<<it->first<<" "<</*smart_print*/(it->second.Zv_kaon_a)<<endl;
  results_out<<"&"<<endl;  
  for(auto it=results.begin();it!=results.end();it++)
    results_out<<it->first<<" "<</*smart_print*/(it->second.Zv_kaon_b)<<endl;

  for(auto it=results.begin();it!=results.end();it++)
    cout<<it->first<<
      "\t"<<smart_print(it->second.Zv_pion)<<
      "\t"<<smart_print(it->second.Zv_kaon_a)<<
      "\t"<<smart_print(it->second.Zv_kaon_b)<<endl;

  
  return 0;
}
