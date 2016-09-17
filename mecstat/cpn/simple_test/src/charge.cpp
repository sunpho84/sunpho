#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <array>
#include <iostream>

#include "action.hpp"
#include "charge.hpp"
#include "data.hpp"
#include "geometry.hpp"
#include "lambda.hpp"
#include "zeta.hpp"

meta_pars_t chrono_charge("charge_potential");

//compute the charge (real and imaginary part)
dcomplex charge(dcomplex *z,dcomplex *l,int n,int mu)
{
  double resr=0,resi=0;
#pragma omp parallel for reduction(+:resr,resi)
  for(int site=0;site<V;site++)
    {
      int site_up=neighup(site,mu);
      dcomplex res=2.0*(conj(z[site_up*N+n])*z[site*N+n]*l[site*NDIMS+mu]);
      resr+=res.real();
      resi+=res.imag();
    }
  
  return dcomplex(resr,resi);
}

//compute the chargedynamical potential derivative using past history
double compute_meta_charge_pot_der(dcomplex *z,dcomplex *l)
{return chrono_charge.compute_pot_der(charge(z,l).imag());}

//compute the chargedynamical potential using past history
double compute_charge_pot(dcomplex *z,dcomplex *l)
{return chrono_charge.compute_pot(charge(z,l).imag());}

//compute the two pieces
void get_ch_pot_meta_pot(double &pot,double &meta_pot,dcomplex *z,dcomplex *l)
{
  if(use_charge_pot==2) meta_pot=compute_meta_charge_pot_der(z,l);
  else meta_pot=0;
  
  pot=-(cosh(ch_pot/L)-1)/g;
}

//compute the weights for staple fw and bw
void get_ch_staple_weight(dcomplex &w1,dcomplex &w2,dcomplex *z,dcomplex *l)
{
  double pot,meta_pot;
  get_ch_pot_meta_pot(pot,meta_pot,z,l);
  w1=2.0*dcomplex(-pot,-meta_pot);
  w2=2.0*dcomplex(-pot,+meta_pot);
}

//compute lambda force due to charge contribution
void sum_charge_lambda_force(dcomplex *z,dcomplex *l)
{
  double pot,meta_pot;
  get_ch_pot_meta_pot(pot,meta_pot,z,l);
  dcomplex w={-meta_pot,pot};
  
#pragma omp parallel for
  for(int site=0;site<V;site++)
    {
      int mu0=ch_pot_dir;
      int n0=ch_pot_n;
#ifdef DEBUG_HMC
      for(int mu=0;mu<NDIMS;mu++) fomega[site*NDIMS+mu]=0;
      fomega[site*NDIMS+mu0]=
#else
      fomega[site*NDIMS+mu0]+=
#endif
	2*get_lambda_real_scalprod(w*conj(l[site*NDIMS+mu0])*z[neighup(site,mu0)*N+n0],z[site*N+n0]);
    }
  
#ifdef DEBUG_HMC
  int site=0;
  double eps=1.e-6;
  double pre_act=charge_action(z,l);
  for(int mu=0;mu<NDIMS;mu++)
    {
      dcomplex pre=l[site*NDIMS+mu];
      l[site*NDIMS+mu]*=dcomplex(cos(eps),sin(eps));
      
      double post_act=charge_action(z,l);
      l[site*NDIMS+mu]=pre;
      
      double f=-(post_act-pre_act)/eps;
      std::cout<<"charge_lambda_force mu: "<<mu<<" fnu: -("<<post_act<<"-"<<pre_act<<")/"<<eps<<"="<<-(post_act-pre_act)<<"/"<<eps<<"="<<
	f<<" fan: "<<fomega[site*NDIMS+mu]<<" num: "<<f<<std::endl;
    }
#endif
}

//compute zeta force due to charge contribution
void sum_charge_zeta_force(dcomplex *z,dcomplex *l)
{
  //zeta orthogonalized spin projection
  int n0=ch_pot_n;
  int mu0=ch_pot_dir;
  dcomplex w1,w2;
  get_ch_staple_weight(w1,w2,z,l);
  
#pragma omp parallel for
  for(int site=0;site<V;site++)
    {
      //put staples to 0
      dcomplex s[N];
      for(int n=0;n<N;n++) s[n]=0;
      
      int site_up=neighup(site,mu0);
      int site_dw=neighdw(site,mu0);
      s[n0]+=w1*z[site_up*N+n0]*conj(l[site*NDIMS+mu0])+w2*z[site_dw*N+n0]*l[site_dw*NDIMS+mu0];
      
      zeta_orthogonalize_with(s,zeta+site*N);
      
#ifdef DEBUG_HMC
      for(int n=0;n<N;n++) fpi[site*N+n]=
#else
      for(int n=0;n<N;n++) fpi[site*N+n]+=
#endif
			     s[n];
    }
  
#ifdef DEBUG_HMC
  int site=0;
  double eps=1.e-8;
  double pre_act=charge_action(zeta,lambda);
  dcomplex pre_val[N];
  for(int m=0;m<N;m++) pre_val[m]=zeta[site*N+m];
  for(int n=0;n<N;n++)
    {
      for(int m=0;m<N;m++) zeta[site*N+m]=pre_val[m];
      zeta[site*N+n].imag(pre_val[n].imag()+eps);
      zeta_unitarize(zeta+site*N);
      double post_act=charge_action(zeta,lambda);
      
      double f=-(post_act-pre_act)/eps;
      
      cout<<"charge_zeta_force n: "<<n<<" fan: "<<fpi[site*N+n].imag()<<" num: -("<<post_act<<"-"<<pre_act<<")/"<<eps<<"="<<f<<endl;
    }
  for(int m=0;m<N;m++) zeta[site*N+m]=pre_val[m];
#endif

}
