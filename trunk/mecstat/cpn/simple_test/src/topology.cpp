#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "data.hpp"
#include "geometry.hpp"
#include "parameters.hpp"
#include "types.hpp"
#include "zeta.hpp"

//return the geometric definition of topology
double geometric_topology_simplified()
{
  double topo=0;
  for(int n=0;n<V;n++)
    {
      int mu=0,nu=1;
      int nmu=neighup(n,mu);
      int nnu=neighup(n,nu);
      int nmu_nu=neighup(nmu,nu);
      
      topo+=
	arg(get_zeta_scalprod(zeta(nmu_nu),zeta(n))*
	    get_zeta_scalprod(zeta(nmu),zeta(nmu_nu))*
	    get_zeta_scalprod(zeta(n),zeta(nmu)))+
	arg(get_zeta_scalprod(zeta(nnu),zeta(n))*
	    get_zeta_scalprod(zeta(nmu_nu),zeta(nnu))*
	    get_zeta_scalprod(zeta(n),zeta(nmu_nu)));
    }
  
  return topo/(2*M_PI);
}

//alternative edition
double geometric_topology()
{
  int mu=0,nu=1;
  double topo=0;
  for(int n=0;n<V;n++)
    {
      dcomplex P1[N*N],P2[N*N],P3[N*N];
      get_zeta_P(P1,zeta(n));
      get_zeta_P(P3,zeta(neighup(neighup(n,mu),nu)));

      dcomplex c;
      
      c=0;
      get_zeta_P(P2,zeta(neighup(n,mu)));
      for(int i=0;i<N;i++)
	for(int j=0;j<N;j++)
	  for(int k=0;k<N;k++)
	    c+=P3[i*N+j]*P2[j*N+k]*P1[k*N+i];
      topo+=log(c).imag();
      
      c=0;
      get_zeta_P(P2,zeta(neighup(n,nu)));
      for(int i=0;i<N;i++)
	for(int j=0;j<N;j++)
	  for(int k=0;k<N;k++)
	    c+=P2[i*N+j]*P3[j*N+k]*P1[k*N+i];
      topo+=log(c).imag();
    }
  
  return topo/(2*M_PI);
}

//gauge version
double topology()
{
  int mu=0,nu=1;
  double topo=0;
  for(int n=0;n<V;n++)
    topo+=(lambda(n)[mu]*lambda(neighup(n,mu))[nu]*conj(lambda(neighup(n,nu))[mu]*lambda(n)[nu])).imag();
  
  return topo;
}
