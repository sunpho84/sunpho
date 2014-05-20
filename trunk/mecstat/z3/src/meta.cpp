#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "action.hpp"
#include "geometry.hpp"
#include "macros.hpp"
#include "parameters.hpp"

#include <iostream>
#include <cmath>

using namespace std;

//increase the meta-potential relative to a certain occupation
double compute_meta_potential(int occ_N,int occ_N0,int itraj)
{
  double P=0;
  
  if(meta_coeff!=0)
    {
      double pref=meta_coeff/(meta_sigma_N*meta_sigma_N0*sqrt(2*M_PI));
      for(int jtraj=0;jtraj<itraj;jtraj++)
	{
	  double dN=(occ_N-data_N[jtraj].N)/meta_sigma_N;
	  double dN0=(occ_N0-data_N[jtraj].N0)/meta_sigma_N0;
	  
	  double d2=dN*dN+dN0*dN0;
	  if(d2<9) P+=exp(-d2/2);
	}
      P*=pref;
    }
  
  return P;
}
