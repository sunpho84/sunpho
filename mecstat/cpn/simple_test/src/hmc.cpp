#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "geometry.hpp"
#include "parameters.hpp"

//perform a hybid monte carlo update
void hmc_update()
{
  //allocate momenta
  dcomplex *pi=new dcomplex[V*N];
  double *omega=new double[V*NDIMS];
  
  //draw momenta
  for(int site=0;site<V;site++)
    {
      
    }
  
  delete[] pi;
  delete[] omega;
}
