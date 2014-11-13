#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "geometry.hpp"
#include "types.hpp"

//compute polyakov loop
dcomplex polyakov(dcomplex *lambda)
{
  dcomplex p=0;
  
  for(int x=0;x<L;x++)
    {
      dcomplex loc_p={1,0};
      for(int t=0;t<L;t++)
	{
	  coords c={x,t};
	  loc_p*=lambda[site_of_coords(c)*NDIMS+1];
	}
      p+=loc_p;
    }
  p/=L;
  
  return p;
}
