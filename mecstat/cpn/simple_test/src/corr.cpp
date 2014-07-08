#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "data.hpp"
#include "geometry.hpp"
#include "parameters.hpp"
#include "types.hpp"

#include <iostream>

using namespace std;

//compute correlation function
void compute_corr(double *out,dcomplex *z)
{
  int n=0;
  for(int x=0;x<L;x++)
    {
      double res=0;
#pragma omp parallel for reduction(+:res)
      for(int y=0;y<L;y++)
	for(int dx=0;dx<L;dx++)
	  {
	    coords c={(x+dx)%L,y};
	    coords c0={dx,0};
	    int s=site_of_coords(c);
	    int s0=site_of_coords(c0);
	    n++;
	    
	    dcomplex t=0.0;
	    for(int n=0;n<N;n++) t+=conj(z[s*N+n])*z[s0*N+n];
	    res+=norm(t);
	  }
      out[x]=res/L;
    }
  cout<<n<<endl;
}
