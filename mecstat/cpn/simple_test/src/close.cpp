#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "data.hpp"
#include "geometry.hpp"
#include "random.hpp"
#include "staples.hpp"

//close the code
void close()
{
  delete[] lambda;
  delete[] zeta;
  
  delete[] lambda_old;
  delete[] zeta_old;
  
  delete[] fpi;
  delete[] fomega;
  
  delete[] pi;
  delete[] omega;
  
  delete[] topo_staples_supp_data;
  delete[] topo_staples_data;
  
  for(int istout_lev=1;istout_lev<=nstout_lev;istout_lev++) delete[] lambda_stout[istout_lev];
  delete[] lambda_stout;
  
  delete[] neigh_data;
}
