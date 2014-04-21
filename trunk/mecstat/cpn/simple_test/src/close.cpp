#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "data.hpp"
#include "geometry.hpp"
#include "random.hpp"

//close the code
void close()
{
  delete[] lambda_data;
  delete[] zeta_data;
  
  delete[] neigh_data;
  
#ifdef GOOD_GENERATOR
  delete dis;
  delete gen;
  delete rd;
#endif
}
