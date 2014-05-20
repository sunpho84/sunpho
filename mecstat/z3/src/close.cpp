#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "data.hpp"
#include "random.hpp"

//close the code
void close()
{
  delete[] phi;
  
#ifdef GOOD_GENERATOR
  delete dis;
  delete gen;
  delete rd;
#endif
}
