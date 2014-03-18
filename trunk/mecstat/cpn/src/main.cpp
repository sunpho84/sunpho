#include <iostream>

#include "geometry.hpp"
#include "simul.hpp"

using namespace std;

//internal main
void in_main(int narg,char **arg)
{
  GET_THREAD_ID();
  
  //test allocate and deallocate  
  double *v=ALLOCATE("v",20,double);
  FREE(v);
  
  geometry_t *geometry=CAST_PTR_FROM_MASTER_THREAD(new geometry_t(10));
  
  if(IS_MASTER_THREAD)
    {
      //
    }
  THREAD_BARRIER();

}

int main(int narg,char**arg)
{
  simul=new simul_t(narg,arg,in_main);
  
  return 0;
}
