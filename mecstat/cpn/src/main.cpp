#include <iostream>

#include "simul.hpp"

using namespace std;

//internal main
void in_main(int narg,char **arg)
{
  printf("ciao from in_main!\n");
}

int main(int narg,char**arg)
{
  simul=new simul_t(narg,arg,in_main);
  
  //geometry_t geometry(10);
  
  return 0;
}
