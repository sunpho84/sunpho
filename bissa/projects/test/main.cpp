#include "bissa.hpp"

using namespace bissa;

void in_main(int narg,char **arg)
{
}

int main(int narg,char **arg)
{
  init_bissa_threaded(narg,arg,in_main);
  close_bissa();
  
  return 0;
}
