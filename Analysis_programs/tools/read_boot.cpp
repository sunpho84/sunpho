#include "include.h"

int main(int narg,char **arg)
{
  if(narg<5)
    {
      crash("Use: %s file nboot njacks ientr raw[=false]",arg[0]);
      exit(1);
    }
  
  int nboot=atoi(arg[2]);
  int njack=atoi(arg[3]);
  boot a(nboot,njack);
  
  a.load(arg[1],atoi(arg[4]));
  
  if(narg<6 || atoi(arg[5])==0) cout<<smart_print(a)<<endl;
  else cout<<raw_print(a)<<endl;
  
  return 0;
}
