#include "include.h"

int main(int narg,char **arg)
{
  if(narg<4)
    {
      crash("Use: %s file njacks ientr raw[=false]",arg[0]);
      exit(1);
    }
  
  int njack=atoi(arg[2]);
  jack a(njack);
  
  a.load(arg[1],atoi(arg[3]));
  
  if(narg<5 || atoi(arg[4])==0) cout<<smart_print(a)<<endl;
  else cout<<raw_print(a)<<endl;
  
  return 0;
}
