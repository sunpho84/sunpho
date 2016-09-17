#include "include.h"

int main(int narg,char **arg)
{
  if(narg<5) crash("use %s path njacks icorr nel [eff_mass]",arg[0]);
  
  char *path=arg[1];
  int njack=atoi(arg[2]);
  int icorr=atoi(arg[3]);
  int nel=atoi(arg[4]);
  
  jvec a(nel,njack);
  
  a.load(path,icorr);
  cout<<a<<endl;
  
  if(narg==6)
    {
      ofstream ef(arg[5]);
      ef<<effective_mass(a.simmetrized(1))<<endl;
    }
  
  return 0;
}
