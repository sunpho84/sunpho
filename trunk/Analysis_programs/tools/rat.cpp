#include <include.h>
#include <iostream>

using namespace std;

int main(int narg,char **arg)
{
  int nboot=100,njack=16;
  boot A(nboot,njack),B(nboot,njack),C(nboot,njack);

  A.load(arg[1]);
  B.load(arg[2]);
  C.load(arg[3]);
  
  cout<<A/B/C<<endl;

  return 0;
}
