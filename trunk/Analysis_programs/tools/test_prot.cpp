#include <include.h>
#include <iostream>

using namespace std;

int main()
{
  int T=48;
  int njack=15;
  
  jvec a=jvec_load("/tmp/ave",T,njack,0);
  jvec b=-jvec_load("/tmp/ave",T,njack,1).simmetric();
  
  a.print_to_file("/tmp/o");
  b.print_to_file("/tmp/p");
  
  return 0;
}
