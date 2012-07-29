#include <include.h>
#include <iostream>
#include <sstream>

using namespace std;

int T=48;
int njack=16;

int main(int narg,char **arg)
{
  jvec evn=jvec_load("/tmp/bub/bubble_evn_P5",T,njack,1);
  jvec odd=jvec_load("/tmp/bub/bubble_odd_P5",T,njack,1);
  jvec prd=jvec_load("/tmp/bub/bubble_prd_P5",T,njack,1);

  jvec tra=evn*0;
  for(int dt=0;dt<T;dt++)
    {
      tra[dt]=0;
      for(int t1=0;t1<T;t1++)
	{
	  int t2=(t1+dt)%T;
	  tra[dt]+=evn[t1]*odd[t2];
	}
      tra[dt]/=T;
    }
	
  cout<<(prd-tra).simmetrized(1)<<endl;
  
  return 0;
}
