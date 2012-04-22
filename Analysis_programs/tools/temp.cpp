#include <include.h>
#include <iostream>

using namespace std;

int main()
{
  int i;
  double x,ex;
  ifstream in("/tmp/a");
  
  boot oo;
  
  int nboot=100,njack=20;
  for(int l=0;l<24;l++)
    {
      in>>i>>x>>ex;
      boot o(nboot,njack);
      o.fill_gauss(x,ex,l);
      
      o=-log(o);
      if(l>0)
	cout<<l-1<<" "<<o-oo<<endl;
      
      oo=o;
    }
  in.close();
  
  return 0;
}
