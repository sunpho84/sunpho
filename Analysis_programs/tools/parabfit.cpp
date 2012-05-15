#include <include.h>
#include <iostream>

using namespace std;

int count(char *path)
{
  int n=0;
  double dum;
  
  ifstream in(path);
  while(in>>dum>>dum>>dum) n++;
  
  return n;
}

int main(int narg,char **arg)
{
  if(narg<2) crash("use %s filein [fileout]",arg[0]);
  
  int n=count(arg[1]);
  cout<<n<<endl;
  
  int nboot=1000;
  int njack=100;
  bvec Y(n,nboot,njack);
  double x[n],xmin,xmax;
  
  ifstream in(arg[1]);
  for(int i=0;i<n;i++)
    {
      double t,dt;
      in>>x[i]>>t>>dt;

      if(i==0) xmin=xmax=x[0];
      else
	{
	  xmin=min(xmin,x[i]);
	  xmax=max(xmax,x[i]);
	}

      Y[i].fill_gauss(t,dt,i+4);
      
      cout<<x[i]<<" "<<Y[i]<<endl;
    }
  
  bvec pars_fit=poly_fit(x,Y,2,xmin,xmax);
  cout<<pars_fit<<endl;
  
  if(narg>2)
    {
      ofstream out(arg[2]);
      out<<"@type xy"<<endl<<
	write_poly_with_error(pars_fit,0,xmax*1.1)<<endl<<
	"&"<<endl<<
	"@type xydy"<<endl;
	
      for(int i=0;i<n;i++)
	out<<x[i]<<" "<<Y[i]<<endl;
    }

  return 0;
}
