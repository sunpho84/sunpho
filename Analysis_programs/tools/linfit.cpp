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

      //find min max for graph
      if(i==0) xmin=xmax=x[0];
      else
	{
	  xmin=min(xmin,x[i]);
	  xmax=max(xmax,x[i]);
	}
      
      //put data
      Y[i].fill_gauss(t,dt,i+4);
      
      cout<<"Data "<<i<<": "<<x[i]<<" "<<smart_print(Y[i])<<endl;
    }
  
  //fit
  cout<<"Fitting with: a+b*x"<<endl;
  bvec pars_fit=poly_fit(x,Y,1,xmin,xmax);
  cout<<"a: "<<smart_print(pars_fit[0])<<endl;
  cout<<"b: "<<smart_print(pars_fit[1])<<endl;
  
  //plot
  if(narg>2)
    {
      ofstream out(arg[2]);
      out<<"@type xy"<<endl<<
	write_poly_with_error(pars_fit,0,xmax*1.1)<<endl
	 <<"&"<<endl
	 <<"@s1 line color 1"<<endl
	 <<"@type xydy"<<endl
	 <<"@s2 linestyle 0"<<endl
	 <<"@s2 symbol 1"<<endl
	 <<"@s2 symbol color 2"<<endl
	 <<"@s2 errorbar color 2"<<endl
	;
	
      for(int i=0;i<n;i++)
	out<<x[i]<<" "<<Y[i]<<endl;
    }

  return 0;
}
