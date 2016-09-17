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
  if(narg<3) crash("use %s filein degree [fileout] [min] [max]",arg[0]);
  
  int n=count(arg[1]);
  int deg=atoi(arg[2]);
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
  
  if(narg>4) sscanf(arg[4],"%lg",&xmin);
  if(narg>5) sscanf(arg[5],"%lg",&xmax);
  
  //fit
  switch(deg)
    {
    case 0: cout<<"Fitting with: a"<<endl;break;
    case 1: cout<<"Fitting with: a+b*x"<<endl;break;
    default: cout<<"Fitting with: a+b*x";
      for(int i=2;i<=deg;i++) cout<<"+"<<(char)('a'+i)<<"*x^"<<i;
      cout<<endl;
    }
  
  bvec pars_fit=poly_fit(x,Y,deg,xmin,xmax);
  for(int ipar=0;ipar<=deg;ipar++)
    cout<<(char)('a'+ipar)<<": "<<smart_print(pars_fit[ipar])<<endl;
  
  //compute chi2
  cout<<"Chi2: "<<smart_print(chi2_poly_fit(x,Y,deg,xmin,xmax,pars_fit))<<endl;
  
  //plot
  if(narg>3)
    {
      ofstream out(arg[3]);
      out<<"@type xy"<<endl<<
	write_poly_with_error(pars_fit,0,xmax*1.1)<<endl<<
	"&"<<endl<<
	"@type xydy"<<endl;
	
      for(int i=0;i<n;i++)
	out<<x[i]<<" "<<Y[i]<<endl;
    }

  return 0;
}
