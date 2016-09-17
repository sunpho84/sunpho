#include <math.h>
#include <iostream>

int L=24;

using namespace std;

double sqr(double x)
{return x*x;}

double relativistic_energy(double M,double th0)
{
  double qi=M_PI/L*th0;
  //return sqrt(M*M+3*qi*qi);
  return 2*asinh(sqrt(3*sqr(sin(qi/2))+sqr(sinh(M/2))));
}

int main()
{
  cout<<relativistic_energy(1.2028,1.26);
}
