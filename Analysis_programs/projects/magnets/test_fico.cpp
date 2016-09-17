#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void cart_to_pol(double &r,double &th,double &ph,double x,double y,double z)
{
  r=sqrt(x*x+y*y+z*z);
  th=atan2(y,x);
  ph=asin(z/r);
}

void pol_to_cart(double &x,double &y,double &z,double r,double th,double ph)
{
  x=r*cos(ph)*cos(th);
  y=r*cos(ph)*sin(th);
  z=r*sin(ph);
}

void print_arr(double x1,double y1,double z1,double x2,double y2,double z2)
{
  printf("%lg %lg %lg\n",x1,y1,z1);
  printf("%lg %lg %lg\n",x2,y2,z2);
  printf("\n\n");
}

int main()
{
  const int nth=13,nph=9;
  
  for(double u=0;u<2*M_PI;u+=2*M_PI/nth)
    for(double v=0;v<2*M_PI;v+=2*M_PI/nph)
      {
	double x1=cos(u)*(1.0*cos(v)+3),y1=sin(u)*(1.0*cos(v)+3),z1=1.0*sin(v);
	double x2=cos(u)*(1.3*cos(v)+3),y2=sin(u)*(1.3*cos(v)+3),z2=1.3*sin(v);
	
	for(int i=0;i<10;i++)
	  {
	    print_arr(x1,y1,z1,x2,y2,z2);
	    
	    double r,th,ph;
	    cart_to_pol(r,th,ph,x2-x1,y2-y1,z2-z1);
	    
	    double x3,y3,z3;
	    pol_to_cart(x3,y3,z3,r,th+((double)rand()/RAND_MAX-0.5)/3,ph+((double)rand()/RAND_MAX-0.5)/3);
	    x1=x2;
	    y1=y2;
	    z1=z2;
	    x2+=x3;
	    y2+=y3;
	    z2+=z3;
	  }
	
	//printf("set arrow %d from %lg,%lg,%lg to %lg,%lg,%lg\n",iarr++,x1,y1,z1,x2,y2,z2);
      }
    
  return 0;
}
