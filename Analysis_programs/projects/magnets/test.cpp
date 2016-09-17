#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void cart_to_pol(double &r,double &th,double &ph,double x,double y,double z)
{
  r=sqrt(x*x+y*y+z*z);
  ph=asin(z/r);
  th=atan2(y,x);
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
  const int nth=14,nph=14;
  const double r=1;
  const double l=0.6;
  
  for(double u=-M_PI+M_PI/nth;u<M_PI;u+=2*M_PI/nth)
  for(double v=-M_PI+M_PI/nph;v<M_PI;v+=2*M_PI/nph)
    //for(int narr=0;narr<400;narr++)
    {
      //double u=M_PI*(2.0*rand()/RAND_MAX-1);
      //double v=M_PI*(2.0*rand()/RAND_MAX-1);
      if(v>-M_PI/3)
	{
	  double x1=cos(u)*(r*cos(v)+3),y1=sin(u)*(r*cos(v)+3),z1=r*sin(v);
	  double x2=cos(u)*((r+l)*cos(v)+3),y2=sin(u)*((r+l)*cos(v)+3),z2=(r+l)*sin(v);
	  print_arr(x1,y1,z1,x2,y2,z2);
	  
	  double r,th,ph;
	  cart_to_pol(r,th,ph,x1-x2,y1-y2,z1-z2);
	  
	  for(int i=0;i<2;i++)
	    {
	      double x3,y3,z3;
	      pol_to_cart(x3,y3,z3,r/4,th+i*M_PI/6,ph+(1-i)*M_PI/6);
	      print_arr(x2,y2,z2,x2+x3,y2+y3,z2+z3);
	      
	      double x4,y4,z4;
	      pol_to_cart(x4,y4,z4,r/4,th-i*M_PI/6,ph-(1-i)*M_PI/6);
	      print_arr(x2,y2,z2,x2+x4,y2+y4,z2+z4);
	    }
	  //printf("set arrow %d from %lg,%lg,%lg to %lg,%lg,%lg\n",iarr++,x1,y1,z1,x2,y2,z2);
	}
    }
    
  return 0;
}
