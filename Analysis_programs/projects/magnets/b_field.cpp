#include <complex>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>

using namespace std;

typedef complex<double> cmplex;
typedef cmplex bi_u1[2];
cmplex I(0,1);

const int L=4;

int nup(int x)
{
  if(x==L-1) return 0;
  else return x+1;
}

int ndw(int x)
{
  if(x==0) return L-1;
  else return x-1;
}

void ph(bi_u1 &u,int x,int y,double b)
{
  const double ar=2*b*M_PI/L/L;
  if(x!=L-1) u[0]=1;
  else       u[0]=exp(-I*ar*(double)(L*y));
  u[1]=exp(I*ar*(double)x);
}

void plot(double b)
{
  bi_u1 u[L*L];
  
  for(int x=0;x<L;x++)
    for(int y=0;y<L;y++)
      ph(u[y*L+x],x,y,b);
  
  char bstr[10];
  sprintf(bstr,"%1.5f",b);
  
  ofstream out("out_video.gp");
  out<<"set terminal pdf fontscale 2\n"
    "set output\"out_"<<bstr<<".pdf\"\n"
    "set pm3d explicit at s depthorder hidden3d 1\n"
    "set hidden3d front\n"
    "set style fill transparent solid 0.65\n"
    "set palette rgb 9,9,3\n"
    "unset border\n"
    "unset tics\n"
    "unset key\n"
    "unset colorbox\n"
    "set xrange [0:"<<L<<"]\n"
    "set yrange [0:"<<L<<"]\n"
    "set zrange [0:1.6]\n"
    "set isosamples 5,5\n"
    "set style line 1 lc rgb '#157545' lt 1 lw 0 # --- green\n"
    "set style line 2 lc rgb '#157545' lt 1 lw 2 # --- green\n"
    "set view 60,30,1.2,1\n"
    "set title \"b="<<bstr<<"\"\n";
  int iarr=1;
  for(int y=L-1;y>=0;y--)
    for(int x=0;x<L;x++)
      {
	int nx=nup(x);
	int ny=nup(y);
	
	cmplex a=u[y*L+x][0]*u[y*L+nx][1]*conj(u[ny*L+x][0])*conj(u[y*L+x][1]);
	
	out<<"set arrow "<<iarr++<<" from "<<x+0.5<<","<<y+0.5<<",0 to "<<x+0.5<<","<<y+0.5<<","<<1-a.real()<<
	  "  lc rgb \"red\" lw 1 head front "<<endl;
      }
  out<<"splot 0 w pm3d"<<endl;
  
  system("gnuplot out_video.gp");
}

int main()
{
  for(double b=4;b<=5;b+=0.01)
    plot(b);
  return 0;
}
