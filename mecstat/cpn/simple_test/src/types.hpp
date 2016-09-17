#ifndef _TYPES_HPP
#define _TYPES_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <complex>
#include <fstream>
//#include <iostream>
#include <vector>

#include "macros.hpp"
#include "routines.hpp"
#include "tools.hpp"

using namespace std;

typedef complex<double> dcomplex;
typedef int coords[NDIMS];

namespace meta_pars
{
  const int symm=1;
  const int not_symm=0;
}

class meta_pars_t
{
public:
  int after;
  int each;
  double coeff;
  double width;
  double barr_left;
  double barr_right;
  int force_out_is_const;
  double force_out_left;
  double force_out_right;
  double well_tempering;
  double bend;
  
  int ngrid;
  vector<double> grid;
  
  string path;
  meta_pars_t(const char *in_path)
  {path=in_path;}
  
  //update the history-dependent potential
  void update(int isweep,double x)
  {
    if(isweep>=after && (isweep-after)%each==0)
      {
	int igrid=floor((x-barr_left)/width);
	double alpha=x/width;
	alpha=alpha-floor(alpha);
	if(igrid>=0 && igrid<=ngrid) grid[igrid]+=(1-alpha)*coeff;
	if(igrid+1>=0 && igrid+1<=ngrid) grid[igrid+1]+=alpha*coeff;
	//cout<<"Adding x="<<x<<", igrid="<<igrid<<", corresponding to x'="<<width*igrid+barr_left<<", alpha="<<alpha<<endl;
      }
  }
  
  //compute the derivative of the potential
  double compute_pot_der(double x)
  {
    //take igrid
    int igrid=floor((x-barr_left)/width);
    //cout<<"  pot_der for x="<<x<<", igrid="<<igrid<<", corresponding to x'="<<width*igrid+barr_left<<endl;
    //inside the barriers
    if(igrid>=0 && igrid<ngrid)
      return (grid[igrid+1]-grid[igrid])/width;
    else
      if(igrid<0)
	if(force_out_is_const) return -force_out_left;
	else return force_out_left*(x-barr_left);
      else
	if(force_out_is_const) return force_out_right;
	else return force_out_right*(x-barr_right);
  }
  
  //compute the potential using past history
  double compute_pot(double x)
  {
    //take igrid
    int igrid=floor((x-barr_left)/width);
    //cout<<" pot for x="<<x<<", igrid="<<igrid<<", corresponding to x'="<<width*igrid+barr_left<<endl;
    
    //inside the barriers
    if(igrid>=0 && igrid<ngrid)
      {
	//interpolate
	double x0=igrid*width+barr_left;
	double m=(grid[igrid+1]-grid[igrid])/width;
	double q=grid[igrid]-m*x0;
	return q+m*x;
      }
    else
      if(igrid<0)
	if(force_out_is_const) return force_out_left*fabs(x-barr_left)+grid[0];
	else return force_out_left*sqr(x-barr_left)/2+grid[0];
      else
	if(force_out_is_const) return force_out_right*fabs(x-barr_right)+grid[ngrid];
	else return force_out_right*sqr(x-barr_right)/2+grid[ngrid];
  }
  
  //write
  void save()
  {
    ofstream fout(path);
    fout.precision(16);
    for(int i=0;i<=ngrid;i++) fout<<barr_left+i*width<<" "<<grid[i]<<endl;
    fout.close();
  }
  
  //read
  void load()
  {
    ifstream fin(path);
    if(!fin.good()) crash("opening \"%s\"",path.c_str());
    for(int igrid=0;igrid<=ngrid;igrid++)
      {
	double xread;
	fin>>xread>>grid[igrid];
	if(!fin.good()) crash("reading line %d of \"%s\"",igrid,path.c_str());
	int jgrid=floor((xread-barr_left+width/2)/width);
	if(igrid!=jgrid) crash("found %d (%lg) when expecting %d",jgrid,xread,igrid);
      }
    fin.close();
  }
  
  //draw the chronological force
  void draw_force(const char *force_path)
  {
    double x_min=barr_left-(barr_right-barr_left)*.1124124;
    double x_max=barr_right+(barr_right-barr_left)*.1124124;
    double x_diff=x_max-x_min;
    int n=ceil(x_diff/width*10);
    if(n==0) n=1;
    double dx=x_diff/n;
    
    //compute
    double *xy=new double[n+1];
    double *xz=new double[n+1];
#pragma omp parallel for
    for(int i=0;i<=n;i++)
      {
	xy[i]=compute_pot_der(x_min+i*dx);
	xz[i]=(compute_pot(x_min+i*dx+dx/10)-compute_pot(x_min+i*dx-dx/10))/(dx/5);
      }
    
    //write
    ofstream fout(force_path);
    for(int i=0;i<=n;i++) fout<<x_min+i*dx<<" "<<xy[i]<<endl;
    fout<<"&"<<endl;
    for(int i=0;i<=n;i++) fout<<x_min+i*dx<<" "<<xz[i]<<endl;
    fout.close();
    
    delete[] xy;
    delete[] xz;
  }
  
  //initialize
  void init()
  {
    ngrid=(barr_right-barr_left+width/2)/width;
    grid.resize(ngrid+1);
    for(auto &g : grid) g=0;
  }
  
  //read all the parameters
  void read_pars(ifstream &input,const char *tag)
  {
    read(after,input,combine("Chrono%sAfter",tag));
    read(each,input,combine("Chrono%sEach",tag));
    read(coeff,input,combine("Chrono%sCoeff",tag));
    read(width,input,combine("Chrono%sWidth",tag));
    read(barr_left,input,combine("Chrono%sBarrLeft",tag));
    read(barr_right,input,combine("Chrono%sBarrRight",tag));
    read(force_out_left,input,combine("Chrono%sForceOutLeft",tag));
    read(force_out_right,input,combine("Chrono%sForceOutRight",tag));
    read(force_out_is_const,input,combine("Chrono%sForceOutIsConst",tag));
    read(bend,input,combine("Chrono%sBend",tag));
    read(well_tempering,input,combine("Chrono%sWellTempering",tag));
  }
};

#endif
