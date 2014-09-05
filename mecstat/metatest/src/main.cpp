#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <signal.h>
#include <omp.h>
#include "random_gen.hpp"

#define GAUSS

using namespace std;

int itraj;
int do_not_exit;

double A=1; //curvature of the pot
double B=0; //weight
double C=1; //place-dependant weight

const double SKIP_OUT_FROM=50;
int nacc=0;
const int debug=0;
const double chrono_barr=2.2;//8;
const double chrono_width=0.25;//W/4;
const double harm_weight=10;
const double chrono_coeff=0.01;
double *chrono_hist;
double *chrono_weight;

inline double sqr(double x)
{return x*x;}

//exit when called
void handle_signal(int sig)
{
  cout<<"ciao"<<endl;
  do_not_exit=0;
}

double x_action(double x)
{return A*x*x/2;}

double p_action(double p)
{return p*p/2;}

double meta_action(double X,bool ave=false)
{
  double chrono_hist_potential=0;
  double harm_pot=0;
  
  if(chrono_coeff!=0)
    {
      if(X<-chrono_barr) harm_pot=harm_weight*sqr(-X-chrono_barr)/2;
      if(X>+chrono_barr) harm_pot=harm_weight*sqr(+X-chrono_barr)/2;
      
      if(X<-chrono_barr) X=-chrono_barr;
      if(X>+chrono_barr) X=+chrono_barr;
      
      int end=itraj;
      int start=(ave?(end/2):0);
      int size=end-start;
#pragma omp parallel for reduction(+:chrono_hist_potential)
      for(int it=0;it<end;it++)
	{
	  int aw=(ave?((it>start)?(end-it):(end-start)):1);
	  double x=chrono_hist[it];
	  double diff=X-x,f=diff/chrono_width;
	  if(fabs(f)<SKIP_OUT_FROM)
	    {
#ifdef GAUSS
	      double cont=aw*exp(-f*f/2+C*X*X/2);
#else
	      double cont=aw/(1+f*f);
#endif
	      
	      chrono_hist_potential+=cont;
	    }
	}
      chrono_hist_potential*=chrono_coeff/(ave?size:1);
    }
  
  return chrono_hist_potential+harm_pot;
}

double action(double x,double p)
{return x_action(x)+p_action(p)+meta_action(x);}

double action_force(double x)
{return -A*x;}

double meta_force(double X)
{
  double chrono_hist_potential=0;
  double harm_pot=0;
  
  if(chrono_coeff!=0)
    {
      double pref=chrono_coeff/sqr(chrono_width);
      if(X>-chrono_barr && X<+chrono_barr)
	{
#pragma omp parallel for reduction(+:chrono_hist_potential)
	  for(int it=0;it<itraj;it++)
	    {
	      double x=chrono_hist[it];
	      double w=chrono_weight[it];
	      double diff=X-x,f=diff/chrono_width;
	      if(fabs(f)<SKIP_OUT_FROM)
		{
#ifdef GAUSS
		  double cont=(pref*diff-chrono_coeff*C*X)*w*exp(-f*f/2+C*X*X/2);
#else
		  double cont=2*diff*w/sqr(1+f*f);
#endif
		  chrono_hist_potential+=cont;
		}
	    }
	}
      //chrono_hist_potential*=pref;
      
      if(X<-chrono_barr) harm_pot=+harm_weight*(-X-chrono_barr);
      if(X>+chrono_barr) harm_pot=-harm_weight*(+X-chrono_barr);
    }
  
  return chrono_hist_potential+harm_pot;
}

double p_force(double x)
{return action_force(x)+meta_force(x);}

void update_x(double &x,double p,double dt)
{x+=p*dt;}

void update_p(double &p,double x,double dt)
{p+=p_force(x)*dt;}

void omelyan(double &x,double &p,double tl,int nt)
{ 
  const double OMELYAN_LAMBDA=0.1931833;
  double dt=tl/nt/2,dth=dt/2,ldt=dt*OMELYAN_LAMBDA,l2dt=2*OMELYAN_LAMBDA*dt,uml2dt=(1-2*OMELYAN_LAMBDA)*dt;

  update_p(p,x,ldt);
    
  for(int it=0;it<nt;it++)
    {
      //cerr<<x<<endl;
      update_x(x,p,dth);
      update_p(p,x,uml2dt);
      update_x(x,p,dth);
      //cerr<<x<<endl;
      update_p(p,x,(it==(nt-1))?ldt:l2dt);
    }
}

bool hmc(double &x,int nt)
{
  double p=get_gauss_double(0.0,1.0);
  double init_act=action(x,p);
  
  double x_old=x;
  
  double tl=1.0;
  omelyan(x,p,tl,nt);
  
  double final_act=action(x,p);
  double delta_act=final_act-init_act;
  double pacc=exp(-delta_act);
  double est=get_unif_double(0,1);
  
  if(est<pacc) 
    {
      if(debug) cout<<x<<" accepted ("<<delta_act<<")"<<endl;
      nacc++;
      
      return 1;
    }
  else
    {
      if(debug) cout<<x<<" rejected ("<<delta_act<<")"<<endl;
      x=x_old;
      
      return 0;
    }
}

int main()
{
  const int ntraj=100000;
  chrono_hist=new double[ntraj];
  chrono_weight=new double[ntraj];
  
  const int nt=10;
  double x=1;
  
  ofstream pot("pot");
  for(double e=-1.5*chrono_barr;e<chrono_barr*1.5;e+=0.02)
    pot<<e<<" "<<-4+x_action(e)<<endl;
  
  ofstream force("force");
  for(double e=-1.5*chrono_barr;e<chrono_barr*1.5;e+=0.02)
    force<<e<<" "<<action_force(e)<<endl;
  
  signal(SIGINT,handle_signal);
  itraj=0;
  do_not_exit=1;
  do
    {
      hmc(x,nt);
      chrono_hist[itraj]=x;
      chrono_weight[itraj]=exp(B*x*x/2);
      //cerr<<x<<endl;
      itraj++;
      if(itraj%(ntraj/100)==0) cout<<(itraj*100/ntraj)<<"%"<<endl;
    }
  while(itraj<ntraj && do_not_exit);
  
  ofstream meta_pot_file("meta_pot");
  double meta_ave_0=meta_action(0,true);
  for(double e=-1.5*chrono_barr;e<chrono_barr*1.5;e+=0.1)
    meta_pot_file<<e<<" "<<-(meta_action(e,true)-meta_ave_0)<<endl;

  ofstream meta_force_file("meta_force");
  for(double e=-1.5*chrono_barr;e<chrono_barr*1.5;e+=0.02)
    meta_force_file<<e<<" "<<meta_force(e)<<endl;
  meta_force_file<<"&"<<endl;
  double la=0.001;
  for(double e=-1.5*chrono_barr;e<chrono_barr*1.5;e+=0.02)
    meta_force_file<<e<<" "<<-((meta_action(e+la)-meta_action(e-la)))/la/2<<endl;
  
  ofstream meta_histo_file("meta_histo");
  for(int it=0;it<itraj;it++) meta_histo_file<<chrono_hist[it]<<endl;
  
  ofstream meta_weight_file("meta_weight");
  for(int it=0;it<itraj;it++) meta_weight_file<<chrono_weight[it]<<endl;
  
  cout<<"acceptance: "<<nacc/(double)itraj<<endl;  

  delete[] chrono_hist;
  delete[] chrono_weight;
  
  return 0;
}
