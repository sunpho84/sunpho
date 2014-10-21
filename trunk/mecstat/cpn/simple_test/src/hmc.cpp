#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "action.hpp"
#include "data.hpp"
#include "geometry.hpp"
#include "lambda.hpp"
#include "parameters.hpp"
#include "random.hpp"
#include "routines.hpp"
#include "staples.hpp"
#include "topology.hpp"
#include "tools.hpp"
#include "zeta.hpp"

#include <iostream>

//draw momenta
void generate_momenta()
{
  for(int site=0;site<V;site++)
    {
      //generate the lambda momenta, gaussianly
      for(int mu=0;mu<NDIMS;mu++) omega[site*NDIMS+mu]=get_gauss_double();
      
      //generate zeta momenta
      for(int n=0;n<N;n++)
	{
	  pi[site*N+n].real(get_gauss_double());
	  pi[site*N+n].imag(get_gauss_double());
	}
      
      //orthogonalize
      zeta_orthogonalize_with(pi+site*N,zeta+site*N);
    }
}

//compute the action for zeta momenta
double zeta_momenta_action()
{
  double act=0;
#pragma omp parallel for reduction(+:act)
  for(int site=0;site<V;site++) act+=get_zeta_real_scalprod(pi+site*N,pi+site*N);

  return act/2;
}

//compute the action for lambda momenta
double lambda_momenta_action()
{
  double act=0;
#pragma omp parallel for reduction(+:act)
  for(int site=0;site<V;site++) for(int mu=0;mu<NDIMS;mu++) act+=sqr(omega[site*NDIMS+mu]);
  
  return act/2;
}

//compute the action of momenta
double momenta_action()
{return zeta_momenta_action()+lambda_momenta_action();}

//compute lambda force
void compute_lambda_forces()
{
#pragma omp parallel for
  for(int site=0;site<V;site++)
    for(int mu=0;mu<NDIMS;mu++)
      {
        //reset
        fomega[site*NDIMS+mu]=0;
        for(int n=0;n<N;n++)
          {
            dcomplex t=dcomplex(0,1)*conj(lambda[site*NDIMS+mu])*zeta[neighup(site,mu)*N+n];
            fomega[site*NDIMS+mu]+=get_lambda_real_scalprod(zeta[site*N+n],t);
          }
        fomega[site*NDIMS+mu]*=-2*beta*N;
      }
  
#ifdef DEBUG_HMC
  int site=0;
  double eps=1.e-6;
  double pre_act=action(zeta,lambda);
  for(int mu=0;mu<NDIMS;mu++)
    {
      double pre_val=arg(lambda[site*NDIMS+mu]);
      lambda[site*NDIMS+mu]=exp(dcomplex(0,1)*(pre_val+eps));
      
      double post_act=action(zeta,lambda);
      lambda[site*NDIMS+mu]=exp(dcomplex(0,1)*pre_val);
      
      double f=-(post_act-pre_act)/eps;
      cout<<"mu: "<<mu<<" f: "<<f<<" "<<f/fomega[site*NDIMS+mu]<<endl;
    }
#endif
}

//compute zeta forces
void compute_zeta_forces()
{
  //zeta orthogonalized spin projection
#pragma omp parallel for
  for(int site=0;site<V;site++)
    {
      site_staple(fpi+site*N,zeta,lambda,site);
      zeta_orthogonalize_with(fpi+site*N,zeta+site*N);
      for(int n=0;n<N;n++) fpi[site*N+n]*=beta*N;
    }
  
#ifdef DEBUG_HMC
  int site=0;
  double eps=1.e-8;
  double pre_act=action(zeta,lambda);
  dcomplex pre_val[N];
  for(int m=0;m<N;m++) pre_val[m]=zeta[site*N+m];
  for(int n=0;n<N;n++)
    {
      for(int m=0;m<N;m++) zeta[site*N+m]=pre_val[m];
      zeta[site*N+n].imag(pre_val[n].imag()+eps);
      zeta_unitarize(zeta+site*N);    
      double post_act=action(zeta,lambda);
      
      double f=-(post_act-pre_act)/eps;
      
      cout<<"n: "<<n<<" f: "<<f<<" "<<f/fpi[site*N+n].imag()<<endl;
    }
  for(int m=0;m<N;m++) zeta[site*N+m]=pre_val[m];
#endif
}

//update pi momenta
void update_zeta_momenta(double eps)
{
  //compute the zeta forces
  compute_zeta_forces();

  //update zeta momenta
#pragma omp parallel for
  for(int site=0;site<V;site++) for(int n=0;n<N;n++) pi[site*N+n]+=fpi[site*N+n]*eps;
}

//update omega momenta
void update_lambda_momenta(double eps)
{
  //compute the lambda forces
  compute_lambda_forces();

  //update momenta
#pragma omp parallel for
  for(int site=0;site<V;site++) for(int mu=0;mu<NDIMS;mu++) omega[site*NDIMS+mu]+=fomega[site*NDIMS+mu]*eps;
  
  if(use_topo_pot)
    {
      //compute topological lambda forces
      compute_topological_force(fomega,stout_rho,nstout_lev,lambda);
      //update momenta
#pragma omp parallel for
      for(int site=0;site<V;site++) for(int mu=0;mu<NDIMS;mu++) omega[site*NDIMS+mu]+=fomega[site*NDIMS+mu]*eps;
    }
}

//update both momenta
void update_momenta(double eps)
{
  update_zeta_momenta(eps);
  update_lambda_momenta(eps);
}

//update the zetas
void update_zeta_positions(double eps)
{
#pragma omp parallel for
  for(int site=0;site<V;site++)
    {
      //compute the norm of pi, to form the matrix
      double pi_norm=get_zeta_norm(pi+site*N);
      
      //compute parameters of rotating matrix
      double al=eps*pi_norm;
      double cal=cos(al),sal=sin(al);
  
      //update zeta according to ortho-bound
      for(int n=0;n<N;n++)
	{
	  //get old values of coord and momenta for z
	  dcomplex x=zeta[site*N+n];
	  dcomplex p=pi[site*N+n];
	  
	  //rotate
	  zeta[site*N+n]=cal*x+sal/pi_norm*p;
	  pi[site*N+n]=-pi_norm*sal*x+cal*p;
	}
    }
}

//update the lambdas
void update_lambda_positions(double eps)
{
  //update lambda
#pragma omp parallel for
  for(int site=0;site<V;site++)
    for(int mu=0;mu<NDIMS;mu++)
      lambda[site*NDIMS+mu]*=dcomplex(cos(eps*omega[site*NDIMS+mu]),sin(eps*omega[site*NDIMS+mu]));
}

//update zetas and lambdas
void update_positions(double eps)
{
  update_zeta_positions(eps);
  update_lambda_positions(eps);
}

//integrate equation of motion
void hmc_integrate(double tl)
{
  double dt=tl/nhmc_steps/2,dth=dt/2,ldt=dt*OMELYAN_LAMBDA,l2dt=2*OMELYAN_LAMBDA*dt,uml2dt=(1-2*OMELYAN_LAMBDA)*dt;

  //     Compute H(t+lambda*dt) i.e. v1=v(t)+a[r(t)]*lambda*dt (first half step)
  update_momenta(ldt);
    
  //         Main loop
  for(int istep=0;istep<nhmc_steps;istep++)
    {
      //cout<<" Omelyan step "<<istep+1<<"/"<<nhmc_steps<<endl;
        
      //decide if last step is final or not
      double last_dt=(istep==(nhmc_steps-1)) ? ldt : l2dt;
        
      //     Compute U(t+dt/2) i.e. r1=r(t)+v1*dt/2
      update_positions(dth);
      //     Compute H(t+(1-2*lambda)*dt) i.e. v2=v1+a[r1]*(1-2*lambda)*dt
      update_momenta(uml2dt);
      //     Compute U(t+dt/2) i.e. r(t+dt)=r1+v2*dt/2
      update_positions(dth);
      //     Compute H(t+dt) i.e. v(t+dt)=v2+a[r(t+dt)]*lambda*dt (at last step) or *2*lambda*dt
      update_momenta(last_dt);
    }

  //check_lambda_conf_unitarity(lambda);
  //check_zeta_conf_unitarity(zeta);
  //normalize the configuration
#pragma omp parallel for
  for(int site=0;site<V;site++)
    {
      zeta_unitarize(zeta+site*N);
      for(int mu=0;mu<NDIMS;mu++)
	lambda_unitarize(lambda[site*NDIMS+mu]);
    }
}

//perform a hybid monte carlo update
void hmc_update(bool skip_test=false)
{
  //copy configuration
  copy_zeta_conf(zeta_old,zeta);
  copy_lambda_conf(lambda_old,lambda);
  
  //generate momenta and compute action
  generate_momenta();
  double start_mom_action=momenta_action();
  double start_theo_action=action(zeta,lambda);
  double start_topo_action=(use_topo_pot?topo_action(stout_rho,nstout_lev,lambda):0);
  double start_action=start_mom_action+start_theo_action+start_topo_action;
  //cout<<"Action: mom="<<start_mom_action<<", coord="<<start_theo_action<<", topo: "<<start_topo_action<<endl;
  
  //integrate for unitary length
  hmc_integrate(1.0);
  
  //compute final action
  double final_mom_action=momenta_action();
  double final_theo_action=action(zeta,lambda);
  double final_topo_action=(use_topo_pot?topo_action(stout_rho,nstout_lev,lambda):0);
  double final_action=final_mom_action+final_theo_action+final_topo_action;
  //cout<<"Action: mom="<<final_mom_action<<", coord="<<final_theo_action<<", topo: "<<final_topo_action<<endl;

  //compute difference of action and print it
  double diff_action=final_action-start_action;
  
  //make metropolis test
  double pacc=std::min(exp(-diff_action),1.0);
  double estr=get_unif_double(1);
  bool acc=estr<pacc;
  
  //copy back old configuration
  if(!skip_test && !acc)
    {
      copy_zeta_conf(zeta,zeta_old);
      copy_lambda_conf(lambda,lambda_old);
    }
  
  const char acc_flag[2][4]={"rej","acc"};
  cout<<"diff action: "<<final_action<<"-"<<start_action<<"="<<diff_action<<", "<<acc_flag[acc]<<endl;
}
