#include <iostream>
#include <cstring>

#include "cpn.hpp"

using namespace std;

dcomplex *pi_old;
double *omega_old;

void copy_pi_momenta(dcomplex *dest,dcomplex *source)
{
#pragma omp parallel for
  for(int site=0;site<V;site++)
    for(int n=0;n<N;n++)
      dest[site*N+n]=source[site*N+n];
}

void copy_omega_momenta(double *dest,double *source)
{
#pragma omp parallel for
  for(int site=0;site<V;site++)
    for(int mu=0;mu<NDIMS;mu++)
      dest[site*NDIMS+mu]=source[site*NDIMS+mu];
}

void check_momenta_orthogonality()
{
  double pi_scp=0;
  for(int site=0;site<V;site++)
    {
      double proj=0;
      for(int n=0;n<N;n++)
	proj+=(pi[site*N+n]*conj(zeta[site*N+n])).real();
      pi_scp=max(pi_scp,proj);
    }
  
  cout<<"Pi-Zeta scalar product: "<<pi_scp<<endl;
}

void check_fields_unitarity()
{
  double lambda_nonun=0;
  for(int site=0;site<V;site++)
    for(int mu=0;mu<NDIMS;mu++)
      lambda_nonun=max(lambda_nonun,get_lambda_norm(lambda[site*NDIMS+mu])-1);

  double zeta_nonun=0;
  for(int site=0;site<V;site++)
    zeta_nonun=max(zeta_nonun,get_zeta_norm(zeta+site*N)-1);

  cout<<"Lambda nonun: "<<lambda_nonun<<endl;
  cout<<"Zeta nonun: "<<zeta_nonun<<endl;
}

void check_momenta_variation()
{
  double omega_diff=0;
  for(int site=0;site<V;site++)
    for(int mu=0;mu<NDIMS;mu++)
      omega_diff=max(omega_diff,omega[site*NDIMS+mu]-omega_old[site*NDIMS+mu]);

  double pi_diff=0;
  for(int site=0;site<V;site++)
    for(int n=0;n<N;n++)
      pi_diff=max(pi_diff,sqrt(norm(pi[site*N+n]-pi_old[site*N+n])));

  cout<<"Omega diff: "<<omega_diff<<endl;
  cout<<"Pi diff: "<<pi_diff<<endl;
}

void check_fields_variation()
{
  double lambda_diff=0;
  for(int site=0;site<V;site++)
    for(int mu=0;mu<NDIMS;mu++)
      lambda_diff=max(lambda_diff,sqrt(norm(lambda[site*NDIMS+mu]-lambda_old[site*NDIMS+mu])));

  double zeta_diff=0;
  for(int site=0;site<V;site++)
    for(int n=0;n<N;n++)
      zeta_diff=max(zeta_diff,sqrt(norm(zeta[site*N+n]-zeta_old[site*N+n])));

  cout<<"Lambda diff: "<<lambda_diff<<endl;
  cout<<"Zeta diff: "<<zeta_diff<<endl;
}

void test_reversibility(int nt)
{
  cout<<"----------------------"<<endl;
  cout<<"     nt   =   "<<nt<<endl;
  cout<<"----------------------"<<endl;
  nhmc_steps=nt;
  
  //copy configuration
  copy_zeta_conf(zeta,zeta_old);
  copy_lambda_conf(lambda,lambda_old);
  copy_pi_momenta(pi,pi_old);
  copy_omega_momenta(omega,omega_old);
  
  //compute action
  double start_mom_action=momenta_action();
  double start_theo_action=action(zeta,lambda);
  double start_topo_action=(use_topo_pot?topo_action(stout_rho,nstout_lev,lambda):0);
  double start_action=start_mom_action+start_theo_action+start_topo_action;
  cout<<"Action: mom="<<start_mom_action<<", coord="<<start_theo_action<<", topo: "<<start_topo_action<<endl;
  
  //integrate
  hmc_integrate(+1);

  //intermediate action
  double inter_mom_action=momenta_action();
  double inter_theo_action=action(zeta,lambda);
  double inter_topo_action=(use_topo_pot?topo_action(stout_rho,nstout_lev,lambda):0);
  double inter_action=inter_mom_action+inter_theo_action+inter_topo_action;
  cout<<"Action: mom="<<inter_mom_action<<", coord="<<inter_theo_action<<", topo: "<<inter_topo_action<<endl;
  double diff_action=inter_action-start_action;
  cout<<"Diff action: "<<diff_action<<endl; 

  //integrate back
  hmc_integrate(-1);
  
  check_momenta_variation();
  check_fields_variation();
  check_momenta_orthogonality();
  check_fields_unitarity();
  
  //compute final action
  double final_mom_action=momenta_action();
  double final_theo_action=action(zeta,lambda);
  double final_topo_action=(use_topo_pot?topo_action(stout_rho,nstout_lev,lambda):0);
  double final_action=final_mom_action+final_theo_action+final_topo_action;
  cout<<"Action: mom="<<final_mom_action<<", coord="<<final_theo_action<<", topo: "<<final_topo_action<<endl;

  //compute difference of action and print it
  diff_action=final_action-start_action;
  cout<<"Diff action: "<<diff_action<<endl; 
}

int main()
{
  //read input and initialize
  read_pars_t read_pars;
  read_input(read_pars,"input");
  int base_isweep;
  init(base_isweep,read_pars);
  
  pi_old=new dcomplex[V*N];
  omega_old=new double[V*NDIMS];
  
  if(!read_pars.use_hmc) crash("implicit in the test to use HMC");
  generate_momenta();

  for(int nt=1;nt<18;nt++)
    {
      copy_zeta_conf(zeta_old,zeta);
      copy_lambda_conf(lambda_old,lambda);
      copy_pi_momenta(pi_old,pi);
      copy_omega_momenta(omega_old,omega);
      
      test_reversibility(nt);
    }
  
  delete[] pi_old;
  delete[] omega_old;
  
  return 0;
}
