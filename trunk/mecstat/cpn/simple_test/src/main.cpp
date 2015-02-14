#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <fstream>
#include <iostream>

#include "cpn.hpp"

using namespace std;

inline int get_unif_int(int max,int site)
{
  uniform_int_distribution<int> dis(0,max-1);
  return dis(gen[site]);
}

int acc1=0,acc2=0;
int nestr1=0,nestr2=0;
double w=0.3;
  
void metro_site(int O)
{
  //cout<<"0 "<<get_zeta_norm(zeta+O*N)<<endl;
  int n=get_unif_int(N,O);
  double ebef=site_action(zeta,lambda,O);
  dcomplex zbef=zeta[O*N+n];
  double th=(get_unif_double(2*M_PI,O)-M_PI)*w;
  zeta[O*N+n]*=dcomplex(cos(th),sin(th));
  double eaft=site_action(zeta,lambda,O);
  double diff=eaft-ebef;
  double tresh=exp(-diff);
  double extr=get_unif_double(1,O);
  if(extr<tresh) acc1++;//cout<<"Acc: "<<extr<<"<exp("<<-diff<<")="<<tresh<<endl;
  else
    {
      //cout<<"Rej: "<<extr<<">exp("<<-diff<<")="<<tresh<<endl;
      zeta[O*N+n]=zbef;
    }
  
  //cout<<"1 "<<get_zeta_norm(zeta+O*N)<<endl;
  nestr1++;
}

void metro_site2(int O)
{
  int i=get_unif_int(N,O);
  int j;
  do j=get_unif_int(N,O);
  while(i==j);
  
  double ebef=site_action(zeta,lambda,O);
  dcomplex zibef=zeta[O*N+i];
  dcomplex zjbef=zeta[O*N+j];
  double th=(get_unif_double(2*M_PI,O)-M_PI)*w;
  double cth=cos(th);
  double sth=sin(th);
  zeta[O*N+i]=+cth*zibef-sth*zjbef;
  zeta[O*N+j]=+sth*zibef+cth*zjbef;
  //cout<<"th: "<<th<<" i: "<<i<<", j: "<<j<<endl;
  //cout<<"Pre: "<<zibef<<" "<<zjbef<<" "<<norm(zibef)+norm(zjbef)<<endl;
  //cout<<"Aft: "<<zeta[O*N+i]<<" "<<zeta[O*N+j]<<" "<<norm(zeta[O*N+i])+norm(zeta[O*N+j])<<endl;
  double eaft=site_action(zeta,lambda,O);
  double diff=eaft-ebef;
  double tresh=exp(-diff);
  double extr=get_unif_double(1,O);
  if(extr<tresh) acc2++;//cout<<"Acc2: "<<extr<<"<exp("<<-diff<<")="<<tresh<<endl;
  else
    {
      //cout<<"Rej2: "<<extr<<">exp("<<-diff<<")="<<tresh<<endl;
      zeta[O*N+i]=zibef;
      zeta[O*N+j]=zjbef;
    } 

  nestr2++;
  
  //cout<<get_zeta_norm(zeta+O*N)<<endl;
}

void metro_sweep_new()
{
  for(int site=0;site<V;site++)
    {
      for(int n=0;n<N*(N-1)/2*2;n++)
	{
	  for(int j=0;j<2;j++) metro_site(site);
	  metro_site2(site);
	}
      for(int mu=0;mu<2;mu++) overheat_update_link(site,mu);
    }

  cout<<"Acc1: "<<acc1/(double)nestr1<<endl;
  cout<<"Acc2: "<<acc2/(double)nestr2<<endl;
}

int main()
{
  //read input and initialize
  read_pars_t read_pars;
  read_input(read_pars,"input");
  int base_isweep;
  init(base_isweep,read_pars);
  
  //open output files
  ios::openmode mode=ios::out;
  if(read_pars.start_cond==LOAD) mode|=ios::app;
  ofstream energy_file("energy",mode);
  ofstream polyakov_file("polyakov",mode);
  ofstream topology_file("topology",mode);
  ofstream corr_file("corr",mode);
  ofstream corrd_file("corrd",mode);
  ofstream mom2_file("mom2",mode);
  ofstream mag_file("mag",mode);
  energy_file.precision(PREC);
  polyakov_file.precision(PREC);
  topology_file.precision(PREC);
  corr_file.precision(PREC);
  corrd_file.precision(PREC);
  mag_file.precision(PREC);
  
  //sweep with overheat-micro
  timing_t tot_time,sweep_time,energy_time,geo_topo_time,topo_time,corr_time;
  
  tot_time.start();
  int isweep;
  for(isweep=base_isweep;isweep<base_isweep+read_pars.nsweep;isweep++)
    {
      //compute geometrical topological charge
      geo_topo_time.start();
      double topo_sim=geometric_topology_simplified(zeta);
      geo_topo_time.stop();

      /*
      for(int i=0;i<100;i++)
	{
	  int O=1;
	  int mu=0,nu=1;
	
	  double topo=0;
	  int A=neighup(O,mu);
	  int B=neighup(A,nu);
	  int C=neighup(O,nu);
	  int D=neighdw(O,mu);
	  int E=neighdw(D,nu);
	  int F=neighdw(O,nu);
	  dcomplex *z=zeta;
	  dcomplex AB=get_zeta_compl_scalprod(z+A*N,z+B*N);
	  dcomplex BC=get_zeta_compl_scalprod(z+B*N,z+C*N);
	  dcomplex CD=get_zeta_compl_scalprod(z+C*N,z+D*N);
	  dcomplex DE=get_zeta_compl_scalprod(z+D*N,z+E*N);
	  dcomplex EF=get_zeta_compl_scalprod(z+E*N,z+F*N);
	  dcomplex FA=get_zeta_compl_scalprod(z+F*N,z+A*N);
	  topo+=
	    arg(get_zeta_compl_scalprod(z+B*N,z+O*N)*
		AB*
		get_zeta_compl_scalprod(z+O*N,z+A*N))+
	    arg(get_zeta_compl_scalprod(z+C*N,z+O*N)*
		BC*
		get_zeta_compl_scalprod(z+O*N,z+B*N))+
	    arg(get_zeta_compl_scalprod(z+O*N,z+E*N)*
		get_zeta_compl_scalprod(z+F*N,z+O*N)*
		EF)+
	    arg(DE*
		get_zeta_compl_scalprod(z+O*N,z+D*N)*
		get_zeta_compl_scalprod(z+E*N,z+O*N))+
	    arg(CD*
		get_zeta_compl_scalprod(z+O*N,z+C*N)*
		get_zeta_compl_scalprod(z+D*N,z+O*N))+
	    arg(get_zeta_compl_scalprod(z+O*N,z+F*N)*
		get_zeta_compl_scalprod(z+A*N,z+O*N)*
		FA);
	  topo/=2*M_PI;

	  double tgeo=0;//geometric_topology_simplified(zeta);
	  double ene=0;//energy(zeta,lambda)/V/NDIMS;
	  cout<<"Topo: "<<topo<<" "<<tgeo<<" "<<ene<<endl;
	  
	  metro_site(O);
	  metro_site2(O);
	}
      */
      //compute energy
      energy_time.start();
      energy_file<<isweep<<" "<<energy(zeta,lambda)/V/NDIMS<<endl;
      energy_time.stop();
      
      //compute polyakov loop
      dcomplex poly=polyakov(lambda);
      polyakov_file<<isweep<<" "<<poly.real()<<" "<<poly.imag()<<endl;
      
      //compute topologycal charge
      topo_time.start();
      stout_lambda_whole_stack(lambda_stout,stout_rho,nstout_lev,lambda);
      for(int ilev=0;ilev<=nstout_lev;ilev++)
	{
	  double topo_num=topology(lambda_stout[ilev]);
	  if(use_topo_pot==2 && ilev==nstout_lev && isweep>=chrono_topo_after && (isweep-chrono_topo_after)%chrono_topo_each==0)
	    {
	      update_chrono_potential(+topo_num);
	      if(isweep%DRAW_EACH==0)
		{
		  draw_chrono_topo_potential();
		  draw_chrono_topo_force();
		}
	    }
	  
	  topology_file<<isweep<<" "<<ilev<<" "<<
	    topo_sim<<" "<<
	    topo_num<<endl;
	}
      topo_time.stop();
      
      //compute the correlation function
      if(isweep%compute_corr_each==0)
	{
	  corr_time.start();
	  double mag0,mag1,mom2,corr[L],corrd[L];
	  //compute_corr_alt(corr,zeta);
	  compute_corr(mag0,mag1,mom2,corr,corrd,zeta);
	  for(int i=0;i<=L/2;i++) corr_file<<isweep<<" "<<i<<" "<<corr[i]<<endl;
	  for(int i=0;i<=L/2;i++) corrd_file<<isweep<<" "<<i/sqrt(2)<<" "<<corrd[i]<<endl;
	  mom2_file<<isweep<<" "<<mom2<<endl;
	  mag_file<<isweep<<" "<<mag0<<" "<<mag1<<endl;
	  corr_time.stop();
	}
  
      //sweep
      sweep_time.start();
      switch(read_pars.use_hmc)
	{
	case 1:
	   hmc_update(isweep<read_pars.nterm);
	   break;
	case 0:
	  for(int imicro=0;imicro<read_pars.nmicro;imicro++) micro_sweep();
	  overheat_sweep();
	  break;
	case -1:
	  metro_sweep_new();
	  break;
	default:
	  crash("unkwnown update %d",read_pars.use_hmc);
	}
      sweep_time.stop();
    }
  
  tot_time.stop();
  
  //write lasted time
  cout<<"Acc: "<<nacc/(double)read_pars.nsweep<<endl;
  cout<<"Tot time: "<<tot_time<<endl;
  cout<<"Sweep time: "<<sweep_time<<endl;
  cout<<"Geo topo time: "<<geo_topo_time<<endl;
  cout<<"Topo time: "<<topo_time<<endl;
  cout<<"Energy time: "<<energy_time<<endl;
  cout<<"Corr time: "<<corr_time<<endl;
  
  //write the conf
  write_conf("conf",isweep);
  if(use_topo_pot==2) draw_chrono_topo_potential();
  
  print_rand_stat();
  
  //finalize
  close();
  
  return 0;
}
