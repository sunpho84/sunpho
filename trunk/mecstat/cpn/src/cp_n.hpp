#ifndef _CPN_HPP
#define _CPN_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "field.hpp"
#include "o_n.hpp"
#include "u1.hpp"

template <int N> class CP_N_t
{
public:
  double beta;
  geometry_t *geometry;
  field_t<O_n_t<N> > Zeta;
  field_t<std::tr1::array<U1_t,2> > Lambda;
  
  enum hot_cold{HOT,COLD}; //enumerator for condition to start
  
  //constructor and destructor
  CP_N_t(int L,double beta,hot_cold hc);
  ~CP_N_t();

  //cool down everything or heat it
  void set_to_cold();
  void set_to_hot();
  void set_to(hot_cold hc);
  
  //return g and energy
  double g(){return 1/(N*beta);}
  double get_energy();
  double get_action(){return get_energy()/g();}
  double get_site_energy(int site);
  double get_site_action(int site){return get_site_energy(site)/g();}
  
private:
  CP_N_t(const CP_N_t &in) {}
};

//set as hot
template <int N> void CP_N_t<N>::set_to_hot()
{
  GET_THREAD_ID();
  PARALLEL_FOR_SITES_OF_FIELD(site,Zeta)
    {
      Zeta[site].set_to_rnd(geometry->loc_rnd_gen[site]);
      //for(int n=0;n<N;n++)
      //{
      //double r=Zeta[site][n].real(),i=Zeta[site][n].imag();
      // printf("site: %d %d, %lg %lg\n",geometry->glb_site_of_loc_site(site),n,r,i);
      //}
      for(int mu=0;mu<2;mu++)
	Lambda[site][mu].set_to_rnd(geometry->loc_rnd_gen[site]);
    }
  Zeta.mark_touched();
  Lambda.mark_touched();
  //PARALLEL_FOR_SITES_OF_FIELD(site,Lambda)
  //if(thread_id==0)
  //printf("site: %d, %lg %lg %lg\n",site,Lambda[site].get_norm2(),Lambda[site].real(),Lambda[site].imag());
  //THREAD_BARRIER();
}

//set as cold
template <int N> void CP_N_t<N>::set_to_cold()
{
  GET_THREAD_ID();
  PARALLEL_FOR_SITES_OF_FIELD(site,Zeta)
    {
      Zeta[site].set_to_one();
      for(int mu=0;mu<2;mu++) Lambda[site][mu].set_to_one();
    }
  Zeta.mark_touched();
  Lambda.mark_touched();
}

//determine the correct one
template <int N> void CP_N_t<N>::set_to(hot_cold hc)
{
  switch(hc)
    {
    case HOT: set_to_hot();break;
    case COLD: set_to_cold();break;
    default:CRASH_SOFTLY("unknown condition");
    }
}

//compute the energy
template <int N> double CP_N_t<N>::get_energy()
{
  GET_THREAD_ID();
  
  //sync borders
  Zeta.sync_outer_sites();
  Lambda.sync_outer_sites();
  
  //compute local energy
  double loc_energy=0;
  PARALLEL_FOR_SITES_OF_FIELD(site,Zeta)
    for(size_t mu=0;mu<2;mu++)
      {
	int site_up=Zeta.get_neigh(site,2*mu+1);
	for(int n=0;n<N;n++)
	  loc_energy+=(conj(Zeta[site_up][n])*Zeta[site][n]*Lambda[site][mu]).real();
      }
  return -(2*rank_threads_reduce(loc_energy)-2*geometry->nglb_sites*2);
}

//return the energy of a single site
template <int N> double CP_N_t<N>::get_site_energy(int site)
{
  double res=0;
  for(int mu=0;mu<2;mu++)
    {
      if(Zeta.get_neigh(Zeta.get_neigh(site,2*mu+0),2*mu+1)!=site) CRASH_SOFTLY("%d %d %d",mu,Zeta.get_neigh(Zeta.get_neigh(site,2*mu+0),2*mu+1),site);
      int site_up=Zeta.get_neigh(site,2*mu+1);
      for(int n=0;n<N;n++)
        res+=(conj(Zeta[site_up][n])*Zeta[site][n]*Lambda[site][mu]).real();
      int site_dw=Zeta.get_neigh(site,2*mu+0);
      for(int n=0;n<N;n++)
        res+=(conj(Zeta[site_dw][n])*Zeta[site][n]*conj(Lambda[site_dw][mu])).real();
    }

  return -(res-2*2);
}

//constructor
template <int N> CP_N_t<N>::CP_N_t(int L,double beta,hot_cold hc): beta(beta)
{
  //init the geometry and the random number generator
  geometry=NEW_BLOCKING("geometry") geometry_t(2,L);
  geometry->init_loc_rnd_gen();
  
  //init the field Zeta and Lambda
  Zeta.create("Zeta",geometry->first_neighbors);
  Lambda.create("Lambda",geometry->first_neighbors);
  
  //init to the correct condition
  set_to(hc);
}

//destructor
template <int N> CP_N_t<N>::~CP_N_t()
{
  DELETE_BLOCKING(geometry);
}

#endif
