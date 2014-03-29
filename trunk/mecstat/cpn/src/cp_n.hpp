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
  field_t<U1_t > Lambda;
  
  enum hot_cold{HOT,COLD}; //enumerator for condition to start
  
  //constructor and destructor
  CP_N_t(int L,double beta,hot_cold hc);
  ~CP_N_t();

  //cool down everything or heat it
  void set_to_cold();
  void set_to_hot();
  void set_to(hot_cold hc);
  
  //return the norms
  double get_Zeta_norm();
  double get_Lambda_norm();
private:
  CP_N_t(const CP_N_t &in) {}
};

//set as hot
template <int N> void CP_N_t<N>::set_to_hot()
{
  GET_THREAD_ID();
  PARALLEL_FOR_SITES_OF_FIELD(site,Zeta)
    //if(thread_id==0)
    {
      Zeta[site].set_to_rnd(geometry->loc_rnd_gen[site]);
      Lambda[site].set_to_rnd(geometry->loc_rnd_gen[site]);
      double r=Lambda[site].real(),i=Lambda[site].imag();
      printf("site: %d, %lg %p %lg %lg\n",site,Lambda[site].get_norm2(),&(Lambda[site]),r,i);
    }
  Zeta.mark_touched();
  Lambda.mark_touched();
  PARALLEL_FOR_SITES_OF_FIELD(site,Lambda)
    //if(thread_id==0)
    printf("site: %d, %lg %lg %lg\n",site,Lambda[site].get_norm2(),Lambda[site].real(),Lambda[site].imag());
  THREAD_BARRIER();
}

//set as cold
template <int N> void CP_N_t<N>::set_to_cold()
{
  GET_THREAD_ID();
  PARALLEL_FOR_SITES_OF_FIELD(site,Zeta)
    {
      Zeta[site].set_to_one();
      Lambda[site].set_to_one();
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
