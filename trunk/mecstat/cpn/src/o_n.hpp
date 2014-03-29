#ifndef _O_N_HPP
#define _O_N_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <cmath>
#include <tr1/array>

#include "dcompl.hpp"
#include "random.hpp"

template <int N> class O_n_t: public std::tr1::array<dcompl,N>
{
public:
  //get the norm
  double get_norm2();
  double get_norm() {return sqrt(get_norm2());}
  
  //normalize, set to one or to a random unitary vector
  void normalize(double new_norm=1);
  void set_to_one() {(*this)[0]=1;for(size_t i=1;i<N;i++) (*this)[i]=0;};
  void set_to_rnd(rnd_gen_t &gen);
  
  //summassign
  void operator+=(const O_n_t<N> &in)
  {for(size_t i=0;i<N;i++) (*this)[i]+=in[i];}
};

//extract randomly
template <int N> void O_n_t<N>::set_to_rnd(rnd_gen_t &gen)
{
  //first of all extract their norm in such: ordering
  double w[N+1];
  w[0]=0;
  for(size_t i=1;i<N;i++) w[i]=gen.get_unif(0,1);
  w[N]=1;
  std::sort(w,w+N+1);
  
  //extract a random complex number for each site
  //using the extracted norm
  for(size_t i=0;i<N;i++)
    {
      double nor=sqrt(w[i+1]-w[i]);
      //MASTER_PRINTF("nor[%d]: %lg\n",(int)i,nor);
      double the=gen.get_unif(0,2*M_PI);
      (*this)[i]=nor*dcompl(cos(the),sin(the));
    }
}

//impose the normalization condition
template <int N> void O_n_t<N>::normalize(double new_norm)
{
  double old_norm=get_norm();
  if(old_norm!=0)
    {
      double norm=new_norm/old_norm;
      for(size_t i=0;i<N;i++) (*this)[i]*=norm;
    }
}

//get the square of the norm
template <int N> double O_n_t<N>::get_norm2()
{
  double norm2=0;
  for(size_t i=0;i<N;i++) norm2+=(*this)[i].real()*(*this)[i].real()+(*this)[i].imag()*(*this)[i].imag();
  return norm2;
}

#endif
