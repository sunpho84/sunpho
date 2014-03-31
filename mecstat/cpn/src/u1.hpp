#ifndef _U1_HPP
#define _U1_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "dcompl.hpp"
#include "random.hpp"

class U1_t: public dcompl
{
public:
  void set_to_phase(double phase){real(cos(phase));imag(sin(phase));}      //set to a passed phase
  void set_to_rnd(rnd_gen_t &gen){set_to_phase(gen.get_unif(0,2*M_PI));}   //put to random
  void set_to_one(){(*this)=1;}                                            //set to 1
  
  double get_norm2(){return real()*real()+imag()*imag();}                  //return the squared norm
  void normalize(){(*this)*=1.0/sqrt(get_norm2());}                        //return it to be U1
  
  U1_t(double phase) {set_to_phase(phase);}
  U1_t() {}
};

#endif
