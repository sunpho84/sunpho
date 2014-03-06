#ifndef _SYSTEM_HPP
#define _SYSTEM_HPP

#include "random.hpp"

typedef bool spin_t;

//simulating system
struct system_t
{
  rnd_gen_t glb_rnd_gen,*loc_rnd_gen;                 //random number generators
  int total_flipped;                                  //total number of flipped per copy
  int glb_par_link,glb_up_spins;                      //store the total number of parallel link and up spins
  spin_t *spins;                                      //configuration
  int *next_cluster,*curr_cluster;                    //cluster at the end of the step and at the beginning
  
  void print_spins();
  system_t(int seed,double filling);
  int change_single_cluster();
  double get_energy(){return -glb_par_link*2;}
  double get_magnetization(int nsites){return (glb_up_spins*2-nsites)/(double)(nsites);}
  ~system_t();
};

#endif
