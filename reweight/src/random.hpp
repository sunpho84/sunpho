#ifndef _RANDOM_HPP
#define _RANDOM_HPP

#define RAN2_NTAB 32

#include <tr1/array>

//structure for the random generator
class rnd_gen_t
{
public:
  bool inited;

  void init(int seed);
  double get_unif(double min,double max);
  int get_pm_one();
  
  rnd_gen_t(): inited(0) {} //constructor
  rnd_gen_t(int seed) {init(seed);}
private:
  int idum;
  int idum2;
  std::tr1::array<int,RAN2_NTAB> iv;
  int iy;
};

#endif
