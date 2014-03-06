#ifndef _RANDOM_HPP
#define _RANDOM_HPP

#define RAN2_NTAB 32

//structure for the random generator
struct rnd_gen_t
{
  int idum;
  int idum2;
  int iv[RAN2_NTAB];
  int iy;
  void init(int seed);
  double rnd_get_unif(double min,double max);
  int rnd_get_pm_one();
};

#endif
