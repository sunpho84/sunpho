#ifndef _RANDOM_H
#define _RANDOM_H

namespace bissa
{
  void convert_rnd_gen_to_text(char *text,rnd_gen *gen);
  double rnd_get_unif(rnd_gen *gen,double min,double max);
  int rnd_get_pm_one(rnd_gen *gen);
  void rnd_fill_pm_one_loc_vector(double *v,int nps);
  void rnd_fill_unif_loc_vector(double *v,int dps,double min,double max);
  void start_glb_rnd_gen(int seed);
  void start_loc_rnd_gen(int seed);
  void start_loc_rnd_gen(char *mess);
  void start_rnd_gen(rnd_gen *out,int seed);
  void stop_loc_rnd_gen();
}

#endif
