#ifndef _HMC_HPP
#define _HMC_HPP

void generate_momenta();
double momenta_action();
void hmc_integrate(double tl=1);
void hmc_update(bool skip_test=false);

#endif
