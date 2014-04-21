#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

//control that each element is unitary
void check_system_unitarity(double res)
{
  //check all sites
  for(int site=0;site<V;site++)
    {
      double dev_zeta=check_zeta_unitarity(zeta(site));
      if(dev_zeta>res) CRASH("zeta norm for site %d deviates from 1 by %lg",site,dev_zeta);

      //check all links
      for(int mu=0;mu<NDIMS;mu++)
        {
          double dev_lambda=check_lambda_unitarity(lambda(site)[mu]);
          if(dev_lambda>res) CRASH("lambda norm for site %d mu %d deviates from 1 by %lg",site,mu,dev_lambda);
        }
    }
}
