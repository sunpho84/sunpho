#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "data.hpp"
#include "geometry.hpp"
#include "types.hpp"

//stout the lambda configuration
void stout_lambda(dcomplex *ext_dest,double rho,dcomplex *source)
{
  //if needed allocate temp
  dcomplex *dest;
  if(ext_dest==source) dest=new dcomplex[V*NDIMS];
  else dest=ext_dest;
  
  for(int mu=0;mu<2;mu++)
    {
      int nu=!mu;                         // EDC
      for(int site=0;site<V;site++)       // FAB
	{
	  int A=site,B=neighup(A,nu);//,C=neighup(B,mu);
	  int D=neighup(A,mu),F=neighdw(A,nu),E=neighup(F,mu);
	  
	  //compute the generator
	  double Q=rho*(conj(source[A*NDIMS+mu])*(source[A*NDIMS+nu]*source[B*NDIMS+mu]*conj(source[D*NDIMS+nu])+
						  conj(source[F*NDIMS+nu])*source[F*NDIMS+mu]*source[E*NDIMS+nu]
						  )).imag();
	  
	  //set it
	  dest[site*NDIMS+mu]=dcomplex(cos(Q),sin(Q))*source[site*NDIMS+mu];
	}
    }
  
  //if needed copy back to true dest
  if(ext_dest==source)
    {
      copy_lambda_conf(ext_dest,dest);
      delete[] dest;
    }
}
