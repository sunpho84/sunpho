#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "data.hpp"
#include "geometry.hpp"
#include "topology.hpp"
#include "types.hpp"

//stout the lambda configuration
void stout_lambda(dcomplex *ext_dest,double rho,dcomplex *source)
{
  //if needed allocate temp
  dcomplex *dest;
  if(ext_dest==source) dest=new dcomplex[V*NDIMS];
  else dest=ext_dest;
  
#pragma omp parallel for
  for(int site=0;site<V;site++)
    for(int mu=0;mu<2;mu++)
      {
	int nu=!mu;
	int A=site,B=neighup(A,nu);//,C=neighup(B,mu);
	int D=neighup(A,mu),F=neighdw(A,nu),E=neighup(F,mu);
	
	//compute the generator
	double Q=rho*(conj(source[A*NDIMS+mu])*(source[A*NDIMS+nu]*source[B*NDIMS+mu]*conj(source[D*NDIMS+nu])+
						conj(source[F*NDIMS+nu])*source[F*NDIMS+mu]*source[E*NDIMS+nu]
						)).imag();
	
	//set it
	dest[site*NDIMS+mu]=dcomplex(cos(Q),sin(Q))*source[site*NDIMS+mu];
      }
  
  //if needed copy back to true dest
  if(ext_dest==source)
    {
      copy_lambda_conf(ext_dest,dest);
      delete[] dest;
    }
}

//stout the lambda configuration retaining the whole stack
void stout_lambda_whole_stack(dcomplex **out,double rho,int nstout_lev,dcomplex *in)
{
  //bind level 0 to the original
  out[0]=in;
  for(int ilev=1;ilev<=nstout_lev;ilev++) stout_lambda(out[ilev],rho,out[ilev-1]);
}

//remap the force by a single level
void stout_remap_force(dcomplex *&f_out,dcomplex *&f_in,double rho,dcomplex *l_unsm,dcomplex *l_sm)
{
#pragma omp parallel for
  for(int site=0;site<V;site++)
    for(int mu=0;mu<NDIMS;mu++)
      {
        int nu=!mu;
        
        int A=site,B=neighup(A,nu);//,C=neighup(B,mu);
        int D=neighup(A,mu),F=neighdw(A,nu),E=neighup(F,mu);
        
        double Q=rho*(conj(l_unsm[A*NDIMS+mu])*(l_unsm[A*NDIMS+nu]*l_unsm[B*NDIMS+mu]*conj(l_unsm[D*NDIMS+nu])+
						conj(l_unsm[F*NDIMS+nu])*l_unsm[F*NDIMS+mu]*l_unsm[E*NDIMS+nu]
						)).imag();
        
        dcomplex FW=l_unsm[D*NDIMS+nu]*conj(l_unsm[B*NDIMS+mu]*l_unsm[A*NDIMS+nu]);
        dcomplex BW=conj(l_unsm[E*NDIMS+nu]*l_unsm[F*NDIMS+mu])*l_unsm[F*NDIMS+nu];
        
        double LA_mu=(l_sm[A*NDIMS+mu]*f_in[A*NDIMS+mu]).real();
        double LA_nu=(l_sm[A*NDIMS+nu]*f_in[A*NDIMS+nu]).real();
        double LB_mu=(l_sm[B*NDIMS+mu]*f_in[B*NDIMS+mu]).real();
        double LD_nu=(l_sm[D*NDIMS+nu]*f_in[D*NDIMS+nu]).real();
        double LF_nu=(l_sm[F*NDIMS+nu]*f_in[F*NDIMS+nu]).real();
        double LF_mu=(l_sm[F*NDIMS+mu]*f_in[F*NDIMS+mu]).real();
        double LE_nu=(l_sm[E*NDIMS+nu]*f_in[E*NDIMS+nu]).real();
        
        dcomplex eiQ(cos(Q),sin(Q));
        
        f_out[A*NDIMS+mu]=eiQ*f_in[A*NDIMS+mu]-rho*(FW*(LA_mu-LA_nu-LB_mu+LD_nu)+BW*(LA_mu+LF_nu-LF_mu-LE_nu));
      }
}

//remap a whole stack of smearing
void stout_remap_force(dcomplex *&f,double rho,int nlevls,dcomplex **l)
{
  for(int ilev=nlevls-1;ilev>=0;ilev--)
    {
      stout_remap_force(topo_staples_supp_data,f,rho,lambda_stout[ilev],lambda_stout[ilev+1]);
      swap(topo_staples_supp_data,f);
    }
}
