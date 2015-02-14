#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>
#include <math.h>

#include "communicate/communicate.hpp"
#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "new_types/new_types_definitions.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif
#ifdef BGQ
 #include "bgq/intrinsic.hpp"
#endif

namespace bissa
{
  //set to zero
  void double_vector_init_to_zero(double *a,int n)
  {
    GET_THREAD_ID();
    BISSA_PARALLEL_LOOP(i,0,n) a[i]=0;
    
    set_borders_invalid(a);
  }
  
  //copy
  template <class T1,class T2> void internal_vector_copy(T1 *a,T2 *b,int n)
  {
    GET_THREAD_ID();
    BISSA_PARALLEL_LOOP(i,0,n) a[i]=b[i];
    
    set_borders_invalid(a);
  }
  
  //double to double
  THREADABLE_FUNCTION_3ARG(double_vector_copy, double*,a, double*,b, int,n)
  {internal_vector_copy(a,b,n);}THREADABLE_FUNCTION_END

  //summ
  THREADABLE_FUNCTION_3ARG(double_vector_summassign, double*,out, double*,in, int,n)
  {GET_THREAD_ID();BISSA_PARALLEL_LOOP(i,0,n) out[i]+=in[i];set_borders_invalid(out);}THREADABLE_FUNCTION_END

  //subt
  THREADABLE_FUNCTION_4ARG(double_vector_subt, double*,out, double*,in1, double*,in2, int,n)
  {GET_THREAD_ID();BISSA_PARALLEL_LOOP(i,0,n) out[i]=in1[i]-in2[i];set_borders_invalid(out);}THREADABLE_FUNCTION_END
  
  //prod with double
  THREADABLE_FUNCTION_4ARG(double_vector_prod_double, double*,out, double*,in, double,r, int,n)
  {GET_THREAD_ID();BISSA_PARALLEL_LOOP(i,0,n) out[i]=r*in[i];set_borders_invalid(out);}THREADABLE_FUNCTION_END

  //prod with double of the summ
  THREADABLE_FUNCTION_5ARG(double_vector_prod_the_summ_double, double*,out, double,r, double*,in1, double*,in2, int,n)
  {GET_THREAD_ID();BISSA_PARALLEL_LOOP(i,0,n) out[i]=r*(in1[i]+in2[i]);set_borders_invalid(out);}THREADABLE_FUNCTION_END
  
  //scalar product
  THREADABLE_FUNCTION_4ARG(double_vector_glb_scalar_prod, double*,glb_res, double*,a, double*,b, int,n)
  {
    //perform thread summ
    double loc_thread_res=0;
    GET_THREAD_ID();
    
    BISSA_PARALLEL_LOOP(i,0,n)
      loc_thread_res+=a[i]*b[i];
  }
  THREADABLE_FUNCTION_END

  //summ all points
  THREADABLE_FUNCTION_3ARG(double_vector_glb_collapse, double*,glb_res, double*,a, int,n)
  {
#ifndef REPRODUCIBLE_RUN
    //perform thread summ
    double loc_thread_res=0;
    GET_THREAD_ID();
    BISSA_PARALLEL_LOOP(i,0,n)
      loc_thread_res+=a[i];
    
    (*glb_res)=glb_reduce_double(loc_thread_res);
#else
    //perform thread summ
    float_128 loc_thread_res={0,0};
    GET_THREAD_ID();
    BISSA_PARALLEL_LOOP(i,0,n)
      float_128_summassign_64(loc_thread_res,a[i]);
    
    float_128 temp;
    glb_reduce_float_128(temp,loc_thread_res);
    (*glb_res)=temp[0];
#endif
  }
  THREADABLE_FUNCTION_END

  //put the passed vector to the new norm, returning the reciprocal of normalizating factor
  THREADABLE_FUNCTION_5ARG(double_vector_normalize, double*,ratio, double*,out, double*,in, double,norm, int,n)
  {
    //compute current norm
    double old_norm;
    double_vector_glb_scalar_prod(&old_norm,in,in,n);
    
    //compute normalizing factor
    double fact=sqrt(norm/old_norm);
    double_vector_prod_double(out,in,fact,n);
    
    if(ratio!=NULL) (*ratio)=1/fact;
    
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END

  //a[]=b[]+c[]*d
  THREADABLE_FUNCTION_6ARG(double_vector_summ_double_vector_prod_double, double*,a, double*,b, double*,c, double,d, int,n, int,OPT)
  {
    GET_THREAD_ID();

    BISSA_PARALLEL_LOOP(i,0,n) a[i]=b[i]+c[i]*d;
    
    if(!(OPT&DO_NOT_SET_FLAGS)) set_borders_invalid(a);
  }
  THREADABLE_FUNCTION_END
  
  //a[]=b[]*c+d[]*e
  THREADABLE_FUNCTION_7ARG(double_vector_linear_comb, double*,a, double*,b, double,c, double*,d, double,e, int,n, int,OPT)
  {
    GET_THREAD_ID();
    
    BISSA_PARALLEL_LOOP(i,0,n) a[i]=b[i]*c+d[i]*e;
    
    if(!(OPT&DO_NOT_SET_FLAGS)) set_borders_invalid(a);
  }
  THREADABLE_FUNCTION_END

  THREADABLE_FUNCTION_3ARG(parallel_memcpy,void*,out, void*,in, int,n)
  {
#ifdef USE_THREADS
    GET_THREAD_ID();
    
    BISSA_CHUNK_WORKLOAD(start,chunk_load,end,0,n,thread_id,NACTIVE_THREADS);
    memcpy((char*)out+start,(char*)in+start,chunk_load);
    end++;//to avoid warning
    THREAD_BARRIER();
#else
    memcpy(out,in,n);
#endif
  }
  THREADABLE_FUNCTION_END
}
