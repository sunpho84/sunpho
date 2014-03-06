#ifndef _VECTORS_HPP
#define _VECTORS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define NISSA_VECT_STRING_LENGTH 30
#define NISSA_VECT_ALIGNMENT 64

//element vector
struct vector_el_t
{
  int nel,size_per_el,line
  char tag[VECTOR_STRING_LENGTH],type[VECTOR_STRING_LENGTH],file[NISSA_VECTOR_STRING_LENGTH];
  nissa_vect *prev,*next;
  uint32_t flag;
  
  //padding to keep memory alignment
  char pad[VECTOR_ALIGNMENT-(3*sizeof(int)+2*sizeof(void*)+3*VECTOR_STRING_LENGTH+sizeof(uint32_t))%VECTOR_ALIGNMENT];
  
  vector_el_t();
  vector_el_t(const char *tag,int nel,int size_per_el,const char *type,const char *file,int line,vector_el_t *prev);
  ~vector_el_t();
};

//structure hosting vectors
struct vectors_t
{
  vectors_t();
  ~vectors_t();
  int total_memory_usage();
};

#endif
