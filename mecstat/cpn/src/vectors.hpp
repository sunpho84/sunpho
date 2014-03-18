#ifndef _VECTORS_HPP
#define _VECTORS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define VECTOR_STRING_LENGTH 20
#define VECTOR_ALIGNMENT 16

#define ALLOCATE(A,B,C) (C*)simul->vectors->allocate(#A,B,sizeof(C),#C,__FILE__,__LINE__)
#define FREE(A) simul->vectors->deallocate((void**)&(A),__FILE__,__LINE__)

//element vector
struct vector_el_t
{
  int nel,size_per_el,line;
  char tag[VECTOR_STRING_LENGTH],type[VECTOR_STRING_LENGTH],file[VECTOR_STRING_LENGTH];
  vector_el_t *prev,*next;
  uint32_t flag;
  
  //padding to keep memory alignment
  char pad[VECTOR_ALIGNMENT-(3*sizeof(int)+2*sizeof(void*)+3*VECTOR_STRING_LENGTH+sizeof(uint32_t))%VECTOR_ALIGNMENT];
  
  void *get_pointer_to_data(){return (void*)(this+1);}
  void print_content();
  void mark_as(const char *_tag,int _nel,int _size_per_el,const char *_type,const char *_file,int _line);
  vector_el_t();
  ~vector_el_t();
};

//structure hosting vectors
struct vectors_t
{
  int required_memory,max_required_memory;
  vector_el_t *main_vector,*last_vector;
  
  void *allocate(const char *tag,int nel,int size_per_el,const char *type,const char *file,int line);
  void deallocate(void **arr,const char *file,int line);
  void print_all_contents(); //print all allocated vectors
  vectors_t();
  int total_memory_usage();
  ~vectors_t();
};

#endif
