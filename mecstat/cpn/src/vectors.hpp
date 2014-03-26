#ifndef _VECTORS_HPP
#define _VECTORS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "threads.hpp"

#define VECTOR_STRING_LENGTH 20
#define VECTOR_ALIGNMENT 16

#define BLOCKING true
#define NON_BLOCKING false

#define NEW(NAME,BL_FLAG) new(BL_FLAG,NAME,__LINE__,__FILE__)
#define NEW_BLOCKING(NAME) NEW(NAME,BLOCKING)
#define NEW_NON_BLOCKING(NAME) NEW(NAME,NON_BLOCKING)

#define DELETE(A,BL_FLAG) destroy(A,BL_FLAG,__LINE__,__FILE__)
#define DELETE_BLOCKING(A) DELETE(A,BLOCKING)
#define DELETE_NON_BLOCKING(A) DELETE(A,NON_BLOCKING)

//element vector
struct vector_el_t
{
  int size,line;
  char name[VECTOR_STRING_LENGTH],file[VECTOR_STRING_LENGTH];
  vector_el_t *prev,*next;
  uint32_t flag;
  
  //padding to keep memory alignment
  char pad[VECTOR_ALIGNMENT-(2*sizeof(int)+2*sizeof(void*)+2*VECTOR_STRING_LENGTH+sizeof(uint32_t))%VECTOR_ALIGNMENT];
  
  void *get_pointer_to_data(){return (void*)(this+1);}
  void print_content();
  void mark_as(int _size,const char *_name,int _line,const char *_file);
  vector_el_t();
  ~vector_el_t();
};

//structure hosting vectors
struct vectors_t
{
  int required_memory,max_required_memory;
  vector_el_t *main_vector,*last_vector;
  
  void *allocate(int size,bool blocking,const char *name,int line,const char *file);
  void deallocate(void **arr,bool blocking,int line,const char *file);
  void print_all_contents(); //print all allocated vectors
  vectors_t();
  int total_memory_usage();
  ~vectors_t();
};

//class to catch error
class new_error
{
public:
  new_error(size_t size,bool blocking,const char *name,int line,const char *file)
  {
    const char flag[2][20]={"non_blocking","blocking"};
    internal_crash(false,line,file,"something failed creating %s object \"%s\" of size %d",flag[blocking],name,(int)size);}
};

//create object shared among threads
void *operator new(size_t size,bool blocking,const char *name,int line,const char *file);
void *operator new[](size_t size,bool blocking,const char *name,int line,const char *file);

//delete object
void deallocate(void **ptr,bool blocking,int line,const char *file);
template<class T> void destroy(T *ptr,bool blocking,int line,const char *file)
{
  GET_THREAD_ID();
  
  if(blocking==BLOCKING) THREAD_BARRIER();
  if(IS_MASTER_THREAD)
    if(ptr!=NULL)
      {
	ptr->~T();                                //call the destructor
	deallocate((void**)&ptr,false,line,file); //deallocate
      }
  if(blocking==BLOCKING) THREAD_BARRIER();
}

#endif
