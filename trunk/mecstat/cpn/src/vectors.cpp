#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <stdlib.h>

#include "debug.hpp"
#include "simul.hpp"
#include "threads.hpp"
#include "utils.hpp"
#include "vectors.hpp"

//print the content of a vect
void vector_el_t::print_content()
{
  if(IS_MASTER_RANK)
    {
      printf("\"%s\" ",name);
      printf("of %d bytes ",(int)size);
      printf("allocated in line %d of file %s\n",line,file);
    }
}

//initialize the first vector element
vector_el_t::vector_el_t()
{
  sprintf(name,"base");
  sprintf(file,__FILE__);
  prev=next=NULL;
  size=0;
  line=__LINE__;
}

//initialize from a certain line
void vector_el_t::mark_as(int _size,const char *_name,int _line,const char *_file)
{
  line=_line;
  size=_size;
  flag=0;
  take_last_characters(file,_file,VECTOR_STRING_LENGTH);
  take_last_characters(name,_name,VECTOR_STRING_LENGTH);
}

//////////////////////////////////////////////////////

//print all vector
void vectors_t::print_all_contents()
{
  vector_el_t *curr=main_vector;
    do
      {  
        curr->print_content();
        curr=curr->next;
      }
    while(curr!=NULL);
}

//initialize the first vector
vectors_t::vectors_t()
{
  //reset the memory asked
  required_memory=max_required_memory=0;
  
  //point last vector to main one and create it
  last_vector=main_vector=new vector_el_t;
}

//alocate a vector
void* vectors_t::allocate(int size,bool shared,const char *name,int line,const char *file)
{
  GET_THREAD_ID();
  if(IS_MASTER_THREAD)
    {
      //try to allocate the new vector
      vector_el_t *nv=(vector_el_t*)malloc(size+sizeof(vector_el_t));
      if(nv==NULL)
	CRASH_HARDLY("could not allocate vector named \"%s\" of size %d request on line %d of file %s",
		     name,size,line,file);
      
      //append the new vector to the list
      nv->mark_as(size,name,line,file);
      nv->next=NULL;
      nv->prev=last_vector;
      last_vector->next=nv;
      last_vector=nv;
        
      if(VERBOSITY_LV3)
	{
	  MASTER_PRINTF("Allocated vector ");
	  last_vector->print_content();
	}

      //define returned pointer and check for its alignement
      simul->returned_malloc_ptr=last_vector->get_pointer_to_data();
      int offset=((long long int)(simul->returned_malloc_ptr))%VECTOR_ALIGNMENT;
      if(offset!=0) CRASH_HARDLY("memory alignment problem, vector %s has %d offset",name,offset);
        
      //Update the amount of required memory
      required_memory+=size;
      max_required_memory=std::max(max_required_memory,required_memory);
      
      cache_flush();
    }
    
  //sync so we are sure that master thread allocated
  if(shared) THREAD_BARRIER();
  void *res=simul->returned_malloc_ptr;
    
  //resync so all threads return the same pointer
  if(shared) THREAD_BARRIER();
    
  return res;
}

//release a vector
void vectors_t::deallocate(void **arr,bool shared,int line,const char *file)
{
  GET_THREAD_ID();

  //sync so all thread are not using the vector
  if(shared) THREAD_BARRIER();

  if(IS_MASTER_THREAD)
    {
      if(arr!=NULL)
	{
	  vector_el_t *vect=(vector_el_t*)(*arr)-1;
	  vector_el_t *prev=vect->prev;
	  vector_el_t *next=vect->next;
            
	  if(VERBOSITY_LV3)
	    {
	      MASTER_PRINTF("At line %d of file %s freeing vector ",line,file);
	      vect->print_content();
	    }
            
	  //detach from previous
	  prev->next=next;
            
	  //if not last element
	  if(next!=NULL) next->prev=prev;
	  else last_vector=prev;
            
	  //update the required memory
	  required_memory-=(vect->size);
	  
	  //really free
	  free(vect);
	}
      
      //put to zero the array and flush the cache
      *arr=NULL;
      cache_flush();
    }
    
  //sync so all threads see that have deallocated
  if(shared) THREAD_BARRIER();
}

//compute current memory usage
int vectors_t::total_memory_usage()
{
  int tot=0;
  vector_el_t *curr=main_vector;
  do
    {  
      tot+=curr->size;
      curr=curr->next;
    }
  while(curr!=NULL);
  
  return tot;
}

//create object shared among threads
void *operator new(size_t size,bool shared,const char *name,int line,const char *file)
{
  void *ptr=simul->vectors->allocate(size,shared,name,line,file);
  if(ptr==NULL) throw new_error(size,shared,name,line,file);
  return ptr;
}
void *operator new[](size_t size,bool shared,const char *name,int line,const char *file)
{
  void *ptr=simul->vectors->allocate(size,shared,name,line,file);
  if(ptr==NULL) throw new_error(size,shared,name,line,file);
  return ptr;
}

//here to solve a problem
void deallocate(void **ptr,bool shared,int line,const char *file)
{if(*ptr!=NULL)simul->vectors->deallocate(ptr,shared,line,file);}

