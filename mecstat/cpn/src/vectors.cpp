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
      printf("\"%s\" ",tag);
      printf("of %d elements of type \"%s\" (%d bytes) ",nel,type,nel*size_per_el);
      printf("allocated in file %s line %d\n",file,line);
    }
}

//initialize the first vector element
vector_el_t::vector_el_t()
{
  sprintf(tag,"base");
  sprintf(type,"(null)");
  prev=next=NULL;
  nel=0;
  size_per_el=0;
}

//initialize from a certain line
void vector_el_t::mark_as(const char *_tag,int _nel,int _size_per_el,const char *_type,const char *_file,int _line)
{
  line=_line;
  nel=_nel;
  size_per_el=_size_per_el;
  flag=0;
  take_last_characters(file,_file,VECTOR_STRING_LENGTH);
  take_last_characters(tag,_tag,VECTOR_STRING_LENGTH);
  take_last_characters(type,_type,VECTOR_STRING_LENGTH);

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
void* vectors_t::allocate(const char *tag,int nel,int size_per_el,const char *type,const char *file,int line)
{
  GET_THREAD_ID();
  if(IS_MASTER_THREAD)
    {
      int size=nel*size_per_el;
      //try to allocate the new vector
      vector_el_t *nv=(vector_el_t*)malloc(size+sizeof(vector_el_t));
      if(nv==NULL)
	CRASH("could not allocate vector named \"%s\" of %d elements of type %s (total size: %d bytes) "
	      "request on line %d of file %s",tag,nel,type,size,line,file);
      
      //append the new vector to the list
      nv->mark_as(tag,nel,size_per_el,type,file,line);
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
      if(offset!=0) CRASH("memory alignment problem, vector %s has %d offset",tag,offset);
        
      //Update the amount of required memory
      required_memory+=size;
      max_required_memory=std::max(max_required_memory,required_memory);
      
      cache_flush();
    }
    
  //sync so we are sure that master thread allocated
  THREAD_BARRIER();
  void *res=simul->returned_malloc_ptr;
    
  //resync so all threads return the same pointer
  THREAD_BARRIER();
    
  return res;
}

//release a vector
void vectors_t::deallocate(void **arr,const char *file,int line)
{
  //sync so all thread are not using the vector
  THREAD_BARRIER();
    
  GET_THREAD_ID();
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
	  required_memory-=(vect->size_per_el*vect->nel);
	  
	  //really free
	  free(vect);
	}
      else CRASH("Error, trying to delocate a NULL vector on line: %d of file: %s\n",line,file);
      
      //put to zero the array and flush the cache
      *arr=NULL;
      cache_flush();
    }
    
  //sync so all threads see that have deallocated
  THREAD_BARRIER();
}

//compute current memory usage
int vectors_t::total_memory_usage()
{
  int tot=0;
  vector_el_t *curr=main_vector;
  do
    {  
      tot+=curr->nel*curr->size_per_el;
      curr=curr->next;
    }
  while(curr!=NULL);
  
  return tot;
}
