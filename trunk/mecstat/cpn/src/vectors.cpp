#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "system.hpp"
#include "vectors.hpp"

//initialize the first vector element
vector_el_t::vector_el_t()
{
  sprintf(main_vect.tag,"base");
  sprintf(main_vect.type,"(null)");
  prev=main_vect.next=NULL;
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

//initialize the first vector
vectors_t::vectors_t()
{
  //reset the memory asked
  required_memory=max_required_memory=0;
  
  //point last vector to main one and create it
  last_vect=main_vect=new vector_el_t;
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
	system->crash("could not allocate vector named \"%s\" of %d elements of type %s (total size: %d bytes) "
		      "request on line %d of file %s",tag,nel,type,size,line,file);
      
      //append the new vector to the list
      nv->mark_as(tag,nel,size_per_el,type,file,line);
      nv->next=NULL;
      nv->prev=last_vect;        
      last_vect->next=nv;
      last_vect=nv;
        
      if(VERBOSITY_LV3)
	{
	  master_printf("Allocated vector ");
	  last_vect->printf();
	}

      //define returned pointer and check for its alignement
      system->return_malloc_ptr=last_vect->get_pointer_to_data;
      int offset=((long long int)(system_return_malloc_ptr))%VECTOR_ALIGNMENT;
      if(offset!=0) CRASH("memory alignment problem, vector %s has %d offset",tag,offset);
        
      //Update the amount of required memory
      required_memory+=size;
      max_required_memory=std::max(max_required_memory,required_memory);
      
      cache_flush();
    }
    
  //sync so we are sure that master thread allocated
  THREAD_BARRIER();
  void *res=return_malloc_ptr;
    
  //resync so all threads return the same pointer
  THREAD_BARRIER();
    
  return res;
}

//compute current memory usage
int vectors_t::total_memory_usage()
{
  int tot=0;
  vector_el_t *curr=main_vect;
  do
    {  
      tot+=curr->nel*curr->size_per_el;
      curr=curr->next;
    }
  while(curr!=NULL);
  
  return tot;
}
