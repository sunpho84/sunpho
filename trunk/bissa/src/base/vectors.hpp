#ifndef _VECTORS_H
#define _VECTORS_H

#include <stdio.h>
#include "new_types/new_types_definitions.hpp"

namespace bissa
{
  char *get_vect_name(void *v);
  int check_borders_allocated(void *data);
  int check_borders_communicated_at_least_once(void *data);
  int check_borders_valid(void *data);
  int compute_vect_memory_usage();
  int get_vect_flag(void *v,unsigned int flag);
  bissa_vect* get_vect(void *v);
  void *internal_bissa_malloc(const char *tag,int nel,int size_per_el,const char *type,const char *file,int line);
  void crash_if_borders_not_allocated(void *v);
  void ignore_borders_communications_warning(void *data);
  void initialize_main_vect();
  void internal_bissa_free(char **arr,const char *file,int line);
  void vector_copy(void *a,void *b);
  void vector_reset(void *a);
  void last_vect_content_printf();
  void vect_content_fprintf(FILE *fout,bissa_vect *vect);
  void vect_content_printf(bissa_vect *vect);
  void print_all_vect_content();
  void reorder_vector(char *vect,int *order,int nel,int sel);
  void set_borders_invalid(void *data);
  void set_borders_valid(void *data);
  void set_vect_flag(void *v,unsigned int flag);
  void set_vect_flag_non_blocking(void *v,unsigned int flag);
  void unset_vect_flag(void *v,unsigned int flag);
  void unset_vect_flag_non_blocking(void *v,unsigned int flag);
  void vect_content_fprintf(FILE *f,void *vec);
  void vect_content_printf(void *vec);
}
#endif
  
