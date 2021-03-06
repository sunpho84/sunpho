#ifndef _IOS_H
#define _IOS_H

#include <string>
#include "new_types/new_types_definitions.hpp"

namespace bissa
{
  int count_substrings(const char *str,const char *sub);
  FILE* open_file(const char *outfile,const char *mode);
  FILE* open_text_file_for_output(const char *outfile);
  int cd(const char *path);
  int cp(char *path_out,char *path_in);
  int create_dir(char *path);
  int master_fprintf(FILE *stream,const char *format,...);
  //int rm(const char *path);
  std::string combine(const char *format,...);
  void close_file(FILE *file);
  void fprintf_friendly_filesize(FILE *fout,int quant);
  void fprintf_friendly_units(FILE *fout,int quant,int orders,const char *units);
  void take_last_characters(char *out,const char *in,int size);
  int count_file_lines(const char *path);
  int get_file_size(const char *path);
}

#endif
