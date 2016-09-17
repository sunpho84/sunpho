#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <cstdlib>
#include <cstdio>
#include <cstdarg>

#include "utils.hpp"

//crash reporting the expanded error message
void crash(int line,const char *file,const char *templ,...)
{
  fflush(stdout);
  fflush(stderr);
  
  fprintf(stderr,"ERROR on line %d of file \"%s\", message error: \"",line,file);
  va_list ap;
  va_start(ap,templ);
  vfprintf(stderr,templ,ap);
  va_end(ap);
  fprintf(stderr,"\"\n");
  
  exit(1);
}
