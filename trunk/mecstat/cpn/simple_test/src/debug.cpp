#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

//crash
void internal_crash(int line,const char *file,const char *templ,...)
{
  //expand error message
  char mess[1024];
  va_list ap;
  va_start(ap,templ);
  vsprintf(mess,templ,ap);
  va_end(ap);
  
  fprintf(stderr,"ERROR on line %d of file \"%s\", message error: \"%s\".\n",line,file,mess);
  exit(1);
}
