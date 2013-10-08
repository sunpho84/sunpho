#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

//crash reporting the expanded error message
void crash(int line,const char *file,const char *templ,...)
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
