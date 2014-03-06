#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <mpi.h>
#include <stdio.h>
#include <unistd.h>

#include "system.hpp"

//crash reporting the expanded error message
void internal_crash(int line,const char *file,const char *templ,...)
{
  //flush everything
  fflush(stdout);
  fflush(stderr);
    
  //give time to master thread to crash, if possible
  GET_THREAD_ID();
  if(!IS_MASTER_THREAD) sleep(1);
    
  if(system->rank==0)
    {
      //expand error message
      char mess[1024];
      va_list ap;
      va_start(ap,templ);
      vsprintf(mess,templ,ap);
      va_end(ap);
        
      fprintf(stderr,"ERROR on line %d of file \"%s\", message error: \"%s\".\n",line,file,mess);
      system->print_backtrace_list();
      system->abort(0);
    }
}
