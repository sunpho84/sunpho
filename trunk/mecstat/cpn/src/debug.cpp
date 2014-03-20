#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <mpi.h>
#include <signal.h>
#include <stdio.h>
#include <unistd.h>

#include "debug.hpp"
#include "simul.hpp"

//print a message with line
void internal_shout(int line,const char *file,const char *templ,...)
{
  //give time to master thread to crash, if possible
  GET_THREAD_ID();

  if(IS_MASTER_THREAD && IS_MASTER_RANK)
    {
      //expand message
      char mess[1024];
      va_list ap;
      va_start(ap,templ);
      vsprintf(mess,templ,ap);
      va_end(ap);

      printf("line %d of file \"%s\" said: \"%s\".\n",line,file,mess);
    }
}

//called when signal received
void signal_handler(int sig)
{
  char name[100];
  switch(sig)
    {
    case SIGSEGV: sprintf(name,"segmentation violation");break;
    case SIGFPE: sprintf(name,"floating-point exception");break;
    case SIGXCPU: sprintf(name,"cpu time limit exceeded");simul->verbosity_lv=3;break;
    default: sprintf(name,"unassociated");break;
    }
  simul->vectors->print_all_contents();
  simul->print_backtrace_list();
  if(sig!=SIGXCPU) CRASH_HARDLY("signal %d (%s) detected, exiting",sig,name);
}

//crash reporting the expanded error message
void internal_crash(bool strongly,int line,const char *file,const char *templ,...)
{
  //flush everything
  fflush(stdout);
  fflush(stderr);
  
  //give time to master thread to crash, if possible
  GET_THREAD_ID();
  if(!IS_MASTER_THREAD) sleep(1);
    
  if(simul->rank_id==0)
    {
      //expand error message
      char mess[1024];
      va_list ap;
      va_start(ap,templ);
      vsprintf(mess,templ,ap);
      va_end(ap);
        
      fprintf(stderr,"ERROR on line %d of file \"%s\", message error: \"%s\".\n",line,file,mess);
      if(strongly) simul->print_backtrace_list();
      simul->abort(0);
    }
}
