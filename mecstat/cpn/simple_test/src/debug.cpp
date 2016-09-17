#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace std;

//crash promptin error message
void internal_crash(int line,const char *file,const char *temp,...)
{
  char buffer[1024];
  va_list args;
  
  va_start(args,temp);
  vsprintf(buffer,temp,args);
  va_end(args);
  
  cerr<<"ERROR at line "<<line<<" of file "<<file<<": "<<buffer<<endl;
  exit(1);
}
