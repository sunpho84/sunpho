#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

//take the last characters of the passed string
void take_last_characters(char *out,const char *in,int size)
{
  int len=strlen(in)+1;
  int copy_len=(len<=size)?len:size;
  const char *str_init=(len<=size)?in:in+len-size;
  memcpy(out,str_init,copy_len);
}
