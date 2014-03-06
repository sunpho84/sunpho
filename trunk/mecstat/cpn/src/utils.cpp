#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "utils.hpp"

//crash reporting the expanded error message
void crash(const char *templ,...)
{
  //expand error message
  char mess[1024];
  va_list ap;
  va_start(ap,templ);
  vsprintf(mess,templ,ap);
  va_end(ap);
  
  fprintf(stderr,"ERROR: \"%s\".\n",mess);
  exit(1);
}

//open a file
FILE *open_file(const char* path,const char* mod)
{
  FILE *out=fopen(path,mod);
  if(out==NULL) crash("opening file '%s' in mode '%s'",path,mod);
  
  return out;
}

//////////////////////////////////////// geometry /////////////////////////////////////

//get coords of site
void coords_of_site(coords_t coords,int site)
{
  for(int mu=0;mu<nmu;mu++)
    {
      coords[mu]=site%L;
      site/=L;
    }
}

//return site of coords
int site_of_coords(coords_t x)
{
  int site=0;
  for(int mu=0;mu<nmu;mu++) site=site*L+x[mu];
  return site;
}
