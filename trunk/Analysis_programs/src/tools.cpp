#pragma once

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

string combine(const char *format,...)
{
  char buffer[1024];
  va_list args;
  
  va_start(args,format);
  vsprintf(buffer,format,args);
  va_end(args);
  
  return string(buffer);
}

double det3(double *a,double *b,double *c)
{
  double d;
  d= a[0]*(b[1]*c[2]-b[2]*c[1]);
  d+=b[0]*(c[1]*a[2]-c[2]*a[1]);
  d+=c[0]*(a[1]*b[2]-a[2]*b[1]);
  
  return d;
}

void crash(const char *temp,...)
{
  char buffer[1024];
  va_list args;

  va_start(args,temp);
  vsprintf(buffer,temp,args);
  va_end(args);

  cerr<<"ERROR: "<<buffer<<endl;
  exit(1);
}

int file_exists(const char *path)
{
  int status=0;
  
  FILE *f=fopen(path,"r");
  if(f!=NULL)
    {
      status=1;
      fclose(f);
    }
  
  return status;
}