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

template<class T> T det3(T *a,T *b,T *c)
{
  T d;
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

string smart_print(double m,double e)
{
  double orim=m;
  int s=0;
  ostringstream o;
  
  if(e==0) o<<m<<"(0)";
  else
    if(e<1)
      {
        //print 1 digit of error or 2 if err starts with 1 or 2
        while(int(e)<=3)
          {
            e*=10;
            m*=10;
            s--;
          }
        
        char t[1000];
        sprintf(t,"%.*f",(unsigned int)abs(s),orim);
        
        o<<t<<"("<<int(e+0.5)<<")";
      }
    else
      {
        if(e>=3)
          {
            //count the numbr of digits to truncate
            int s=0;
            /*
              while(int(e)>=30)
              {
              m/=10;
              e/=10;
              s++;
              }
	    */
            o<<(int)(m+0.5)*pow(10,s)<<"("<<(int)(e+0.5)*pow(10,s)<<")";
          }
        else
          {
            char t[1000];
            sprintf(t,"%.1f(%.1f)",m,e);
            o<<t;
          }
      }
  
  return o.str();
}
