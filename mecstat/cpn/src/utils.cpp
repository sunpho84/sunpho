#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <stdio.h>
#include <string.h>

#include "simul.hpp"

//return the log2 of N
int log2N(int N)
{
  int log2N=0;
    
  do log2N++;
  while ((2<<log2N)<N);
    
  return log2N;
}

//factorize a number
int factorize(int *list,int N)
{
  int nfatt=0;
  int fatt=2;
    
  while(N>1)
    {
      int div=N/fatt;
      int res=N-div*fatt;
      if(res!=0) fatt++;
      else 
	{
	  N=div;
	  list[nfatt]=fatt;
	  nfatt++;
	}
    }
    
  return nfatt;
}

//take the last characters of the passed string
void take_last_characters(char *out,const char *in,int size)
{
  int len=strlen(in)+1;
  int copy_len=(len<=size)?len:size;
  const char *str_init=(len<=size)?in:in+len-size;
  memcpy(out,str_init,copy_len);
}

//print a number in some kind of familiar units
void fprintf_friendly_units(FILE *fout,int quant,int orders,const char *units)
{
  const char units_prefix[6][2]={"","K","M","G","T","P"};
  int iord=0;
  double temp=quant;
    
  while(temp>orders)
    {
      temp/=orders;
      iord++;
    }
    
  quant=(int)(temp+0.5);
    
  simul->master_fprintf(fout,"%d %s%s",quant,units_prefix[iord],units);
}
void fprintf_friendly_filesize(FILE *fout,int quant)
{fprintf_friendly_units(fout,quant,1024,"Bytes");}
