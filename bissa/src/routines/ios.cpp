#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>
#include <mpi.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <sys/stat.h>

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "io/input.hpp"

#include "mpi_routines.hpp"
#ifdef USE_THREADS
 #include "thread.hpp"
#endif

namespace bissa
{
  //return the number of occurency of "sub" inside "str"
  int count_substrings(const char *str,const char *sub)
  {
    int rc=0;
    for(str=(str!=NULL)?strstr(str,sub):str;str!=NULL;str=strstr(str+1,sub)) rc++;

    return rc;
  }
  
  //only master rank and thread print
  int master_fprintf(FILE *stream,const char *format,...)
  {
    GET_THREAD_ID();
    int ret=0;
    
    if(rank==0 && IS_MASTER_THREAD)
      {
	va_list ap;
	va_start(ap,format);
	ret=vfprintf(stream,format,ap);
	va_end(ap);
      }
    
    return ret;
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
    
    master_fprintf(fout,"%d %s%s",quant,units_prefix[iord],units);
  }
  void fprintf_friendly_filesize(FILE *fout,int quant)
  {fprintf_friendly_units(fout,quant,1024,"Bytes");}
  
  //create a dir
  int create_dir(char *path)
  {
    int res=(rank==0) ? mkdir(path,480) : 0;
    MPI_Bcast(&res,1,MPI_INT,0,MPI_COMM_WORLD);
    if(res!=0)
      master_printf("Warning, failed to create dir %s, returned %d. Check that you have permissions and that parent dir exists.\n",path,res);
    
    return res;
  }
  
  //copy a file
  int cp(char *path_out,char *path_in)
  {
    int rc=0;
    if(rank==0)
      {
	char command[1024];
	sprintf(command,"cp %s %s",path_in,path_out);
	rc=system(command);
	if(rc!=0) crash("cp failed!");
      }
    
    return master_broadcast(rc);
  }
  
  //pass to the folder
  int cd(const char *path)
  {
    int rc=0;
    if(rank==0)
      {
	char command[1024];
	sprintf(command,"cd %s",path);
	rc=system(command);
	if(rc!=0) crash("cd failed!");
      }
    
    return master_broadcast(rc);
  }
  
  //Open a file checking it
  FILE* open_file(const char *outfile,const char *mode)
  {
    FILE *fout=NULL;
    
    if(rank==0)
      {
	fout=fopen(outfile,mode);
	if(fout==NULL) crash("Couldn't open the file: %s for mode: %s",outfile,mode);
      }
    
    return fout;
  }
  
  //Open a text file for output
  FILE* open_text_file_for_output(const char *outfile)
  {return open_file(outfile,"w");}
  
  //Open a text file for input
  FILE* open_text_file_for_input(const char *infile)
  {return open_file(infile,"r");}
  
  //close an open file
  void close_file(FILE *file)
  {if(rank==0) fclose(file);}
  
  //count the number of lines in a file
  int count_file_lines(const char *path)
  {
    //return -1 if file does not exist
    if(!file_exists(path)) return -1;
    
    //scan the file
    FILE *fin=open_text_file_for_input(path);
    int n=0;
    if(rank==0)
      {
	int ch;
	while(EOF!=(ch=getchar())) if (ch=='\n') n++;
      }
    
    //close file and broadcast n
    close_file(fin);
    
    return master_broadcast(n);
  }
  
  //get the size of a file
  int get_file_size(const char *path)
  {
    //return -1 if file does not exist
    if(!file_exists(path)) return -1;
    
    //scan the file
    FILE *fin=open_text_file_for_input(path);
    int file_size=0;
    if(rank==0)
      {
	if(fseek(fin,0,SEEK_END)) crash("while seeking");
	file_size=ftell(fin);
      }
    
    return master_broadcast(file_size);
  }
  
  //take the last characters of the passed string
  void take_last_characters(char *out,const char *in,int size)
  {
    int len=strlen(in)+1;
    int copy_len=(len<=size)?len:size;
    const char *str_init=(len<=size)?in:in+len-size;
    memcpy(out,str_init,copy_len);
  }
  
  //combine arguments in a single string
  std::string combine(const char *format,...)
  {
    char buffer[1024];
    va_list args;
    
    va_start(args,format);
    vsprintf(buffer,format,args);
    va_end(args);
    
    return std::string(buffer);
  }
}
