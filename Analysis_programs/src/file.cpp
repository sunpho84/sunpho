#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void expect_string_from_file(FILE *fin,const char *exp)
{
  char read[1024];
  int nscan=fscanf(fin,"%s",read);
  if(nscan!=1||strcmp(exp,read)!=0)
    {
      if(nscan==EOF) fprintf(stderr,"Error,reached file end while waiting for '%s'\n",exp);
      else           fprintf(stderr,"Error, read '%s' instead than '%s'\n",read,exp);
      exit(1);
    }
}

void read_formatted_from_file(char *out,FILE *fin,const char *what,const char *varname)
{
  int nscan=fscanf(fin,what,out);
  if(nscan!=1)
    {
      if(nscan==EOF) fprintf(stderr,"Error,reached file end while reading '%s'\n",varname);
      else           fprintf(stderr,"Error, not enough data while reading '%s'\n",varname);
      exit(1);
    }
}

void read_formatted_from_file_expecting(char *out,FILE *fin,const char *what,const char *varname)
{
  expect_string_from_file(fin,varname);
  read_formatted_from_file(out,fin,what,varname);
}
int skip_line(FILE *fin)
{return fscanf(fin,"%*[^\n]\n");}

FILE *open_file(string path,const char* mod)
{
  FILE *out=fopen(path.c_str(),mod);

  if(out==NULL)
    {
      perror(combine("Error opening file '%s' in mode '%s'",path.c_str(),mod).c_str());
      exit(1);
    }
  
  return out;
}

int get_file_size(char *path)
{
  FILE *fin=open_file(path,"r");
  
  fseek(fin,0L,SEEK_END);
  int size=ftell(fin);
  
  fclose(fin);
  
  return size;
}
