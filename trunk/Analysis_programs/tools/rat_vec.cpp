#include <stdio.h>
#include <stdlib.h>

using namespace std;

FILE *open_file(char *path,const char *mode)
{
  FILE *fil=fopen(path,mode);

  if(fil==NULL)
    {
      fprintf(stderr,"Error opening file %s\n",path);
      exit(1);
    }
  
  return fil;
}

int get_nel(char *path)
{
  FILE *fin=open_file(path,"r");
  
  fseek(fin,0L,SEEK_END);
  int size=ftell(fin);
  
  fclose(fin);
  
  return size/sizeof(double);
}

void load_data(double *buf,char *path,int nel)
{
  FILE *fin=open_file(path,"r");
  int nr=fread(buf,sizeof(double),nel,fin);
  if(nr!=nel)
    {
      fprintf(stderr,"Error reading from file %s\n",path);
      exit(1);
    }
  
  fclose(fin);
}

void write_data(char *path,int nel,double *buf)
{
  FILE *fout=open_file(path,"w");
  int nw=fwrite(buf,sizeof(double),nel,fout);
  if(nw!=nel)
    {
      fprintf(stderr,"Error writing to file %s\n",path);
      exit(1);
    }
  
  fclose(fout);
}

int main(int narg,char **arg)
{
  if(narg<4)
    {
      fprintf(stderr,"Error, use %s file_num file_den file_frac\n",arg[0]);
      exit(1);
    }
  
  int nel1=get_nel(arg[1]);
  int nel2=get_nel(arg[2]);
  
  if(nel1!=nel2)
    {
      fprintf(stderr,"Error, file num and den have different number of elements!\n");
      exit(1);
    }
  
  double *num=(double*)malloc(sizeof(double)*nel1);
  double *den=(double*)malloc(sizeof(double)*nel2);
  double *fra=(double*)malloc(sizeof(double)*nel2);
  
  load_data(num,arg[1],nel1);
  load_data(den,arg[2],nel2);

  for(int iel=0;iel<nel1;iel++)
    fra[iel]=-num[iel]/den[iel];
  
  write_data(arg[3],nel2,fra);
  
  free(fra);
  free(den);
  free(num);
  
  return 0;
}
