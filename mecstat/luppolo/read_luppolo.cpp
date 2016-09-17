#include <cstdio>
#include <cstdlib>

int main()
{
  FILE *fin=fopen("/tmp/Luppoli","r");
  if(fin==NULL)
    {
      fprintf(stderr,"error opening\n");
      exit(1);
    }

  double data[12][12][12][2];
  int rc=fread(data,sizeof(data),1,fin);
  if(rc!=1)
    {
      fprintf(stderr,"erorr loading\n");
      exit(1);
    }

  for(int i=0;i<12;i++)
    printf("%lg\tf%lg\t%lg\n",data[i][0][0][0],data[0][i][0][0],data[0][0][i][0]);
  
  return 0;
}
