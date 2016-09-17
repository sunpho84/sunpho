#include "include.h"

const size_t njack=16;

int main(int narg,char **arg)
{
  if(narg<3) crash("use %s in out",arg[0]);
  
  size_t nbytes=get_file_size(arg[1]);
  size_t ndoubles=nbytes/sizeof(double);
  if(ndoubles*sizeof(double)!=nbytes) crash("file size %lu is not a multiple of 8",nbytes);
  size_t nel=ndoubles/(njack+1);
  if((njack+1)*nel!=ndoubles) crash("number of double %lu is not a multiple of number of jacknives+1 %lu",
				    ndoubles,njack+1);
  
  double *data=(double*)malloc(nbytes);
  FILE *fin=open_file(arg[1],"r");
  size_t nbytes_read=fread(data,1,nbytes,fin);
  if(nbytes_read!=nbytes) crash("error reading file, obtained %lu instead of %lu bytes",nbytes_read,nbytes);
  fclose(fin);
  
  if(!check_endianess()) doubles_to_doubles_changing_endianess(data,data,ndoubles);

  FILE *fout=open_file(arg[2],"w");
  for(size_t iel=0;iel<nel;iel++)
    {
      for(size_t ijack=0;ijack<=njack;ijack++) fprintf(fout,"%16.16lg ",data[iel*(njack+1)+ijack]);
      fprintf(fout,"\n");
    }
  fclose(fout);
  
  return 0;
}
