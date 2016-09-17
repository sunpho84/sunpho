#include "include.h"

int njack=15;

char mess_pt1[]="               GROUPING: jacknife\n               VALUES OF K_WILSON: 0000 0000\n               CORRELATION(S): cor2\n               FROM CONF. 1001 TO CONF. 1016\n               BY CLUSTERING OR DECIMATING    1 CONFS\n  MA_P5P5           ";

char mess_pt2[]="  ZA_P5P5           ";

int remap(int ijack)
{
  if(ijack==0) return njack;
  else return ijack-1;
}

int main(int narg,char **arg)
{
  if(narg<2) crash("use: %s input_file",arg[0]);
  
  FILE *in=open_file(arg[1],"r");
  
  int nens;
  read_formatted_from_file_expecting((char*)&nens,in,"%d","nens");
  
  FILE *list_file=open_file("lista","w");
  
  for(int iens=0;iens<nens;iens++)
    {
      char path_ens_dir[100];
      read_formatted_from_file(path_ens_dir,in,"%s","path_ens");
      int ibeta;
      read_formatted_from_file((char*)&ibeta,in,"%d","ibeta");
      int L;
      read_formatted_from_file((char*)&L,in,"%d","L");
      double am;
      read_formatted_from_file((char*)&am,in,"%lg","am");
      
      char path_ens[100];
      sprintf(path_ens,"%s/MZ2",path_ens_dir);
      
      jack M(njack),Z2(njack);
      M.load(path_ens,0);
      Z2.load(path_ens,1);
      cout<<"Opened file: "<<path_ens<<endl;
      
      string outdir=combine("%s/%d/%2d/%1.4f/","lattice_data",ibeta+1,L,am);
      system(combine("mkdir -p %s",outdir.c_str()).c_str());
      string outpath=outdir+"P5P5";
      
      FILE *fout=open_file(outpath.c_str(),"w");
      cout<<"Writing to file: "<<outpath<<endl;
      fprintf(fout,"%s\n",mess_pt1);
      for(int ijack=0;ijack<=njack;ijack++) fprintf(fout,"   %09.09lg\t\t%d\n",M.data[remap(ijack)],ijack);
      fprintf(fout,"%s\n",mess_pt2);
      for(int ijack=0;ijack<=njack;ijack++) fprintf(fout,"   %09.09lg\t\t%d\n",Z2.data[remap(ijack)],ijack);
      fclose(fout);
      
      fprintf(list_file,"1 %1.4f %2d %d\n%s\n",am,L,ibeta+1,outdir.c_str());
    }
  
  fclose(in);  
  
  return 0;
}
