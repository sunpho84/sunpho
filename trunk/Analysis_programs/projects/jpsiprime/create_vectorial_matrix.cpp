#include "include.h"

void read_file(double *out,const char *path_in,int nel)
{
  FILE *fin=open_file(path_in,"r");
  int stat=fread((void*)out,sizeof(double),nel,fin);
  if(stat!=nel) crash("obtained %d while reading",stat);
  fclose(fin);
}

void write_file(const char *path_out,double *in,const char *mode,int nel)
{
  FILE *fout=open_file(path_out,mode);
  int stat=fwrite((void*)in,sizeof(double),nel,fout);
  if(stat!=nel) crash("obtained %d while writing",stat);
  fclose(fout);
}

int main()
{
  int fs=get_file_size("2pts_TKTK");
  int nel=fs/sizeof(double);
  
  double TKTK[nel];
  double TKVK[nel];
  double VKVK[nel];
  
  read_file(VKVK,"2pts_VKVK",nel);
  read_file(TKVK,"2pts_TKVK",nel);
  read_file(TKTK,"2pts_TKTK",nel);
  
  //change signs
  for(int iel=0;iel<nel;iel++)
    {
      TKTK[iel]*=-1;
      if(iel%48<24) TKVK[iel]*=-1;
    }
  
  int njacks=16;
  int T=48;
  int ntheta=2;
  int nr=4;
  int nri=2;
  int ndoub_per_combo=T*(njacks+1)*ntheta*nr*nri;
  int ncombo_lev=nel/ndoub_per_combo;
  int nlev=sqrt(ncombo_lev);
  if(nlev*nlev*ndoub_per_combo!=nel) crash("check");
  cout<<nlev<<endl;
  
  char out_path[]="2pts_VEVE";
  FILE *fout=open_file(out_path,"w");
  for(int itv1=0;itv1<2;itv1++)
    for(int ilev1=0;ilev1<nlev;ilev1++)
      for(int itv2=0;itv2<2;itv2++)
	for(int ilev2=0;ilev2<nlev;ilev2++)
	  {
	    double *b[4]={VKVK,TKVK,TKVK,TKTK};
	    
	    int itvc=itv1*2+itv2;
	    
	    for(int iel=0;iel<ndoub_per_combo;iel++)
	      {
		int iin=(ilev1*2+ilev2)*ndoub_per_combo+iel;
		double o=b[itvc][iin];
		
		int stat=fwrite((void*)&o,sizeof(double),1,fout);
		if(stat!=1) crash("obtained %d while writing",stat);
	      }
	  }

  return 0;
}
