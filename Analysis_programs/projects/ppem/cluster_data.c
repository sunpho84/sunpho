#include "include.h"

int njacks,L,TH,T;

void read_file(void *out,int nbytes,FILE *in)
{
  int nr=fread(out,1,nbytes,in);
  if(nr!=nbytes) crash("could not read %d bytes, obtained %d",nbytes,nr);
}

void write_file(FILE *out,int nbytes,void *in)
{
  int nw=fwrite(in,1,nbytes,out);
  if(nw!=nbytes) crash("could not write %d bytes, obtained %d",nbytes,nw);
}

int get_int(FILE *fin)
{
  int out;
  read_file(&out,sizeof(int),fin);
  
  return out;
}

double get_double(FILE *fin)
{
  double out;
  read_file(&out,sizeof(double),fin);
  
  return out;
}

void get_doubles(double *out,int ndoub,FILE *fin)
{read_file(out,ndoub*sizeof(double),fin);}

void put_doubles(FILE *fout,int ndoub,double *in)
{write_file(fout,ndoub*sizeof(double),in);}

int get_remaining_size(FILE *fin)
{
  int curr_pos=ftell(fin);
  int o1=fseek(fin,0,SEEK_END);
  int size=ftell(fin)-curr_pos;
  int o2=fseek(fin,curr_pos,SEEK_SET);
  if(o1||o2) crash("could not get file size");
  
  return size;
}

int main(int narg,char **arg)
{
  FILE *flist=open_file("input","r");
  read_formatted_from_file_expecting((char*)&L,flist,"%d","L");
  T=2*L;
  TH=L;
  read_formatted_from_file_expecting((char*)&njacks,flist,"%d","njacks");
  int nfiles;
  read_formatted_from_file_expecting((char*)&nfiles,flist,"%d","nfiles");
  
  for(int ifile=0;ifile<nfiles;ifile++)
    {
      char pathin[1024],pathout[1024];
      
      read_formatted_from_file(pathin,flist,"%s","pathin");
      read_formatted_from_file(pathout,flist,"%s","pathout");
      
      FILE *fin=open_file(pathin,"r");
      
      int twist=get_int(fin);
      printf("twist: %d\n",twist);
      
      int nf=get_int(fin);
      printf("nf: %d\n",nf);
      
      int nsrc=get_int(fin);
      printf("nsrc: %d\n\n",nsrc);
      
      int T=get_int(fin);
      printf("T: %d\n",T);
      
      int L1=get_int(fin);
      printf("L1: %d\n",L1);
      
      int L2=get_int(fin);
      printf("L2: %d\n",L2);
      
      int L3=get_int(fin);
      printf("L3: %d\n\n",L3);
      
      int nk=get_int(fin);
      printf("nk: %d\n",nk);
      
      int nmoms=get_int(fin);
      printf("nmoms: %d\n\n",nmoms);
      
      double beta=get_double(fin);
      printf("beta: %lg\n",beta);
      
      double ksea=get_double(fin);
      printf("ksea: %lg\n",ksea);
      
      double musea=get_double(fin);
      printf("musea: %lg\n",musea);
      
      double csw=get_double(fin);
      printf("csw: %lg\n\n",csw);
      
      for(int ik=0;ik<nk;ik++)
	{
	  double k=get_double(fin);
	  printf("k: %lg\n",k);
	}
      printf("\n");
      
      for(int ik=0;ik<nk;ik++)
	{
	  double mu_k=get_double(fin);
	  printf("mu_k: %lg\n",mu_k);
	}
      printf("\n");
      
      for(int imom=0;imom<nmoms;imom++)
	{
	  printf("mom_%d:",imom);
	  for(int idir=0;idir<4;idir++)
	    {
	      double mom=get_double(fin);
	      printf(" %lg",mom);
	    }
	  printf("\n\n");
	}
      
      int ndoubles=get_int(fin);
      printf("ndoubles: %d\n\n",ndoubles);
      int ndoubles_exp=2*nk*nk*nmoms*nmoms*T;
      int ndoubles_loop_exp=2*nk*nk*nk*nmoms*nmoms*T;
      if(ndoubles!=ndoubles_exp&&ndoubles!=ndoubles_loop_exp) crash("ndoubles read (%d) do not agree with expectation (%d or %d)",ndoubles,ndoubles_exp,ndoubles_loop_exp);
      
      int k3_max;
      if(ndoubles==ndoubles_loop_exp) k3_max=nk;
      else k3_max=1;
      
      double *data_read=(double*)malloc(ndoubles*sizeof(double));
      double *clusters[njacks];
      for(int iclust=0;iclust<njacks;iclust++)
	{
	  clusters[iclust]=(double*)malloc(ndoubles*sizeof(double));
	  memset(clusters[iclust],0,sizeof(double)*ndoubles);
	}
      int nconfs_teo=get_remaining_size(fin)/sizeof(double)/ndoubles;
      printf("File contains %d confs\n",nconfs_teo);
      int clust_size=nconfs_teo/njacks;
      printf("Cluster size: %d\n",clust_size);
      int nconfs=clust_size*njacks;
      printf("We will use %d confs\n",nconfs);
      
      for(int iconf=0;iconf<nconfs;iconf++)
	{
	  int iclust=iconf/clust_size;
	  
	  int id_conf=get_int(fin);
	  printf("id_conf: %d, clust: %d\n",id_conf,iclust);
	  
	  get_doubles(data_read,ndoubles,fin);
	  int pos=0;
	  for(int k3=0;k3<k3_max;k3++)
	    for(int k2=0;k2<nk;k2++)
	      for(int k1=0;k1<nk;k1++)
		for(int imom2=0;imom2<nmoms;imom2++)
		  for(int imom1=0;imom1<nmoms;imom1++)
		    for(int t=0;t<T;t++)
		      for(int ri=0;ri<2;ri++)
			{
			  clusters[iclust][pos]+=data_read[pos];
			  pos++;
			}
	}
      
      fclose(fin);
      
      FILE *fout=open_file(pathout,"w");
      
      for(int k3=0;k3<k3_max;k3++)
	for(int k2=0;k2<nk;k2++)
	  for(int k1=0;k1<nk;k1++)
	    for(int imom2=0;imom2<nmoms;imom2++)
	      for(int imom1=0;imom1<nmoms;imom1++)
		for(int ri=0;ri<2;ri++)
		  for(int t=0;t<T;t++)
		    {
		      int pos_in=ri+2*(t+T*(imom1+nmoms*(imom2+nmoms*(k1+nk*(k2+nk*k3)))));
		      double jacks[njacks+1];
		      jacks[njacks]=0;
		      for(int ijack=0;ijack<njacks;ijack++)
			{
			  jacks[ijack]=clusters[ijack][pos_in];
			  jacks[njacks]+=clusters[ijack][pos_in];
			}
		      for(int ijack=0;ijack<njacks;ijack++) jacks[ijack]=(jacks[njacks]-jacks[ijack])/(nconfs-clust_size);
		      jacks[njacks]/=nconfs;
		      
		      put_doubles(fout,njacks+1,jacks);
		      
		      if(ri==0&&imom2==0&&imom2==0&&k1==0&&k2==0&&k3==0) printf("%lg %d %d %d\n",jacks[njacks],t,pos_in,k3_max);
		    }
      
      fclose(fout);
  
      for(int iclust=0;iclust<njacks;iclust++) free(clusters[iclust]);
      free(data_read);
    }

  return 0;
}
