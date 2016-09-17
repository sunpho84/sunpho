#include "include.h"

void read(int &nconfs,double *&ll_dag,double *&ll,double *&l,int L,const char *path)
{
  //define pars
  int LH=L/2,VH=(LH+1)*(LH+1)*(LH+1);
  int ll_nel=2*VH;
  int l_nel=2;
  int tot_nel=l_nel+ll_nel;
  
  //get file size
  FILE *fin=open_file("Luppoli","r");
  fseek(fin,0,SEEK_END);
  int rc=ftell(fin);
  rewind(fin);
  
  //count nconfs
  int tot_size=2*(tot_nel*sizeof(double)+4*sizeof(int));
  if(rc%tot_size!=0) crash("file size %d is not a multiple of size_per_conf = %d",rc,tot_size);
  nconfs=rc/tot_size;
  cout<<"Nconfs: "<<nconfs<<endl;
  
  //allocate ll, ll_dag and l
  ll=new double[nconfs*ll_nel];
  ll_dag=new double[nconfs*ll_nel];
  l=new double[nconfs*l_nel];
  
  //scan it
  for(int iconf=0;iconf<nconfs;iconf++)
    for(int il=0;il<2;il++)
      {
	int temp_int[4];
	double *ll_list[2]={ll_dag,ll};
	if((int)fread(temp_int,sizeof(int),4,fin)!=4) crash("loading 4 int for conf=%d, il=%d",iconf,il);
	if((int)fread(l+2*iconf,sizeof(double),2,fin)!=2) crash("poly for conf=%d, il=%d",iconf,il);
	if((int)fread(ll_list[il]+ll_nel*iconf,sizeof(double),ll_nel,fin)!=ll_nel) crash("loading poly_dag_corr iconf %d",iconf);
      }
  
  //check eof
  char c;
  if(fread(&c,1,1,fin)==1) crash("should not have read...");
}
