#include "include.h"

const int nm=6;

struct header1_t
{
  int twist,nf,nsrc,T,L1,L2,L3,nk,nmoms;
};

struct header2_t
{
  double beta,ksea,musea,csw;
  double k[nm];
  double mu[nm];
  double mom[4];
};

struct header3_t
{
  int size;
};

ostream& operator<<(ostream &out,header1_t &h)
{
  out<<"twist: "<<h.twist<<endl;
  out<<"nf: "<<h.nf<<endl;
  out<<"nsrc: "<<h.nsrc<<endl;
  out<<"T: "<<h.T<<endl;
  out<<"L1: "<<h.L1<<endl;
  out<<"L2: "<<h.L2<<endl;
  out<<"L3: "<<h.L3<<endl;
  out<<"nk: "<<h.nk<<endl;
  out<<"nmoms: "<<h.nmoms<<endl;
  return out;
}

ostream& operator<<(ostream &out,header2_t &h)
{
  out<<"beta: "<<h.beta<<endl;
  out<<"ksea: "<<h.ksea<<endl;
  out<<"musea: "<<h.musea<<endl;
  out<<"csw: "<<h.csw<<endl;
  for(int im=0;im<nm;im++) out<<"k["<<im<<"]: "<<h.k[im]<<endl;
  for(int im=0;im<nm;im++) out<<"mu["<<im<<"]: "<<h.mu[im]<<endl;
  for(int mu=0;mu<4;mu++) out<<"mom["<<mu<<"]: "<<h.mom[mu]<<endl;
  return out;
}

ostream& operator<<(ostream &out,header3_t &h)
{
  out<<"size: "<<h.size<<endl;
  return out;
};

void read_file(const char *path)
{
  header1_t header1;
  header2_t header2;
  header3_t header3;
  
  FILE *fin=open_file(path,"r");
  
  if(fread(&header1,sizeof(header1_t),1,fin)!=1) crash("reading header1");
  if(fread(&header2,sizeof(header2_t),1,fin)!=1) crash("reading header2");
  if(fread(&header3,sizeof(header3_t),1,fin)!=1) crash("reading header3");
  
  cout<<header1<<header2<<header3<<endl;
  
  double *data=new double[header3.size];
  
  int T=48,njacks=15,nconfs=150;
  int clust_size=nconfs/njacks;
  jvec PP(T,njacks);
  PP=0;
  
  for(int iconf=0;iconf<nconfs;iconf++)
    {
      int jconf;
      if(fread(&jconf,sizeof(int),1,fin)!=1) crash("reading jconf");
      if(fread(data,sizeof(double),header3.size,fin)!=header3.size) crash("reading data");
      
      int iclust=iconf/clust_size;
      for(int t=0;t<T;t++) PP[t][iclust]+=data[2*t];
    }
  
  PP.clusterize(clust_size);
  
  cout<<-PP/header1.L1/header1.L2/header1.L3<<endl;
  
  delete[] data;
  
  fclose(fin);
}

int main(int narg,char **arg)
{
  if(narg<3) crash("use: %s filout filein",arg[0]);
  
  read_file(arg[2]);
  
  return 0;
}
