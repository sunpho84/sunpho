#include "include.h"

const int nflavs=3;
const int ndirs=1;
int T,L;

int main(int narg,char **arg)
{
  if(narg<8) crash("use %s file nconfs nterm njacks T L outfile",arg[0]);
  char *path=arg[1];
  int nconfs=atoi(arg[2]);
  int nterm=atoi(arg[3]);
  int njacks=atoi(arg[4]);
  T=atoi(arg[5]);
  L=atoi(arg[6]);
  const char *outpath=arg[7];
  int clust_size=(nconfs-nterm)/njacks;
  
  ifstream in(path);
  if(in.good()!=true) crash("opening");
  
  jvec data[ndirs*nflavs*(nflavs+1)/2];
  for(int icombo=0;icombo<ndirs*nflavs*(nflavs+1)/2;icombo++) data[icombo]=jvec(T,njacks)*0;
  
  for(int iconf=0;iconf<nconfs;iconf++)
    {
      int icombo=0;
      int ijack=(iconf-nterm)/clust_size;

      cout<<"iconf: "<<iconf<<", ijack: "<<ijack<<endl;
            
      for(int dir=0;dir<ndirs;dir++)
	for(int iflav=0;iflav<nflavs;iflav++)
	  for(int jflav=0;jflav<=iflav;jflav++)
	    {
	      for(int ih=0;ih<22;ih++)
		{
		  char ash[100];
		  if(!(in>>ash))
		    {
		      perror("boh");
		      crash("reading ash %d for conf %d",ih,iconf);
		    }
		}
	      
	      int NP;
	      if(dir==0) NP=T;
	      else       NP=L;
	      
	      if(iconf==0) data[icombo]=jvec(NP,njacks)*0;
	      
	      for(int t=0;t<NP;t++)
		{
		  int it;
		  double temp;
		  if(!(in>>it>>temp)||it!=t) crash("error while reading (icombo: %d, t: %d, it: %d, corr=%lg)",icombo,t,it,temp);
		  if(ijack>=0) data[icombo][t][ijack]+=temp;
		}
	      icombo++;
	    }
    }
  in.close();
  
  /////////////////////////////
  
  int icombo=0;
  for(int dir=0;dir<ndirs;dir++)
    for(int iflav=0;iflav<nflavs;iflav++)
      for(int jflav=0;jflav<=iflav;jflav++)
	{
	  data[icombo].clusterize(clust_size);
	  
	  ofstream fout(combine("/tmp/combo_flv_%d_flv_%d_dir_%d.xmg",iflav,jflav,dir));
	  fout<<"@type xydy"<<endl;
	  fout<<effective_mass(data[icombo].simmetrized(1));

	  if(icombo==0) data[icombo].write_to_binfile(outpath);
	  else data[icombo].append_to_binfile(outpath);
	  
	  icombo++;	  
	}  
  
  return 0;
}
