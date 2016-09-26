#include "include.h"
#include <complex>

typedef complex<double> dcompl;

const double N_t=8;
const int ncopies=256;
const int nflavs=3;
const int nconfs=6,njacks=nconfs;
FILE *fin;

//all base traces
const int nbase_traces=9;

//name them
enum itrace_t{M_dM,M_d2M,M_dM_M_dM,M_d3M,M_dM_M_d2M,M_dM_M_dM_M_dM,M_d2M_M_d2M,M_d2M_M_dM_M_dM,M_dM_M_dM_M_dM_M_dM};

//give access to the traces
inline int ind(int iflav,int itrace,int icopy)
{return icopy+ncopies*(itrace+nbase_traces*iflav);}

//hold data
const int nentr=ind(nflavs-1,nbase_traces-1,ncopies-1)+1;
vector<dcompl> base_Tr(nentr);

//read all data
void read_flav_per_conf()
{
  for(int icopy=0;icopy<ncopies;icopy++)
    {
      //read the label
      int conf_lab;
      if(fscanf(fin,"%d",&conf_lab)!=1) crash("reading iconf: %d",conf_lab);
      
      for(int iflav=0;iflav<nflavs;iflav++)
	for(int itrace=0;itrace<nbase_traces;itrace++)
	  {
	    double re,im;
	    if(fscanf(fin,"%lg %lg",&re,&im)!=2) crash("reading icopy %d, iflav %d, itrace %d",icopy,iflav,itrace);
	    base_Tr[ind(iflav,itrace,icopy)]=dcompl(re,im);
	  }
    }
}

//compute the triangular partial sum of the list of traces passed
dcompl Tr(vector<int> itrace_list,int iflav)
{
  dcompl res=0;
  int nperm=0;
  sort(itrace_list.begin(),itrace_list.end());
  
  do
    {
      //store the partial sums
      vector<dcompl> P(ncopies+1,1);
      
      for(size_t itr=0;itr<itrace_list.size();itr++)
	{
	  //used to init to 1 the Partial sum of order 0
	  dcompl prev=0;
	  swap(prev,P[0]);
	  
	  //use "prev" as a temporary storage for the result
	  //then swap it so the next round it will contain the previous partial sum
	  //put here the combinatorial factor
	  double fact=1.0/(ncopies-itr);
	  for(int icopy=1;icopy<=ncopies;icopy++)
	    {
	      prev=P[icopy-1]+base_Tr[ind(iflav,itrace_list[itr],icopy-1)]*prev*fact;
	      swap(prev,P[icopy]);
	    }
	}
      
      res+=P[ncopies];
      nperm++;
    }
  while(next_permutation(itrace_list.begin(),itrace_list.end()));
  
  return res/(double)nperm;
}

int main()
{
  int iflav=0;
  
  fin=open_file("rende_new","r");
  
  for(int iconf=0;iconf<nconfs;iconf++)
    {
      read_flav_per_conf();
      
      dcompl Tr_M_dM=Tr({M_dM},iflav);
      dcompl Tr_M_d2M=Tr({M_d2M},iflav);
      dcompl Tr_M_dM_M_dM=Tr({M_dM_M_dM},iflav);
      dcompl Tr_M_d3M=Tr({M_d3M},iflav);
      dcompl Tr_M_dM_M_d2M=Tr({M_dM_M_d2M},iflav);
      dcompl Tr_M_dM_M_dM_M_dM=Tr({M_dM_M_dM_M_dM},iflav);
      dcompl Tr_M_d2M_M_d2M=Tr({M_d2M_M_d2M},iflav);
      dcompl Tr_M_d2M_M_dM_M_dM=Tr({M_d2M_M_dM_M_dM},iflav);
      dcompl Tr_M_dM_M_dM_M_dM_M_dM=Tr({M_dM_M_dM_M_dM_M_dM},iflav);
      dcompl Tr_M_dM_Tr_M_dM=Tr({M_dM,M_dM},iflav);
      dcompl Tr_M_dM_Tr_M_dM_M_dM=Tr({M_dM,M_dM_M_dM},iflav);
      dcompl Tr_M_dM_Tr_M_dM_Tr_M_dM=Tr({M_dM,M_dM,M_dM},iflav);
      dcompl Tr_M_dM_Tr_M_dM_M_d2M=Tr({M_dM,M_dM_M_d2M},iflav);
      dcompl Tr_M_d2M_Tr_M_dM_M_dM=Tr({M_d2M,M_dM_M_dM},iflav);
      dcompl Tr_M_dM_M_dM_Tr_M_dM_M_dM=Tr({M_dM_M_dM,M_dM_M_dM},iflav);
      dcompl Tr_M_dM_M_dM_M_dM_Tr_M_dM=Tr({M_dM_M_dM_M_dM,M_dM},iflav);
      dcompl Tr_M_d2M_Tr_M_d2M=Tr({M_d2M,M_d2M},iflav);
      dcompl Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM=Tr({M_dM_M_dM,M_dM_M_dM,M_dM_M_dM},iflav);
      dcompl Tr_M_d2M_Tr_M_dM_Tr_M_dM=Tr({M_d2M,M_dM,M_dM},iflav);
      dcompl Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM=Tr({M_dM,M_dM,M_dM,M_dM},iflav);
      
      cout<<Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM<<endl;
    }
  
  fclose(fin);
  
  return 0;
}
