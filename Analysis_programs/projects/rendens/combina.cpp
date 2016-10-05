#include "include.h"
#include <complex>

typedef complex<double> dcompl;

const double N_t=8;
const int ncopies=256;
const int nflavs=3;
//const int nconfs=124;
const int nconfs_poss=176;
const int njacks=20;
const int clust_size=nconfs_poss/njacks;
const int nconfs=clust_size*njacks;
const int Ns=32,Nt=8;
const int V4=Ns*Ns*Ns/Nt;
FILE *fin;

//all base traces
const int nbase_traces=9;

//name them
enum ibase_trace_t{M_dM,M_d2M,M_dM_M_dM,M_d3M,M_dM_M_d2M,M_dM_M_dM_M_dM,M_d2M_M_d2M,M_d2M_M_dM_M_dM,M_dM_M_dM_M_dM_M_dM};

//give access to the traces
inline int ind(int iflav,int itrace,int icopy)
{return icopy+ncopies*(itrace+nbase_traces*iflav);}

//hold data
const int nentr=ind(nflavs-1,nbase_traces-1,ncopies-1)+1;
vector<dcompl> base_Tr(nentr);

//read all data
void read_conf()
{
  for(int icopy=0;icopy<ncopies;icopy++)
    {
      //read the label
      int conf_lab;
      if(fscanf(fin,"%d",&conf_lab)!=1) crash("reading iconf: %d for copy %d",conf_lab,icopy);
      
      for(int iflav=0;iflav<nflavs;iflav++)
	for(int itrace=0;itrace<nbase_traces;itrace++)
	  {
	    double re,im;
	    if(fscanf(fin,"%lg %lg",&re,&im)!=2) crash("reading icopy %d, iflav %d, itrace %d",icopy,iflav,itrace);
	    base_Tr[ind(iflav,itrace,icopy)]=dcompl(re,im);
	  }
    }
}

//compute the product of traces summed over all possible copies
dcompl Tr_noperm(vector<int> itrace_list,int iflav)
{
  //Store the partial sums.
  //The last element will contain the full sum.
  //All elements are initialized to 1 so that
  //at first round the product is trivially the
  //sum of all traces
  vector<dcompl> P(ncopies+1,1);
  
  //loop over all traces inside the product
  for(size_t itr=0;itr<itrace_list.size();itr++)
    {
      //combinatorial factor
      double fact=(itr+1.0)/(ncopies-itr);
      
      //init to 0 the partial sum of the first 0 copies,
      //after taking note of the previous value
      dcompl prev=P[0];
      P[0]=0;
      
      //Put the result in a temporary storage, take note of old
      //partial sum of icopy, and replace it with the result, to avoid
      //using different vectors for partial sums
      for(int icopy=1;icopy<=ncopies;icopy++)
	{
	  dcompl res=P[icopy-1]+base_Tr[ind(iflav,itrace_list[itr],icopy-1)]*prev*fact;
	  prev=P[icopy];
	  P[icopy]=res;
	}
    }
  
  //return the last elment, total sum
  return P[ncopies];
}

//compute the triangular partial sum of the list of traces passed
dcompl Tr(vector<int> itrace_list,int iflav)
{
  //result averaged over all permutations
  dcompl res=0;
  
  //sort the list of indices of traces so that the algorithm that
  //generates all permutations find them in the correct order and stop
  //appropriately
  int nperm=0;
  sort(itrace_list.begin(),itrace_list.end());
  
  //loop over all permutations
  do
    {
      res+=Tr_noperm(itrace_list,iflav);
      nperm++;
    }
  while(next_permutation(itrace_list.begin(),itrace_list.end()));
  
  //return the average over all permutations
  return res/(double)nperm;
}

//return the number of permutations
int Tr_get_nperm(vector<int> itrace_list)
{
  int nperm=0;
  sort(itrace_list.begin(),itrace_list.end());
  
  do nperm++;
  while(next_permutation(itrace_list.begin(),itrace_list.end()));
  
  return nperm;
}

//return the multiplicity
int Tr_get_mult(vector<int> itrace_list)
{
  int fact=1;
  for(size_t itr=0;itr<itrace_list.size();itr++) fact*=itr+1.0;
  
  return fact;
}

//name of the full trace
enum ifull_trace_t{
  Tr_M_dM,
  Tr_M_d2M,
  Tr_M_dM_M_dM,
  Tr_M_d3M,
  Tr_M_dM_M_d2M,
  Tr_M_dM_M_dM_M_dM,
  Tr_M_d2M_M_d2M,
  Tr_M_d2M_M_dM_M_dM,
  Tr_M_dM_M_dM_M_dM_M_dM,
  Tr_M_dM_Tr_M_dM,
  Tr_M_dM_Tr_M_dM_M_dM,
  Tr_M_dM_Tr_M_d2M,
  Tr_M_dM_Tr_M_dM_Tr_M_dM,
  Tr_M_dM_Tr_M_dM_M_d2M,
  Tr_M_d2M_Tr_M_dM_M_dM,
  Tr_M_dM_M_dM_Tr_M_dM_M_dM,
  Tr_M_dM_M_dM_M_dM_Tr_M_dM,
  Tr_M_d2M_Tr_M_d2M,
  Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM,
  Tr_M_d2M_Tr_M_dM_Tr_M_dM,
  Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM};
int nfull_Tr=Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM+1;

//prepare a map
#define _CONCAT(X,Y) X##Y
#define _STR(X) #X
#define STR(X) _STR(X)
#define CONCAT(X,Y) _CONCAT(X,Y)
#define CONCAT2(TR1,TR2) CONCAT(TR1,TR2)
#define CONCAT_(TR1,TR2) CONCAT(CONCAT(TR1,_),TR2)
#define NAME_TR_1(TR) CONCAT2(Tr_,TR)
#define NAME_TR_2(TR1,TR2) CONCAT_(NAME_TR_1(TR1),NAME_TR_1(TR2))
#define NAME_TR_3(TR1,TR2,TR3) CONCAT_(NAME_TR_1(TR1),NAME_TR_2(TR2,TR3))
#define NAME_TR_4(TR1,TR2,TR3,TR4) CONCAT_(NAME_TR_1(TR1),NAME_TR_3(TR2,TR3,TR4))
#define DEF_MAP_1(TR1) full_Tr_map[NAME_TR_1(TR1)]=make_full_Tr_map({TR1},STR(NAME_TR_1(TR1)))
#define DEF_MAP_2(TR1,TR2) full_Tr_map[NAME_TR_2(TR1,TR2)]=make_full_Tr_map({TR1,TR2},STR(NAME_TR_2(TR1,TR2)))
#define DEF_MAP_3(TR1,TR2,TR3) full_Tr_map[NAME_TR_3(TR1,TR2,TR3)]=make_full_Tr_map({TR1,TR2,TR3},STR(NAME_TR_3(TR1,TR2,TR3)))
#define DEF_MAP_4(TR1,TR2,TR3,TR4) full_Tr_map[NAME_TR_4(TR1,TR2,TR3,TR4)]=make_full_Tr_map({TR1,TR2,TR3,TR4},STR(NAME_TR_4(TR1,TR2,TR3,TR4)))
inline pair<vector<int>,string> make_full_Tr_map(vector<int> map,string name)
{return make_pair(map,name);}

//host the full traces
jvec full_Tr(nflavs*nfull_Tr,njacks);
int ind_full(int iflav,int ifull)
{return ifull+nfull_Tr*iflav;}

//temporarily fix the flavor
int iflav_glb=0;
inline jack fTr(int ifull)
{return full_Tr[ind_full(iflav_glb,ifull)];}

int main()
{
  cout<<"Effective confs: "<<nconfs<<endl;
  
  fin=open_file("rende_new","r");
  vector<pair<vector<int>,string>> full_Tr_map(nfull_Tr);
  
  //prepare all maps
  DEF_MAP_1(M_dM);
  DEF_MAP_1(M_d2M);
  DEF_MAP_1(M_dM_M_dM);
  DEF_MAP_1(M_d3M);
  DEF_MAP_1(M_dM_M_d2M);
  DEF_MAP_1(M_dM_M_dM_M_dM);
  DEF_MAP_1(M_d2M_M_d2M);
  DEF_MAP_1(M_d2M_M_dM_M_dM);
  DEF_MAP_1(M_dM_M_dM_M_dM_M_dM);
  DEF_MAP_2(M_dM,M_dM);
  DEF_MAP_2(M_dM,M_dM_M_dM);
  DEF_MAP_2(M_dM,M_d2M);
  DEF_MAP_3(M_dM,M_dM,M_dM);
  DEF_MAP_2(M_dM,M_dM_M_d2M);
  DEF_MAP_2(M_d2M,M_dM_M_dM);
  DEF_MAP_2(M_dM_M_dM,M_dM_M_dM);
  DEF_MAP_2(M_dM_M_dM_M_dM,M_dM);
  DEF_MAP_2(M_d2M,M_d2M);
  DEF_MAP_3(M_dM_M_dM,M_dM,M_dM);
  DEF_MAP_3(M_d2M,M_dM,M_dM);
  DEF_MAP_4(M_dM,M_dM,M_dM,M_dM);
  
  //read all flavour, confs and take all trace products conf by conf
  for(int iconf=0;iconf<nconfs;iconf++)
    {
      const int ijack=iconf/clust_size;
      
      read_conf();
      
      for(int iflav=0;iflav<nflavs;iflav++)
  	for(int ifull=0;ifull<nfull_Tr;ifull++)
  	  full_Tr[ind_full(iflav,ifull)][ijack]+=
  	    Tr(full_Tr_map[ifull].first,iflav).real();
      
  //     dcompl a=0;
  //     double nc=0;
  //     for(int icopy=0;icopy<ncopies;icopy++)
  //     	for(int jcopy=icopy+1;jcopy<ncopies;jcopy++)
  //     	   for(int kcopy=jcopy+1;kcopy<ncopies;kcopy++)
  // 	     for(int lcopy=kcopy+1;lcopy<ncopies;lcopy++)
  // 	     {
  // 	       a+=
  // 		 base_Tr[ind(0,M_dM,icopy)]*
  // 		 base_Tr[ind(0,M_dM,jcopy)]*
  // 		 base_Tr[ind(0,M_dM,kcopy)]*
  // 		 base_Tr[ind(0,M_dM,lcopy)];
  // 	       nc+=1;
  //     	      }
      
  //     double ex=a.real()/nc;
  //     double al=full_Tr[ind_full(0,Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM)][iconf];
  //     cout<<"exact: "<<ex<<", algo: "<<al<<", rel diff: "<<(ex-al)/(ex+al)<<", nc: "<<nc<<" "<<((1.0)*(nconfs)*(ncopies)*(ncopies-1)*(ncopies-2)*(ncopies-3))<<endl;
    }
  
  //close and clusterize
  fclose(fin);
  full_Tr.clusterize(clust_size);
  
  for(int ifull=0;ifull<nfull_Tr;ifull++)
    cout<<full_Tr_map[ifull].second<<" ("<<Tr_get_nperm(full_Tr_map[ifull].first)<<" "<<Tr_get_mult(full_Tr_map[ifull].first)<<") = "<<fTr(ifull)<<endl;
  
  //compute susc2
  jack susc2_disc=(fTr(Tr_M_dM_Tr_M_dM)-sqr(fTr(Tr_M_dM)))/(16*V4);
  jack susc2_conn=(fTr(Tr_M_d2M)-fTr(Tr_M_dM_M_dM))/(4*V4);
  jack susc2_tot=susc2_conn+susc2_disc;
  
  jack susc4=
    -6*fTr(Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM)/64
    +1*fTr(Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM)/256
    +6*fTr(Tr_M_d2M_Tr_M_dM_Tr_M_dM)/64
    -3*fTr(Tr_M_dM_Tr_M_dM)*fTr(Tr_M_d2M)/64
    +3*fTr(Tr_M_dM_Tr_M_dM)*fTr(Tr_M_dM_M_dM)/64
    -3*fTr(Tr_M_dM_Tr_M_dM)*fTr(Tr_M_dM_Tr_M_dM)/256
    -6*fTr(Tr_M_d2M_Tr_M_dM_M_dM)/16
    +3*fTr(Tr_M_dM_M_dM_Tr_M_dM_M_dM)/16
    +1*fTr(Tr_M_dM_M_dM_M_dM_Tr_M_dM)/2
    +1*fTr(Tr_M_dM_Tr_M_dM)/16
    +3*fTr(Tr_M_d2M_Tr_M_d2M)/16
    -33*fTr(Tr_M_dM_Tr_M_dM_M_d2M)/64
    -3*fTr(Tr_M_d2M)*fTr(Tr_M_d2M)/16
    +3*fTr(Tr_M_d2M)*fTr(Tr_M_dM_M_dM)/16
    -3*fTr(Tr_M_d2M)*fTr(Tr_M_dM_Tr_M_dM)/64
    -3*fTr(Tr_M_d2M_M_d2M)/4
    +3*fTr(Tr_M_d2M_M_dM_M_dM)/1
    -1*fTr(Tr_M_dM_M_dM)/1
    -3*fTr(Tr_M_dM_M_dM_M_dM_M_dM)/2
    +1*fTr(Tr_M_d2M)/4
    +3*fTr(Tr_M_dM_M_dM)*fTr(Tr_M_d2M)/16
    -3*fTr(Tr_M_dM_M_dM)*fTr(Tr_M_dM_M_dM)/16
    +3*fTr(Tr_M_dM_M_dM)*fTr(Tr_M_dM_Tr_M_dM)/64;
  
  susc4/=V4*Nt*Nt;
  
  cout<<"susc2 = "<<susc2_tot<<endl;
  cout<<"susc4 = "<<susc4<<endl;
  
  return 0;
}
