#include "include.h"
#include <complex>

typedef complex<double> dcompl;

const double N_t=8;
const int ncopies=256;
const int nflavs=3;
const int nconfs=6,njacks=nconfs;

namespace base
{
  //all base traces
  const int ntraces=9;
  
  //name them
  enum itrace_t{M_dM,M_d2M,M_dM_M_dM,M_d3M,M_dM_M_d2M,M_dM_M_dM_M_dM,M_d2M_M_d2M,M_d2M_M_dM_M_dM,M_dM_M_dM_M_dM_M_dM};
  
  //give access to the traces
  inline int ind(int iconf,int iflav,int itrace,int icopy)
  {return icopy+ncopies*(itrace+ntraces*(iflav+nflavs*iconf));}
  
  //hold data
  const int nentr=ind(nconfs-1,nflavs-1,ntraces-1,ncopies-1)+1;
  vector<dcompl> Tr(nentr);
  
  //read all data
  void read()
  {
    FILE *fin=open_file("rende_new","r");
    for(int iconf=0;iconf<nconfs;iconf++)
      for(int icopy=0;icopy<ncopies;icopy++)
	{
	  //read the label
	  int conf_lab;
	  if(fscanf(fin,"%d",&conf_lab)!=1) crash("reading iconf: %d",conf_lab);
	  
	  for(int iflav=0;iflav<nflavs;iflav++)
	    for(int itrace=0;itrace<base::ntraces;itrace++)
	      {
		double re,im;
		if(fscanf(fin,"%lg %lg",&re,&im)!=2) crash("reading iconf %d, icopy %d, iflav %d, itrace %d",iconf,icopy,iflav,itrace);
		Tr[ind(iconf,iflav,itrace,icopy)]=dcompl(re,im);
	      }
	}
  }
}

//compute the triangular partial sum of the list of traces passed
dcompl traces_calculator(vector<int> itrace_list,int iconf,int iflav)
{
  //store the partial sums
  vector<dcompl> P(ncopies+1,1);
  
  for(auto itrace : itrace_list)
    {
      //used to init to 1 the Partial sum of order 0
      dcompl prev=0;
      swap(prev,P[0]);
      
      //use "prev" as a temporary storage for the result
      //then swap it so the next round it will contain the previous partial sum
      for(int icopy=1;icopy<=ncopies;icopy++)
	{
	  prev=P[icopy-1]+base::Tr[base::ind(iconf,iflav,itrace,icopy-1)]*prev;
	  swap(prev,P[icopy]);
	}
    }
  
  return P[ncopies];
}

int main()
{
  int iconf=0,iflav=0;
  
  int ntraces=nbase_traces;
  vector<vector<ibase_tr>> itr_list(nbase_traces);
  
    
  dcompl Tr_M_dM=traces_calculator({M_dM},iconf,iflav);
  dcompl Tr_M_d2M=traces_calculator({M_d2M},iconf,iflav);
  dcompl Tr_M_dM_M_dM=traces_calculator({M_dM_M_dM},iconf,iflav);
  dcompl Tr_M_d3M=traces_calculator({M_d3M},iconf,iflav);
  dcompl Tr_M_dM_M_d2M=traces_calculator({M_dM_M_d2M},iconf,iflav);
  dcompl Tr_M_dM_M_dM_M_dM=traces_calculator({M_dM_M_dM_M_dM},iconf,iflav);
  
  dcompl Tr_M_dM2=
  
  dcompl d=0;
  for(int icopy=0;icopy<ncopies;icopy++)
    for(int jcopy=0;jcopy<icopy;jcopy++)
      d+=
	Tr[itr(iconf,iflav,0,icopy)]*
	Tr[itr(iconf,iflav,0,jcopy)];
  cout<<d<<endl;
  
  return 0;
}


//compute TR_M_dM
Tr_M_dM[iflav]+=     MAT_RES[icopy][iflav*nind_trace+0];

//compute TR_M_d2M
Tr_M_d2M[iflav]+=    MAT_RES[icopy][iflav*nind_trace+1];

//compute TR_M_dM_M_dM
Tr_M_dM_M_dM[iflav]+=MAT_RES[icopy][iflav*nind_trace+2];

//compute TR_M_dM_M_d2M
Tr_M_dM_M_d2M[iflav]+=MAT_RES[icopy][iflav*nind_trace+4];

//compute TR_M_dM_M_dM_M_dM
Tr_M_dM_M_dM_M_dM[iflav]+=MAT_RES[icopy][iflav*nind_trace+5];

//compute (TrM_dM)^2
traces_id2[0]=0;
traces_id2[1]=0;
Tr_M_dM2[iflav]+= traces_calculator(2,traces_id2, MAT_RES, iflav);


			//store the trace computed over one config. Will be used when computing the products of traces of different flavour matrix

		  //compute (TrM_dM)(Tr(M_dM)^2)
		  for(int icopy=0;icopy<ncopies_eff;icopy++)
		    for(int rcopy=icopy+1;rcopy<ncopies_eff;rcopy++)
		      {
			double complex complex1=MAT_RES[icopy][iflav*nind_trace+0];
			double complex complex2=MAT_RES[rcopy][iflav*nind_trace+2];
			Tr_M_dM_Tr_M_dM_M_dM[iflav]+=complex1*complex2;
			//compute the products for rcopy<icopy
			complex1=MAT_RES[rcopy][iflav*nind_trace+0];
			complex2=MAT_RES[icopy][iflav*nind_trace+2];
			Tr_M_dM_Tr_M_dM_M_dM[iflav]+=complex1*complex2;
		      }

//compute (TrM_dM)(TrM_d2M)
		  for(int icopy=0;icopy<ncopies_eff;icopy++)
		    for(int rcopy=icopy+1;rcopy<ncopies_eff;rcopy++)
		      {
			double complex complex1=MAT_RES[icopy][iflav*nind_trace+0];
			double complex complex2=MAT_RES[rcopy][iflav*nind_trace+1];
			Tr_M_dM_Tr_M_d2M[iflav]+=complex1*complex2;
			//compute the products for rcopy<ncopy
			complex1=MAT_RES[rcopy][iflav*nind_trace+0];
			complex2=MAT_RES[icopy][iflav*nind_trace+1];
			Tr_M_dM_Tr_M_d2M[iflav]+=complex1*complex2;

		      }

		  //compute (TrM_dM)^3

                  traces_id3[0]=0; traces_id3[1]=0; traces_id3[2]=0;
                  Tr_M_dM3[iflav]+= traces_calculator(3,traces_id3, MAT_RES, iflav);


		  //Traces at order 4

		  //compute Tr_M_d2M_M_d2M

		  for (int icopy=0;icopy<ncopies_eff;icopy++)
		    {
		      Tr_M_d2M_M_d2M[iflav] += MAT_RES[icopy][iflav*nind_trace + 6];
		    }

		  //compute Tr_M_d2M_M_dM_M_dM

		  for (int icopy=0;icopy<ncopies_eff;icopy++)
		    {
		      Tr_M_dM_M_dM_M_d2M[iflav] += MAT_RES[icopy][iflav*nind_trace + 7];
		    }

		  //compute Tr_M_dM_M_dM_M_dM_M_dM

		  for (int icopy=0;icopy<ncopies_eff;icopy++)
		    {
		      Tr_M_dM_M_dM_M_dM_M_dM[iflav] += MAT_RES[icopy][iflav*nind_trace + 8];
		    }



		  //        compute Tr_M_dM_Tr_M_dM_M_d2M

		  for(int icopy=0;icopy<ncopies_eff;icopy++)
		    for(int rcopy=icopy+1;rcopy<ncopies_eff;rcopy++)
		      {
			double complex complex1=MAT_RES[icopy][iflav*nind_trace+0];
			double complex complex2=MAT_RES[rcopy][iflav*nind_trace+4];
			Tr_M_d2M_Tr_M_dM_M_dM[iflav]+=complex1*complex2;
			//compute the products for rcopy<ncopy
			complex1=MAT_RES[rcopy][iflav*nind_trace+0];
			complex2=MAT_RES[icopy][iflav*nind_trace+4];
			Tr_M_d2M_Tr_M_dM_M_dM[iflav]+=complex1*complex2;

		      }



		  //compute Tr_M_d2M_Tr_M_dM_M_dM

		  for(int icopy=0;icopy<ncopies_eff;icopy++)
		    for(int rcopy=icopy+1;rcopy<ncopies_eff;rcopy++)
		      {
			double complex complex1=MAT_RES[icopy][iflav*nind_trace+1];
			double complex complex2=MAT_RES[rcopy][iflav*nind_trace+2];
			Tr_M_d2M_Tr_M_dM_M_dM[iflav]+=complex1*complex2;
			//compute the products for rcopy<ncopy
			complex1=MAT_RES[rcopy][iflav*nind_trace+1];
			complex2=MAT_RES[icopy][iflav*nind_trace+2];
			Tr_M_d2M_Tr_M_dM_M_dM[iflav]+=complex1*complex2;

		      }
		  //compute Tr_M_dM_M_dM_Tr_M_dM_M_dM

		  for(int icopy=0;icopy<ncopies_eff;icopy++)
		    for(int rcopy=icopy+1;rcopy<ncopies_eff;rcopy++)
		      {
			double complex complex1=MAT_RES[icopy][iflav*nind_trace+2];
			double complex complex2=MAT_RES[rcopy][iflav*nind_trace+2];
			Tr_M_dM_M_dM_Tr_M_dM_M_dM[iflav]+=complex1*complex2;
		      }

		  //compute Tr_M_dM_M_dM_M_dM_Tr_M_dM

		  for(int icopy=0;icopy<ncopies_eff;icopy++)
		    for(int rcopy=icopy+1;rcopy<ncopies_eff;rcopy++)
		      {
			double complex complex1=MAT_RES[icopy][iflav*nind_trace+5];
			double complex complex2=MAT_RES[rcopy][iflav*nind_trace+0];
			Tr_M_dM_M_dM_M_dM_Tr_M_dM[iflav]+=complex1*complex2;
			//compute the products for rcopy<ncopy
			complex1=MAT_RES[rcopy][iflav*nind_trace+5];
			complex2=MAT_RES[icopy][iflav*nind_trace+0];
			Tr_M_dM_M_dM_M_dM_Tr_M_dM[iflav]+=complex1*complex2;

		      }


		  //compute Tr_M_d2M_Tr_M_d2M

		  for(int icopy=0;icopy<ncopies_eff;icopy++)
		    for(int rcopy=icopy+1;rcopy<ncopies_eff;rcopy++)
		      {
			double complex complex1=MAT_RES[icopy][iflav*nind_trace+1];
			double complex complex2=MAT_RES[rcopy][iflav*nind_trace+1];
			Tr_M_d2M_Tr_M_d2M[iflav]+=complex1*complex2;
		      }

		  //compute Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM
                  traces_id3[0]=2; traces_id3[1]=0; traces_id3[2]=0;
                  Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[iflav]+= traces_calculator(3,traces_id3, MAT_RES, iflav);
                  traces_id3[0]=0; traces_id3[1]=2; traces_id3[2]=0;
                  Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[iflav]+= traces_calculator(3,traces_id3, MAT_RES, iflav);
                  traces_id3[0]=0; traces_id3[1]=0; traces_id3[2]=2;
                  Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[iflav]+= traces_calculator(3,traces_id3, MAT_RES, iflav);

		  //compute Tr_M_d2M_Tr_M_dM_Tr_M_dM


		  traces_id3[0]=0; traces_id3[1]=0; traces_id3[2]=1;

                  Tr_M_dM_Tr_M_dM_Tr_M_d2M[iflav]+= traces_calculator(3,traces_id3, MAT_RES, iflav);
                  traces_id3[0]=0; traces_id3[1]=1; traces_id3[2]=0;

                  Tr_M_dM_Tr_M_dM_Tr_M_d2M[iflav]+= traces_calculator(3,traces_id3, MAT_RES, iflav);
                  traces_id3[0]=1; traces_id3[1]=0; traces_id3[2]=0;

                  Tr_M_dM_Tr_M_dM_Tr_M_d2M[iflav]+= traces_calculator(3,traces_id3, MAT_RES, iflav);

		  //compute Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM


                  traces_id4[0]=0; traces_id4[1]=0; traces_id4[2]=0; traces_id4[3]=0;
                  Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM[iflav] += traces_calculator(4,traces_id4, MAT_RES, iflav);


