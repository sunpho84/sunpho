#include "include.h"

const int njacks=16;
int deb;
int T,TH,L,t0;
int parity[100];
int nlevls,nlevls_sto;
int reorder;
jvec *data;
char infile[100];


extern "C" {
  int dpotrf(char *uplo,int *n,double *a,int *lda,int *info);
  int dsygst(int *itype,char *uplo,int *n,double *a,int *lda,double *b,int *ldb,int *info);
  int dsygv(int *itype,char *jobz,char *uplo,int *n,double *a,int *lda,double *b,int *ldb,double *w,double *work,int *lwork,int *info);
  int dggev(char *jobz,char *uplo,int *n,double *a,int *lda,double *b,int *ldb,double *alpahr,double *alphai,double *beta,double *vl,int *nvl,double *vr,int *nvr,double *work,int *lwork,int *info);
}

void matr_copy(double *m,jvec *data,int t,int ijack)
{
  //copy jack of input matrix of t (nb: fortran order)
  for(int ilev=0;ilev<nlevls;ilev++)
    for(int jlev=0;jlev<nlevls;jlev++)
      m[jlev*nlevls+ilev]=data[ilev*nlevls+jlev][t][ijack];
}

void load_raw_data(const char *raw_data_path,int *map)
{
  //open output if path passed
  ofstream raw_data,eff_mass_raw_data;
  if(raw_data_path!=NULL)
    {
      raw_data.open(raw_data_path);
      raw_data<<"@type xydy"<<endl;
      
      eff_mass_raw_data.open(combine("eff_mass_%s",raw_data_path).c_str());
      eff_mass_raw_data<<"@type xydy"<<endl;
    }
  
  //load data and write the file with all the corrs used
  for(int ism_so=0;ism_so<nlevls;ism_so++)
    for(int ism_si=0;ism_si<nlevls;ism_si++)
      {
	//load and put the correct sign
	jvec temp=jvec_load(infile,T,njacks,map[ism_so]*nlevls_sto+map[ism_si]);
	if(temp[1].med()<0) temp=-temp;
	
	//infer correct parity
	parity[ism_so*nlevls+ism_si]=2*((temp[1]*temp[T-2]).med()>0)-1;
	data[ism_so*nlevls+ism_si]=temp.simmetrized(parity[ism_so*nlevls+ism_si]);
	
	//write raw data
	if(raw_data_path!=NULL)
	  if(ism_so==ism_si)
	    {
	      raw_data<<data[ism_so*nlevls+ism_si]<<"&"<<endl;
	      eff_mass_raw_data<<effective_mass(data[ism_so*nlevls+ism_si],TH,parity[ism_so*nlevls+ism_si]==1)
			       <<"&"<<endl;
	    }
      }
  
  //change correlator with the more precise
  int tch=14;
  for(int ism_so=0;ism_so<nlevls;ism_so++)
    for(int ism_si=0;ism_si<nlevls;ism_si++)
      {
	int i1=ism_so*nlevls+ism_si,i2=ism_si*nlevls+ism_so;
	if(data[i1][tch].err()>data[i2][tch].err()) data[i1]=data[i2];
	else data[i2]=data[i1];
    }
  
  //check the norm
  cout<<"Checking correct norm"<<endl;
  for(int ism_so=0;ism_so<nlevls;ism_so++)
    {
      for(int ism_si=0;ism_si<nlevls;ism_si++)
	cout<<smart_print(data[ism_so*nlevls+ism_si][tch]/
			  sqrt(data[ism_si*nlevls+ism_si][tch]*data[ism_so*nlevls+ism_so][tch]))<<" ";
      cout<<endl;
    }
  
  //check signularity
  for(int t=0;t<=TH;t++)
    {
      bool is_sing;

      jack sing(njacks);
      TMatrixD d(nlevls,nlevls);
      
      for(int ijack=0;ijack<=njacks;ijack++)
	{
	  for(int i=0;i<nlevls;i++)
	    for(int j=0;j<nlevls;j++)
	      d(i,j)=data[i*nlevls+j][t][ijack];
	  sing.data[ijack]=d.Determinant();
	}
      is_sing=(3*sing.err()>fabs(sing.med()));
      if(is_sing) cout<<"WARNING data matrix singular at: "<<t<<endl;
    }
  
  //close output
  if(raw_data_path!=NULL)
    {
      raw_data.close();
      eff_mass_raw_data.close();
    }
}

//solve the generalized eigenvalue problem
void gevp(jvec *eig_va,jvec *eig_ve,jvec &L,jvec *data,int nlevls,int t0)
{
  char low[]="L",VEC[]="V";//,no[]="N";
  
  for(int ijack=0,gijack=ijack;ijack<=njacks;ijack++)
    for(int t=0;t<=TH;t++)
      {
	//solve eigenproblem
	int itype=1,info;
	int lwork=(nlevls+2)*nlevls;
	double va[nlevls],work[lwork];
	
	//copy mt and mt0 in
	double mt[nlevls*nlevls],mt0[nlevls*nlevls];
	matr_copy(mt,data,t,ijack);
	matr_copy(mt0,data,t0,gijack);
	
	dsygv(&itype,VEC,low,&nlevls,mt,&nlevls,mt0,&nlevls,va,work,&lwork,&info);
	if(info) crash("info at t=%d ijack=%d after dsygv: %d",t,ijack,info);
	
	//copy eigenvalues, eigenvectors and L
	for(int ilev=0;ilev<nlevls;ilev++)
	  {
	    eig_va[ilev][t].data[ijack]=va[ilev];
	    for(int iop=0;iop<nlevls;iop++)
	      {
		//mt has each eigenvector in fortran-wise column (c row)
		L[ilev*nlevls+iop].data[ijack]=mt0[iop*nlevls+ilev];
		eig_ve[ilev*nlevls+iop][t].data[ijack]=mt[iop*nlevls+ilev];
	      }
	  }
      }
}

//compute the scalar product of vectors eig_ve[ilev1][t1] and eig_ve[ilev2][t2] w.r.t data[t0] metric
double scal_prod(jvec *eig_ve,jvec *data,int nlevls,int t0,int ijack0,int ilev1,int t1,int ijack1,
		 int ilev2,int t2,int ijack2)
{
  double r=0;
  for(int iop=0;iop<nlevls;iop++)
    for(int jop=0;jop<nlevls;jop++)
      r+=eig_ve[iop*nlevls+ilev1][t1][ijack1]*data[iop*nlevls+jop][t0][ijack0]*eig_ve[jop*nlevls+ilev2][t2][ijack2];
  
  return r;
}

//reorder eigenvalues and eigenvectors
void reorder_eig(jvec *eig_va,jvec *eig_ve,int nlevls)
{
  //fix comparison target
  int t2=t0+1;
  
  //compute ordering of ilev1(t) w.r.t ilev2(t-1)
  int flipped[TH+1];
  for(int t=0;t<=TH;t++) flipped[t]=0;
  
  for(int ijack1=njacks;ijack1>=0;ijack1--)
    {
      int ijack0=ijack1;
      int ijack2=ijack1;
      
      int imax[TH+1][nlevls];
      int sign[TH+1][nlevls];
      
      //check other times
      for(int t1=0;t1<=TH;t1++)
	{
	  //compute all scalar products
	  double r[nlevls*nlevls];
	  for(int ilev1=0;ilev1<nlevls;ilev1++)
	    for(int ilev2=0;ilev2<nlevls;ilev2++)
	      r[ilev1*nlevls+ilev2]=scal_prod(eig_ve,data,nlevls,t0,ijack0,ilev1,t1,ijack1,ilev2,t2,ijack2);
	  
	  //reset assignement and max
	  int assigned2[nlevls],assigned1[nlevls];
	  for(int ilev=0;ilev<nlevls;ilev++) assigned1[ilev]=assigned2[ilev]=0;
	  
	  //search the first three max
	  for(int ilev_ext=0;ilev_ext<nlevls;ilev_ext++)
	    {
	      if(deb) printf("Searching for %d max\n",ilev_ext);
	      double max_found=0;
	      int ilev1_found,ilev2_found;
	      
	      //scan the various levels
	      int found=0;
	      for(int ilev1=0;ilev1<nlevls;ilev1++)
		if(!assigned1[ilev1])
		  for(int ilev2=0;ilev2<nlevls;ilev2++)
		    if(!assigned2[ilev2])
		      {
			double locr=r[ilev1*nlevls+ilev2];
			if(deb)
			  printf("Scalar prod between vector %d (jack %d) at time %d and %d (jack %d) at time %d: %lg\n",
				 ilev1,ijack2,t1,ilev2,ijack2,t2,locr);
			
			//compare
			if(fabs(locr)>=fabs(max_found))
			  {
			    ilev2_found=ilev2;
			    ilev1_found=ilev1;
			    max_found=locr;
			    if(deb) printf(" new max found\n");
			    found=1;
			  }
			else if(deb) printf(" no max found\n");
		      }
	      
	      if(found==0) crash("no max found in any level!");
	      
	      //mark as assigned
	      assigned2[ilev2_found]=assigned1[ilev1_found]=1;
	      imax[t1][ilev1_found]=ilev2_found;
	      sign[t1][ilev1_found]=(max_found<0)?-1:+1;
	    }
	}

      //reorder time by time
      for(int t1=0;t1<=TH;t1++)
	{
	  //copy not to overwrite
	  double eig_va_te[nlevls],eig_ve_te[nlevls*nlevls];
	  for(int ilev1=0;ilev1<nlevls;ilev1++)
	    {
	      eig_va_te[ilev1]=eig_va[ilev1][t1][ijack1];
	      for(int ilev2=0;ilev2<nlevls;ilev2++) eig_ve_te[ilev1*nlevls+ilev2]=eig_ve[ilev1*nlevls+ilev2][t1][ijack1];
	    }
	  
	  //really reorder
	  for(int ilev1=0;ilev1<nlevls;ilev1++)
	    {
	      eig_va[ilev1][t1].data[ijack1]=eig_va_te[imax[t1][ilev1]];
	      for(int ilev2=0;ilev2<nlevls;ilev2++) 
		eig_ve[ilev1*nlevls+ilev2][t1].data[ijack1]=sign[t1][ilev2]*eig_ve_te[ilev1*nlevls+imax[t1][ilev2]];
	      
	      if(deb) printf("Lev %d most overlap with lev %d with sign %+d\n",ilev1,imax[t1][ilev1],sign[t1][ilev1]);
	      if(imax[t1][ilev1]!=ilev1)
		{
		  flipped[t1]++;
		  if(deb) printf("!!!!!! Time %d jack %d, ilev %d exchanged to %d\n",t1,ijack1,ilev1,imax[t1][ilev1]);
		}
	    }
	}
    }
  
  //count the number of flipping before and after t0
  int nflipped[2]={0,0};
  for(int t=1;t<=TH;t++)
    if(flipped[t]!=flipped[t-1]) nflipped[t>=t0]++;
  printf("NFlipped, before t0 %d, after or at t0: %d\n",nflipped[0],nflipped[1]);
}

int main(int narg,char **arg)
{
  if(narg<2) crash("use %s input" ,arg[0]);
  FILE *fin=open_file(arg[1],"r");
  
  //size
  read_formatted_from_file_expecting((char*)&T,fin,"%d","T");
  L=TH=T/2;
  
  //file path
  read_formatted_from_file_expecting(infile,fin,"%s","infile");
  
  //read the stored number of nlevls
  read_formatted_from_file_expecting((char*)&nlevls_sto,fin,"%d","nlevls_sto");

  //allocate nlevls-depending stuff
  read_formatted_from_file_expecting((char*)&nlevls,fin,"%d","nlevls");
  int map[nlevls];
  for(int ilev=0;ilev<nlevls;ilev++) read_formatted_from_file((char*)&(map[ilev]),fin,"%d","lev");
  data=(jvec*)calloc(nlevls*nlevls,sizeof(jvec));
  
  //timeslice for normalization
  read_formatted_from_file_expecting((char*)&t0,fin,"%d","t0");
  printf("Using t0: %d\n",t0);
  
  //read wheter reorder or not, and debug
  read_formatted_from_file_expecting((char*)&reorder,fin,"%d","reorder");
  printf("Reordering: %d\n",reorder);
  read_formatted_from_file_expecting((char*)&deb,fin,"%d","debug");
  printf("Debug: %d\n",deb);
  
  //interval to fit ground state
  int tfit_ground_min,tfit_ground_max;
  read_formatted_from_file_expecting((char*)&tfit_ground_min,fin,"%d","tfit_ground");
  read_formatted_from_file((char*)&tfit_ground_max,fin,"%d","tfit_ground");
  
  //interval to fit first state
  int tfit_first_min,tfit_first_max;
  read_formatted_from_file_expecting((char*)&tfit_first_min,fin,"%d","tfit_first");
  read_formatted_from_file((char*)&tfit_first_max,fin,"%d","tfit_first");
  
  fclose(fin);
  
  //load data
  load_raw_data("raw_data.xmg",map);
  
  ////////////////////////////// finished reading input ///////////////////////////
  
  //init
  jvec eig_va[nlevls],eig_ve[nlevls*nlevls];
  jvec L(nlevls*nlevls,njacks);
  for(int ilev=0;ilev<nlevls;ilev++)
    {
      eig_va[ilev]=jvec(TH+1,njacks);
      for(int iop=0;iop<nlevls;iop++) eig_ve[iop*nlevls+ilev]=jvec(TH+1,njacks);
    }
  
  //resolve gevp
  gevp(eig_va,eig_ve,L,data,nlevls,t0);
  
  //check orthogonality with respect to C(t0)  
  for(int t=0;t<=TH;t++)
    for(int ilev1=0;ilev1<nlevls;ilev1++)
      for(int ilev2=0;ilev2<nlevls;ilev2++)
	{
	  int ijack=0;
	  double r=scal_prod(eig_ve,data,nlevls,t0,ijack,ilev1,t,ijack,ilev2,t,ijack);
	  if(fabs(r-(ilev1==ilev2))>1.e-10) crash("r=%lg when expecting %d (err: %lg), t=%d",
						  r,(ilev1==ilev2),fabs(r-(ilev1==ilev2)),t);
	}
  
  //reorder
  if(reorder) reorder_eig(eig_va,eig_ve,nlevls);
  
  for(int t=0;t<TH;t++)
    {
      for(int ilev=0;ilev<nlevls;ilev++) cout<<smart_print(eig_va[ilev][t])<<" ";
      cout<<endl;
    }
  cout<<endl;
  
  {
    ofstream out("gevp.xmg");
    out<<"@type xydy"<<endl;
    for(int ilev=0;ilev<nlevls;ilev++)
      {
	jvec m=effective_mass(eig_va[ilev]);
	out<<m<<endl;
	out<<"&"<<endl;
      }
  }
  
  if(0)
  {
    ofstream out("/tmp/o.xmg");
    out<<"@type xydy"<<endl;
    out<<effective_mass(eig_va[1])/effective_mass(eig_va[0])<<endl;
  }
  
  return 0;
}
