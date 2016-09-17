#ifndef NOROOT
#include <TMatrixD.h>
#endif

const int DEBUG_GEVP=0;

extern "C" {
  int dpotrf(char *uplo,int *n,double *a,int *lda,int *info);
  int dsygst(int *itype,char *uplo,int *n,double *a,int *lda,double *b,int *ldb,int *info);
  int dsygv_(int *itype,char *jobz,char *uplo,int *n,double *a,int *lda,double *b,int *ldb,double *w,double *work,int *lwork,int *info);
  int dggev(char *jobz,char *uplo,int *n,double *a,int *lda,double *b,int *ldb,double *alpahr,double *alphai,double *beta,double *vl,int *nvl,double *vr,int *nvr,double *work,int *lwork,int *info);
  int dsyev_(char *jobz,char *uplo,int *n,double *a,int *lda,double *work,double *w,int *lwork,int *info);
}

class gevp_pars_t
{
public:
  void load_raw_data(const char *raw_data_path,const char *infile,int *map,int nlevls_sto,int off);
  gevp_pars_t(int nlevls,int njacks,int TH,int t0);
  ~gevp_pars_t();
  void reorder_eig();
  void check_norm(int tch);
  void check_singularity(int tch);
  void check_orthogonality();
  double scal_prod(int ijack0,int ilev1,int t1,int ijack1,int ilev2,int t2,int ijack2);
  void matr_copy(double *m,jvec *data,int t,int ijack);
  void gevp();
  void convert_to_full_eig_ve();
  void gevp_preliminary();
  jvec sqrt_data_t0;
  jvec inv_sqrt_data_t0;
  jvec *eig_va,*eig_ve,*full_eig_ve;
  jvec *data;
  int *parity;
  int nlevls;
  int njacks;
  int TH;
  int t0;
};

void gevp_pars_t::load_raw_data(const char *raw_data_path,const char *infile,int *map,int nlevls_sto,int off)
{
  const int T=TH*2;
  
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
  parity=new int[nlevls*nlevls];
  int *sign=new int[nlevls*nlevls];
  for(int ilev_so=0;ilev_so<nlevls;ilev_so++)
    for(int ilev_si=0;ilev_si<nlevls;ilev_si++)
      {
	int i=ilev_so*nlevls+ilev_si;
	//load and put the correct sign
	jvec temp=jvec_load(infile,T,njacks,map[ilev_so]*nlevls_sto+map[ilev_si]+off);
	sign[i]=1;
	if(temp[1].med()<0)
	  {
	    temp=-temp;
	    sign[i]=-1;
	  }
	
	//infer correct parity
	parity[i]=2*((temp[1]*temp[T-2]).med()>0)-1;
	data[i]=temp.simmetrized(parity[i]);
      } 

  cout<<"Parity, sign:"<<endl;
  for(int ilev_so=0;ilev_so<nlevls;ilev_so++)
    {
      cout<<" ";
      for(int ilev_si=0;ilev_si<nlevls;ilev_si++) cout<<parity[ilev_so*nlevls+ilev_si]<<" ";
      cout<<"\t";
      for(int ilev_si=0;ilev_si<nlevls;ilev_si++) cout<<sign[ilev_so*nlevls+ilev_si]<<" ";
      cout<<endl;
    }
  cout<<endl;

  //change correlator with the more precise
  int tch=14;
  for(int ilev_so=0;ilev_so<nlevls;ilev_so++)
    for(int ilev_si=0;ilev_si<nlevls;ilev_si++)
      {
	int i1=ilev_so*nlevls+ilev_si,i2=ilev_si*nlevls+ilev_so;
	if(data[i1][tch].err()>data[i2][tch].err()) data[i1]=data[i2];
	else data[i2]=data[i1];
      }
  
  //write raw data
  if(raw_data_path!=NULL)
    for(int ilev_so=0;ilev_so<nlevls;ilev_so++)
      for(int ilev_si=0;ilev_si<nlevls;ilev_si++)
	//if(ilev_so==ilev_si)
	{
	  raw_data<<data[ilev_so*nlevls+ilev_si]<<"&"<<endl;
	  eff_mass_raw_data<<effective_mass(data[ilev_so*nlevls+ilev_si],TH,
					    (parity[ilev_so*nlevls+ilev_si]==1))<<"&"<<endl;
	}

  //close output
  if(raw_data_path!=NULL)
    {
      raw_data.close();
      eff_mass_raw_data.close();
    }

  check_norm(tch);
  check_singularity(tch);
}

//check the norm
void gevp_pars_t::check_norm(int tch)
{
  cout<<"Checking correct norm"<<endl;
  for(int ilev_so=0;ilev_so<nlevls;ilev_so++)
    {
      for(int ilev_si=0;ilev_si<nlevls;ilev_si++)
	cout<<smart_print(data[ilev_so*nlevls+ilev_si][tch]/
			  sqrt(data[ilev_si*nlevls+ilev_si][tch]*data[ilev_so*nlevls+ilev_so][tch]))<<" ";
      cout<<endl;
    }
}

//check singularity
void gevp_pars_t::check_singularity(int tch)
{
#ifndef NOROOT
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
      is_sing=(3*sing.err()>=fabs(sing.med()));
      if(is_sing) cout<<"WARNING data matrix singular at: "<<t<<": "<<smart_print(sing)<<endl;
    }

#endif
}

gevp_pars_t::gevp_pars_t(int nlevls,int njacks,int TH,int t0) : nlevls(nlevls),njacks(njacks),TH(TH),t0(t0)
{
  eig_va=new jvec[nlevls];
  eig_ve=new jvec[nlevls*nlevls];
  full_eig_ve=new jvec[nlevls*nlevls];
  data=new jvec[nlevls*nlevls];
  sqrt_data_t0=jvec(nlevls*nlevls,njacks);
  inv_sqrt_data_t0=jvec(nlevls*nlevls,njacks);
  for(int ilev=0;ilev<nlevls;ilev++)
    {
      eig_va[ilev]=jvec(TH+1,njacks);
      for(int iop=0;iop<nlevls;iop++)
	{
	  eig_ve[iop*nlevls+ilev]=jvec(TH+1,njacks);
	  full_eig_ve[iop*nlevls+ilev]=jvec(TH+1,njacks);
	  data[iop*nlevls+ilev]=jvec(TH+1,njacks);
	}
    }
}

gevp_pars_t::~gevp_pars_t()
{
  delete [] eig_ve;
  delete [] eig_va;
  delete [] data;
}

void gevp_pars_t::matr_copy(double *m,jvec *d,int t,int ijack)
{
  //copy jack of input matrix of t (nb: fortran order)
  for(int ilev=0;ilev<nlevls;ilev++)
    for(int jlev=0;jlev<nlevls;jlev++)
      m[jlev*nlevls+ilev]=d[ilev*nlevls+jlev][t][ijack];
}

void gevp_pars_t::gevp_preliminary()
{
  char low[]="L",VEC[]="V"; 
  int info;
  int lwork=(nlevls+2)*nlevls;
  double va[nlevls],work[lwork];
  double eig_ve0[nlevls*nlevls],eig_va0[nlevls];
  
  for(int ijack=0;ijack<=njacks;ijack++)
    {
      //copy mt and mt0 in
      int ijack0=ijack;
      double mt0[nlevls*nlevls];
      matr_copy(mt0,data,t0,ijack0);
      
      dsyev_(VEC,low,&nlevls,mt0,&nlevls,va,work,&lwork,&info);
      if(info) crash("info at ijack=%d after dsyev: %d",ijack,info);

	//copy eigenvalues, eigenvectors
	for(int ilev=0;ilev<nlevls;ilev++)
	  {
	    eig_va0[ilev]=va[ilev];
	    for(int iop=0;iop<nlevls;iop++)
	      eig_ve0[ilev*nlevls+iop]=mt0[iop*nlevls+ilev];
	  }
	
	//reconstruct inverse square root
	for(int i=0;i<nlevls;i++)
	  for(int j=0;j<nlevls;j++)
	    {
	      sqrt_data_t0[i*nlevls+j].data[ijack]=0;
	      inv_sqrt_data_t0[i*nlevls+j].data[ijack]=0;
	      for(int k=0;k<nlevls;k++)
		{		
		  sqrt_data_t0[i*nlevls+j].data[ijack]+=sqrt(eig_va0[k])*eig_ve0[j*nlevls+k]*eig_ve0[i*nlevls+k];
		  inv_sqrt_data_t0[i*nlevls+j].data[ijack]+=1/sqrt(eig_va0[k])*eig_ve0[j*nlevls+k]*eig_ve0[i*nlevls+k];
		}
	    }
    }
  
  //test
  if(DEBUG_GEVP)
    {
      cout<<"test inverse square root"<<endl;
      for(int i=0;i<nlevls;i++)
	{
	  for(int l=0;l<nlevls;l++)
	    {
	      jack t(njacks);
	      t*=0;
	      for(int j=0;j<nlevls;j++)
		for(int k=0;k<nlevls;k++)
		  t+=inv_sqrt_data_t0[i*nlevls+j]*inv_sqrt_data_t0[j*nlevls+k]*data[k*nlevls+l][t0];
	      cout<<smart_print(t)<<" ";
	    }
	  cout<<endl;
	}
    }
}

//solve the generalized eigenvalue problem
void gevp_pars_t::gevp()
{
  gevp_preliminary();
  
  char low[]="L",VEC[]="V";//,no[]="N";
  for(int ijack=0;ijack<=njacks;ijack++)
    for(int t=0;t<=TH;t++)
      {
	//solve eigenproblem
	int info;
	int lwork=(nlevls+2)*nlevls;
	double va[nlevls],work[lwork];
	
	//copy mt and mt0 in
	int ijack0=ijack;
	double mt[nlevls*nlevls];

	for(int iop=0;iop<nlevls;iop++)
	  for(int lop=0;lop<nlevls;lop++)
	    {
	      mt[iop*nlevls+lop]=0;
	      for(int jop=0;jop<nlevls;jop++)
		for(int kop=0;kop<nlevls;kop++)
		  mt[iop*nlevls+lop]+=
		    inv_sqrt_data_t0[iop*nlevls+jop][ijack0]*
		    data[jop*nlevls+kop][t][ijack]*
		    inv_sqrt_data_t0[kop*nlevls+lop][ijack0];
	    }

	dsyev_(VEC,low,&nlevls,mt,&nlevls,va,work,&lwork,&info);
	if(info) crash("info at t=%d ijack=%d after dsyev: %d",t,ijack,info);
	
	//copy eigenvalues, eigenvectors
	for(int ilev=0;ilev<nlevls;ilev++)
	  {
	    eig_va[ilev][t].data[ijack]=va[ilev];
	    //mt has each eigenvector in fortran-wise column (c row)
	    //this is indeed how we want it
	    for(int iop=0;iop<nlevls;iop++)
	      eig_ve[iop*nlevls+ilev][t].data[ijack]=mt[ilev*nlevls+iop];
	  }
      }
}

//compute the scalar product of vectors eig_ve[ilev1][t1] and eig_ve[ilev2][t2]
double gevp_pars_t::scal_prod(int ijack0,int ilev1,int t1,int ijack1,int ilev2,int t2,int ijack2)
{
  double r=0;
  for(int iop=0;iop<nlevls;iop++)
      r+=eig_ve[iop*nlevls+ilev1][t1][ijack1]*
	eig_ve[iop*nlevls+ilev2][t2][ijack2];
  
  return r;
}

//reorder eigenvalues and eigenvectors
void gevp_pars_t::reorder_eig()
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
      
      int imax[(TH+1)*nlevls];
      int sign[(TH+1)*nlevls];
      
      //check other times
      for(int t1=0;t1<=TH;t1++)
	{
	  //compute all scalar products
	  double r[nlevls*nlevls];
	  for(int ilev1=0;ilev1<nlevls;ilev1++)
	    for(int ilev2=0;ilev2<nlevls;ilev2++)
	      r[ilev1*nlevls+ilev2]=scal_prod(ijack0,ilev1,t1,ijack1,ilev2,t2,ijack2);
	  
	  //reset assignement and max
	  int assigned2[nlevls],assigned1[nlevls];
	  for(int ilev=0;ilev<nlevls;ilev++) assigned1[ilev]=assigned2[ilev]=0;
	  
	  //search the first max
	  for(int ilev_ext=0;ilev_ext<nlevls;ilev_ext++)
	    {
	      if(DEBUG_GEVP) printf("Searching for %d max\n",ilev_ext);
	      double max_found=0;
	      int ilev1_found=0,ilev2_found=0;
	      
	      //scan the various levels
	      int found=0;
	      for(int ilev1=0;ilev1<nlevls;ilev1++)
		if(!assigned1[ilev1])
		  for(int ilev2=0;ilev2<nlevls;ilev2++)
		    if(!assigned2[ilev2])
		      {
			double locr=r[ilev1*nlevls+ilev2];
			if(DEBUG_GEVP)
			  printf("Scalar prod between vector %d (jack %d) at time %d and %d (jack %d) at time %d: %lg\n",
				 ilev1,ijack2,t1,ilev2,ijack2,t2,locr);
			
			//compare
			if(fabs(locr)>=fabs(max_found))
			  {
			    ilev2_found=ilev2;
			    ilev1_found=ilev1;
			    max_found=locr;
			    if(DEBUG_GEVP) printf(" new max found\n");
			    found=1;
			  }
			else if(DEBUG_GEVP) printf(" no max found\n");
		      }
	      
	      if(found==0) crash("no max found in any level!");
	      
	      //mark as assigned
	      assigned2[ilev2_found]=assigned1[ilev1_found]=1;
	      imax[t1*nlevls+ilev1_found]=ilev2_found;
	      sign[t1*nlevls+ilev1_found]=(max_found<0)?-1:+1;
	    }
	}

      //reorder time by time
      for(int t1=0;t1<=TH;t1++)
	{
	  //copy not to overwrite
	  double eig_va_te[nlevls],eig_ve_te[nlevls*nlevls];
	  for(int ilev=0;ilev<nlevls;ilev++)
	    {
	      eig_va_te[ilev]=eig_va[ilev][t1][ijack1];
	      for(int iop=0;iop<nlevls;iop++) eig_ve_te[iop*nlevls+ilev]=eig_ve[iop*nlevls+ilev][t1][ijack1];
	    }
	  
	  //really reorder
	  for(int ilev=0;ilev<nlevls;ilev++)
	    {
	      int ilev_max=imax[t1*nlevls+ilev];
	      eig_va[ilev][t1].data[ijack1]=eig_va_te[ilev_max];
	      for(int iop=0;iop<nlevls;iop++) 
		eig_ve[iop*nlevls+ilev][t1].data[ijack1]=sign[t1*nlevls+ilev_max]*eig_ve_te[iop*nlevls+ilev_max];
	      
	      if(DEBUG_GEVP) printf("Lev %d most overlap with lev %d with sign %+d\n",
				    ilev,ilev_max,sign[t1*nlevls+ilev_max]);
	      if(imax[t1*nlevls+ilev]!=ilev)
		{
		  flipped[t1]++;
		  if(DEBUG_GEVP) printf("!!!!!! Time %d jack %d, ilev %d exchanged to %d\n",t1,ijack1,ilev,ilev_max);
		}
	    }
	}
    }
  
  //count the number of flipping before and after t0
  int nflipped[2]={0,0};
  for(int t=1;t<=TH;t++)
    if(flipped[t]!=flipped[t-1]) nflipped[t>=t0]++;
  printf("NFlipped, before t0 %d, after or at t0: %d\n",nflipped[0],nflipped[1]);
  
  //perform global reordering according to eff mass
  std::vector<std::pair<double,int> > mord;
  for(int ilev=0;ilev<nlevls;ilev++)
    {
      double m=effective_mass(eig_va[ilev])[2].data[njacks];
      mord.push_back(std::make_pair(m,ilev));
    }
  std::sort(mord.begin(),mord.end());
  
  for(int t1=0;t1<=TH;t1++)
    for(int ijack1=0;ijack1<=njacks;ijack1++)
      {
	//copy not to overwrite
	double eig_va_te[nlevls],eig_ve_te[nlevls*nlevls];
	for(int ilev=0;ilev<nlevls;ilev++)
	  {
	    eig_va_te[ilev]=eig_va[ilev][t1][ijack1];
	    for(int iop=0;iop<nlevls;iop++) eig_ve_te[iop*nlevls+ilev]=eig_ve[iop*nlevls+ilev][t1][ijack1];
	  }
	
	//really reorder
	for(int ilev=0;ilev<nlevls;ilev++)
	  {
	    int ilev_dest=mord[ilev].second;
	    eig_va[ilev][t1].data[ijack1]=eig_va_te[ilev_dest];
	    for(int iop=0;iop<nlevls;iop++) 
	      eig_ve[iop*nlevls+ilev][t1].data[ijack1]=eig_ve_te[iop*nlevls+ilev_dest];
	  }
      }
}

//convert to the eigenvector of gevp
void gevp_pars_t::convert_to_full_eig_ve()
{
  //multiply by inv of sqrt(dacta[t0])
  for(int t=0;t<=TH;t++)
    for(int iop=0;iop<nlevls;iop++)
      for(int ilev=0;ilev<nlevls;ilev++)
	{
	  full_eig_ve[iop*nlevls+ilev][t]=0;
	  for(int lop=0;lop<nlevls;lop++)
	    full_eig_ve[iop*nlevls+ilev][t]+=
	      inv_sqrt_data_t0[iop*nlevls+lop]*
	      eig_ve[lop*nlevls+ilev][t];
	}

  //printing rebuilt
  if(DEBUG_GEVP)
    {
      ofstream test_gevp_reco("gevp_reco.xmg");
      for(int ilev=0;ilev<nlevls;ilev++)
	{
	  jvec a(TH+1,njacks);
	  a*=0;
	  for(int t=0;t<=TH;t++)
	    for(int iop=0;iop<nlevls;iop++)
	      for(int kop=0;kop<nlevls;kop++)
		a[t]+=
		  full_eig_ve[iop*nlevls+ilev][t]*
		  data[iop*nlevls+kop][t]*
		  full_eig_ve[kop*nlevls+ilev][t];
	  
	  test_gevp_reco<<effective_mass(a);
	  test_gevp_reco<<"&"<<endl;
	}
    }

  //printing full eigen
  if(DEBUG_GEVP)
    {
      for(int ilev=0;ilev<nlevls;ilev++)
	{
	  cout<<"testing full eigenvectors for lev "<<ilev<<endl;
	  for(int t=0;t<=TH;t++)
	    {
	      cout<<t<<" ";
	      for(int iop=0;iop<nlevls;iop++) cout<<smart_print(full_eig_ve[iop*nlevls+ilev][t])<<" ";
	      cout<<endl;
	    }
	}
    }
}

//check orthogonality with respect to C(t0)
void gevp_pars_t::check_orthogonality()
{
  for(int t=0;t<=TH;t++)
    for(int ilev1=0;ilev1<nlevls;ilev1++)
      for(int ilev2=0;ilev2<nlevls;ilev2++)
	{
	  int ijack0=njacks,ijack1=ijack0,ijack2=ijack0;
	  double r=scal_prod(ijack0,ilev1,t,ijack1,ilev2,t,ijack2);
	  if(fabs(r-(ilev1==ilev2))>1.e-10) crash("r=%lg when expecting %d (err: %lg), t=%d",
						  r,(ilev1==ilev2),fabs(r-(ilev1==ilev2)),t);
	}
}
