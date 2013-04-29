#include <TMatrixD.h>
#include <TVectorD.h>

jvec invert_square_root(int t,jvec *data,int nlevls)
{
  int L=data[0].nel-1;
  int njacks=data[0].njack;
  
  jvec sol(nlevls*nlevls,njacks);
  
  for(int ijack=0;ijack<=njacks;ijack++)
    if(t>=0&&t<=L)
      {
	TMatrixD toi(nlevls,nlevls);
	
	//upload
	for(int ism=0;ism<nlevls;ism++)
	  for(int jsm=0;jsm<nlevls;jsm++)
	    toi[ism][jsm]=data[ism*nlevls+jsm][t][ijack];
	
	TVectorD eigenValues;
	TMatrixD eigenVectors=toi.EigenVectors(eigenValues);
	
	//take the inverse square
	for(int ism=0;ism<nlevls;ism++)
	  for(int jsm=0;jsm<nlevls;jsm++)
	    {
	      sol[ism*nlevls+jsm].data[ijack]=0;
	      for(int ieig=0;ieig<nlevls;ieig++)
		sol[ism*nlevls+jsm].data[ijack]+=1/sqrt(eigenValues[ieig])*eigenVectors[ism][ieig]*eigenVectors[jsm][ieig];
	    }
    }
    else
      for(int ism=0;ism<nlevls;ism++)
	for(int jsm=0;jsm<nlevls;jsm++)
	  if(ism==jsm) sol[ism*nlevls+jsm]=1;
	  else         sol[ism*nlevls+jsm]=0;
  
  return sol;
}

/*
void diag(double *eigenValues,double *eigenVectors,double *ext_matr,int n)
{
  //transpose
  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
      eigenVectors[i+j*n]=ext_matr[i*n+j];
  
  //Query and allocate the optimal workspace
  int lwork=-1,lda=n,info;
  double wkopt;
  double *work;
  dsyev((char*)"Vectors",(char*)"Upper",&n,eigenVectors,&lda,eigenValues,&wkopt,&lwork,&info);
  lwork=(int)wkopt;
  work=(double*)malloc(lwork*sizeof(double));
  //Solve eigenproblem
  dsyev((char*)"Vectors",(char*)"Upper",&n,eigenVectors,&lda,eigenValues,work,&lwork,&info );
  //Check for convergence
  if(info>0) crash("The algorithm failed to compute eigenvalues");
  free(work);
}
*/

//find the eigenvalues
void find_eigens(int t,jvec &eig_va,jvec &eig_ve,jvec &inv,jvec *data,int nlevls)
{
  int njacks=data[0].njack;
  
  //create the operator
  jvec temp(nlevls*nlevls,njacks);
  for(int ism=0;ism<nlevls;ism++)
    for(int jsm=0;jsm<nlevls;jsm++)
      {
	temp[ism*nlevls+jsm]=0;
	for(int ksm=0;ksm<nlevls;ksm++)
	  temp[ism*nlevls+jsm]+=inv[ism*nlevls+ksm]*data[ksm*nlevls+jsm][t];
      }
  
  jvec op(nlevls*nlevls,njacks);
  for(int ism=0;ism<nlevls;ism++)
    for(int jsm=0;jsm<nlevls;jsm++)
      {
	op[ism*nlevls+jsm]=0;
	for(int ksm=0;ksm<nlevls;ksm++)
	  op[ism*nlevls+jsm]+=temp[ism*nlevls+ksm]*inv[ksm*nlevls+jsm];
      }

  for(int ijack=0;ijack<=njacks;ijack++)
    {
      TMatrixD toi(nlevls,nlevls);
      //double M[nlevls*nlevls];

      //set the matrix
      for(int ism_so=0;ism_so<nlevls;ism_so++)
	for(int ism_si=0;ism_si<nlevls;ism_si++)
	  /*M[ism_so*nlevls+ism_si]=*/toi[ism_so][ism_si]=op[ism_so*nlevls+ism_si][ijack];
      
      //hermitianize
      for(int ism_so=0;ism_so<nlevls;ism_so++)
	for(int ism_si=ism_so+1;ism_si<nlevls;ism_si++)
	  /*M[ism_so*nlevls+ism_si]=*/toi[ism_so][ism_si]=toi[ism_si][ism_so];
      
      //double E[nlevls],F[nlevls*nlevls];
      //diag(E,F,M,nlevls);
      
      //eigenvec
      TVectorD eigenValues;
      TMatrixD eigenVectors=toi.EigenVectors(eigenValues);
      //for(int i=0;i<nlevls;i++)cout<<i<<" "<<eigenValues[i]<<" "<<E[nlevls-1-i]<<endl;
      //download
      for(int ieig=0;ieig<nlevls;ieig++)
	{
	  eig_va[ieig].data[ijack]=/*E[ieig];/*/eigenValues[ieig];
	  for(int ilev=0;ilev<nlevls;ilev++)
	    {
	      eig_ve[ilev*nlevls+ieig].data[ijack]=eigenVectors[ilev][ieig];
	      if(eigenVectors[0][ieig]<0) eig_ve[ilev*nlevls+ieig].data[ijack]*=-1;
	    }
	}
    }
}
//find the matrix that diagonalize the corr
void find_diagonalizing_matrix(double *diag_m,int tinv,int t,jvec *data,int nlevls)
{
  int njacks=data[0].njack;
  
  //compute the inverse
  jvec inv=invert_square_root(tinv,data,nlevls);
  
  //allocate eigenvals
  jvec eig_va(nlevls,njacks);
  
  //allocate eigenvectors for each slice and diagonalize
  jvec eig_ve(nlevls*nlevls,njacks);
  find_eigens(t,eig_va,eig_ve,inv,data,nlevls);
  
  //take the product of inv*eig_ve, which is the matrix 
  //which diagonalize the original corr matrix
  for(int ilev=0;ilev<nlevls;ilev++)
    for(int ieig=0;ieig<nlevls;ieig++)
      {
	diag_m[ilev*nlevls+ieig]=0;
	for(int k=0;k<nlevls;k++)
	  diag_m[ilev*nlevls+ieig]+=(inv[ilev*nlevls+k]*eig_ve[k*nlevls+ieig]).med();
      }
      
  //normalize the matrix
  double norm(njacks);
  for(int ieig=0;ieig<nlevls;ieig++)
    {
      norm=0;
      for(int ilev=0;ilev<nlevls;ilev++)
	norm+=sqr(diag_m[ilev*nlevls+ieig]);
      norm=1/sqrt(norm);
      for(int ilev=0;ilev<nlevls;ilev++)
	diag_m[ilev*nlevls+ieig]*=norm;
    }
  
  //round the matrix
  /*
    for(int ieig=0;ieig<nlevls;ieig++)
    for(int ilev=0;ilev<nlevls;ilev++)
    diag_m[ilev*nlevls+ieig]=round(1000*diag_m[t1][ilev*nlevls+ieig])/1000;
  */
}
