#include "common.cpp"

double **corr_two_pts_fit,**err_two_pts_fit;
int *corr_parity;

template <class T> T fun_migrad_A_fit(T Za,T Zb,T M,double t,int par)
{
  if(par==1) return Za*Zb*exp(-M*TH)*cosh(M*(TH-t))/M;
  else       return Za*Zb*exp(-M*TH)*sinh(M*(TH-t))/M;
}

void chi2_stage_A(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  double M=p[nlevls];
  double *Z=p;
  
  for(int iop=0;iop<nlevls;iop++)
    for(int jop=0;jop<nlevls;jop++)
      for(int t=tfit_PP_min;t<=tfit_PP_max;t++)
          {
            double num=corr_two_pts_fit[iop*nlevls+jop][t];
            double teo=fun_migrad_A_fit(Z[iop],Z[jop],M,t,corr_parity[iop*nlevls+jop]);
            double diff=num-teo;
            double err=err_two_pts_fit[iop*nlevls+jop][t];
            double cont=sqr(diff/err);
            ch+=cont;
          }
}

void stage_A_fit(jack *M,jack *Z,jvec *corr,int *parity,int ist)
{
  //minimizator
  TMinuit minu;
  minu.SetPrintLevel(-1);
  minu.SetFCN(chi2_stage_A);

  //copy the data
  corr_parity=parity;
  corr_two_pts_fit=new double*[nlevls*nlevls];
  err_two_pts_fit=new double*[nlevls*nlevls];
  for(int iop=0;iop<nlevls;iop++)
    for(int jop=0;jop<nlevls;jop++)
      {
        err_two_pts_fit [iop*nlevls+jop]=new double[TH+1];
        corr_two_pts_fit[iop*nlevls+jop]=new double[TH+1];
        for(int t=0;t<=TH;t++) err_two_pts_fit[iop*nlevls+jop][t]=corr[iop*nlevls+jop][t].err();
      }
  
  //define the parameters
  for(int iop=0;iop<nlevls;iop++)
    {
      two_pts_fit(M[ist],Z[ist*nlevls+iop],corr[iop*nlevls+iop],tfit_PP_min,tfit_PP_max,
		  combine("plots/pre_fit_st_%d_m_%02d.xmg",ist,iop).c_str());
      Z[ist*nlevls+iop]=sqrt(Z[ist*nlevls+iop]);
      minu.DefineParameter(iop,combine("Z_%d",iop).c_str(),Z[ist*nlevls+iop].med(),Z[ist*nlevls+iop].err(),0,0);
    }
  
  minu.DefineParameter(nlevls,"M",M[ist].med(),M[ist].err(),0,0);
  
  for(int ijack_fit=0;ijack_fit<=njacks;ijack_fit++)
    {
      //copy jack
      for(int iop=0;iop<nlevls;iop++)
	for(int jop=0;jop<nlevls;jop++)
	  for(int t=0;t<=TH;t++)
	    corr_two_pts_fit[iop*nlevls+jop][t]=corr[iop*nlevls+jop][t][ijack_fit];
      
      //minimize
      minu.Migrad();
      
      //get back the parameters
      double dum;
      minu.GetParameter(nlevls,M[ist].data[ijack_fit],dum);
      for(int iop=0;iop<nlevls;iop++)
	minu.GetParameter(iop,Z[ist*nlevls+iop].data[ijack_fit],dum);
    }
  
  //get pars back
  double par[nlevls+1];
  for(int iop=0;iop<=nlevls;iop++)
    {
      double p,dum;
      minu.GetParameter(iop,p,dum);
      par[iop]=p;
    }
  
  //compute chi2
  double chi2;
  minu.Eval(nlevls+1,NULL,chi2,par,2);
  cout<<"M: "<<smart_print(M[ist])<<", Z: "<<smart_print(Z[ist*nlevls+0])<<", ch2: "<<chi2<<endl;
  
  //write the plots
  for(int iop=0;iop<nlevls;iop++)
    for(int jop=0;jop<nlevls;jop++)
      write_constant_fit_plot(combine("plots/fit_A_st_%d_op_%02d_%02d.xmg",ist,iop,jop).c_str(),
			      effective_mass(corr[iop*nlevls+jop],TH),M[ist],tfit_PP_min,tfit_PP_max);
}

int main(int narg,char **arg)
{
  if(narg<2) crash("use %s input" ,arg[0]);
  scan_input(arg[1]);
  
  jvec corr[nlevls][nlevls*nlevls];
  jack M[nlevls],Z[nlevls*nlevls];
  
  cout<<"Stage A"<<endl;
  for(int ilev=0;ilev<2;ilev++)
    {
      for(int iop=0;iop<nlevls;iop++)
        for(int jop=0;jop<nlevls;jop++)
          if(ilev==0) corr[ilev][iop*nlevls+jop]=g->data[iop*nlevls+jop];
          else
	    {
	      corr[ilev][iop*nlevls+jop]=corr[ilev-1][iop*nlevls+jop];
	      for(int t=0;t<=TH;t++)
		corr[ilev][iop*nlevls+jop][t]-=
		  fun_migrad_A_fit(Z[(ilev-1)*nlevls+iop],Z[nlevls*(ilev-1)+jop],M[ilev-1],t,g->parity[iop*nlevls+jop]);
	    }
      
      stage_A_fit(M,Z,corr[ilev],g->parity,ilev);
    }
  
  return 0;
}
