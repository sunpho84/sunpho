#include "../radiative_charmonium/radiative_common.cpp"

const int nbeta=4;
char **base_corrs_path,**ens_name;
int nens,*L,*ibeta,*nmass;
double *cmass,*lmass;
double *FSE_K_corr,*H1_corr,*H2_corr;
int approx_Pi_FSE,approx_K_FSE;

const double kappa=2.8373;

//read data
bvec dM2K,dM2Pi,MPi,MK,eps_gam;

bvec ml;
int ref_ml_beta[4]={-1,-1,-1,-1};

bvec par_res_fit;

double aem=1.0/137;

const char set_color[nbeta][1024]={"black","blue","red","green4"};
const char set_fill_color[nbeta][1024]={"grey","turquoise","yellow","green"};
const char set_symbol[nbeta][1024]={"square","circle","triangle up","triangle left"};
const char set_legend[nbeta][1024]={"\\xb\\0=3.80","\\xb\\0=3.90","\\xb\\0=4.05","\\xb\\0=4.20"};
const char set_legend_fm[nbeta][1024]={"a = 0.098 fm","a = 0.085 fm","a = 0.067 fm","a = 0.054 fm"};

int plot_iboot;
double ms=ms_phys_med;

double fun_fit(double A,double B,double C,double D,double ml,double ms,double a)
{
  double a1=a/lat[1][plot_iboot];
  
  double db0=db0_med;
  double qpif=4*M_PI*f0_med;
  double chi_ll=db0*ml/qpif;
  double chi_ls=db0*(ml+ms)/2/qpif;
  
  return 4*M_PI*aem*sqr(f0_med)*(
				 4*A/pow(f0_med,4)+
				 -8*A/pow(f0_med,4)*chi_ll*log(db0*ml)+
				 -(3+8*A/pow(f0_med,4))*chi_ls*log(db0*(ml+ms)/2)+
				 B*chi_ll+
				 C*chi_ls+
				 D*a1*a1
				 );
}

double fin_inf_diff(double C,int iens,double M,double dM2L,int PiK)
{
  double pref=aem/sqr(L[iens]*lat[ibeta[iens]][plot_iboot]);
  
  double app=-pref*(2+M*L[iens]*lat[ibeta[iens]][plot_iboot]);
  double full=pref*(-16*C/pow(f0_med,4)*H1_corr[iens]+H2_corr[iens]);
  
  if(PiK==0)
    if(approx_Pi_FSE) return app;
    else return full;
  else
    switch(approx_K_FSE)
      {
      case 0:
	return -FSE_K_corr[iens]; //Silvano put a -
	break;
      case 1:
	return app;
	break;
      case 2:
	return -FSE_K_corr[iens]/(1+FSE_K_corr[iens]/dM2L);
	break;
      case 3:
	return app/(1-app/dM2L);
	break;
      default:
	crash("unknown approx: %d",approx_K_FSE);
	return 0;
	break;
      }
}

void plot_funz_ml(const char *out_path,const char *title,const char *xlab,const char *ylab,bvec &X,bvec Y,bvec Z,bvec &par,double X_phys,boot &chiral_extrap_cont,int plot_PiK)
{
  //setup the plot
  grace out(out_path);
  out.fout<<"@    legend 0.25, 0.8"<<endl;
  out.plot_size(800,600);
  out.plot_title(combine("Chiral extrapolation of %s",title).c_str());
  out.axis_label(xlab,ylab);
  
  //plot the function with error
  int npoint=100;
  double X_pol[npoint];
  bvec Y_pol(npoint,nboot,njack);
  for(int ipoint=0;ipoint<npoint;ipoint++) X_pol[ipoint]=0.0599/(npoint-1)*ipoint+0.0001;
  for(int ib=0;ib<nbeta;ib++)
    {
      bvec Y_pol(npoint,nboot,njack);
      for(int ipoint=0;ipoint<npoint;ipoint++)
	for(int iboot=plot_iboot=0;iboot<nboot+1;plot_iboot=iboot++)
	  {
	    double Pi=fun_fit(par[0][iboot],par[1][iboot],par[2][iboot],par[3][iboot],X_pol[ipoint],X_pol[ipoint],lat[ib][iboot]);
	    double K=fun_fit(par[0][iboot],par[1][iboot],par[2][iboot],par[3][iboot],X_pol[ipoint],ms_phys.med(),lat[ib][iboot]);
	    double eps_gam=K/Pi-1;
	    switch(plot_PiK)
	      {
	      case 0:Y_pol.data[ipoint].data[iboot]=Pi;break;
	      case 1:Y_pol.data[ipoint].data[iboot]=K;break;
	      case 2:Y_pol.data[ipoint].data[iboot]=eps_gam;break;
	      }
	  }
      
      //out.set(1,set_fill_color[ib]);
      //out.polygon(X_pol,Y_pol);
      //out.new_set();
      out.set(1,set_color[ib]);
      out.set_line_size(2);
      out.ave_line(X_pol,Y_pol);
      out.new_set();
    }
  //plot continuum curve
  for(int ipoint=0;ipoint<npoint;ipoint++)
    for(int iboot=plot_iboot=0;iboot<nboot+1;plot_iboot=iboot++)
      {
	double Pi=fun_fit(par[0][iboot],par[1][iboot],par[2][iboot],par[3][iboot],X_pol[ipoint],X_pol[ipoint],0);
	double K=fun_fit(par[0][iboot],par[1][iboot],par[2][iboot],par[3][iboot],X_pol[ipoint],ms_phys.med(),0);
	double eps_gam=K/Pi-1;
	switch(plot_PiK)
	  {
	  case 0:Y_pol.data[ipoint].data[iboot]=Pi;break;
	  case 1:Y_pol.data[ipoint].data[iboot]=K;break;
	  case 2:Y_pol.data[ipoint].data[iboot]=eps_gam;break;
	  }
      }
  out.set(1,"magenta");
  out.set_line_size(3);
  out.ave_line(X_pol,Y_pol);
  out.new_set();
  
  //plot the original data with error  
  for(int ib=0;ib<nbeta;ib++)
    {
      out.set(4,"none",set_symbol[ib],set_color[ib],"filled");
      out.set_legend(set_legend_fm[ib]);
      out.set_line_size(2);
      out.fout<<"@type xydy"<<endl;
      for(int iens=0;iens<nens;iens++) if(ibeta[iens]==ib) out.fout<<X[iens].med()<<" "<<(Y-Z).data[iens]<<endl;
      out.new_set();
    }
  
  //plot the corrected data with error  
  for(int ib=0;ib<nbeta;ib++)
    {
      out.set(4,"none",set_symbol[ib],set_color[ib],"unfilled");
      out.set_line_size(1);
      out.fout<<"@type xydy"<<endl;
      for(int iens=0;iens<nens;iens++) if(ibeta[iens]==ib) out.fout<<X[iens].med()<<" "<<(Y).data[iens]<<endl;
      out.new_set();
    }
  
  //plot the extrapolated point with error
  out.set(4,"none","circle","indigo","filled");
  out.set_line_size(3);
  out.set_legend("Physical point");
  out.print_graph(X_phys,chiral_extrap_cont);
  out.new_set();
}

double *ml_fit;
double *dM2Pi_fit,*err_dM2Pi_fit;
double *dM2K_fit,*err_dM2K_fit;

//calculate the chi square
double chi2(double A,double B,double C,double D,double *a,bool verb=false)
{
  double ch2=0;
  
  for(int iens=0;iens<nens;iens++)
    {
      double teo_dM2Pi=fun_fit(A,B,C,D,ml_fit[iens],ml_fit[iens],a[ibeta[iens]]);
      double teo_dM2K=fun_fit(A,B,C,D,ml_fit[iens],ms,a[ibeta[iens]]);
      
      double uncorr_dM2Pi=dM2Pi_fit[iens];
      double uncorr_dM2K=dM2K_fit[iens];
      
      double correction_dM2Pi=fin_inf_diff(A,iens,MPi[iens][plot_iboot],teo_dM2Pi,0);
      double correction_dM2K=fin_inf_diff(A,iens,MK[iens][plot_iboot],teo_dM2K,1);
      
      double corr_dM2Pi=uncorr_dM2Pi-correction_dM2Pi;
      double corr_dM2K=uncorr_dM2K-correction_dM2K;
      
      double ch2_contr_dM2Pi=pow((corr_dM2Pi-teo_dM2Pi)/err_dM2Pi_fit[iens],2);
      double ch2_contr_dM2K=pow((corr_dM2K-teo_dM2K)/err_dM2K_fit[iens],2);

      ch2+=ch2_contr_dM2Pi+ch2_contr_dM2K;
      
      if(verb) cout<<"Ens: "<<iens<<", contr_Pi: "<<ch2_contr_dM2Pi<<", contr_K: "<<ch2_contr_dM2K<<endl;
    }
  
  return ch2;
}

//wrapper for the calculation of the chi2
void chi2_wr(int &npar,double *fuf,double &ch,double *p,int flag)
{
  double A=p[0],B=p[1],C=p[2],D=p[3];
  double *a=p+4;
  ch=chi2(A,B,C,D,a);
}

void fit(boot &A,boot &B,boot &C,boot &D)
{
  //copy ml
  ml_fit=new double[nens];
  for(int iens=0;iens<nens;iens++) ml_fit[iens]=ml[iens].med();
  
  //alloc dM2Pi and dM2K
  dM2Pi_fit=new double[nens];
  err_dM2Pi_fit=new double[nens];
  dM2K_fit=new double[nens];
  err_dM2K_fit=new double[nens];
  
  TMinuit minu;
  minu.SetPrintLevel(-1);
  
  int npars=4;
  minu.DefineParameter(0,"A",0.0,0.0001,0,0);
  minu.DefineParameter(1,"B",0.0,0.0001,0,0);
  minu.DefineParameter(2,"C",0.0,0.0001,0,0);
  minu.DefineParameter(3,"D",0.0,0.0001,0,0);
  minu.SetFCN(chi2_wr);
  
  double C2;
  for(int iboot=0;iboot<nboot+1;iboot++)
    {
      if(iboot>0) minu.SetPrintLevel(-1);
      
      minu.DefineParameter(4,"a380",lat[0][iboot],0.0001,0,0);
      minu.DefineParameter(5,"a390",lat[1][iboot],0.0001,0,0);
      minu.DefineParameter(6,"a405",lat[2][iboot],0.0001,0,0);
      minu.DefineParameter(7,"a420",lat[3][iboot],0.0001,0,0);
      minu.FixParameter(4);
      minu.FixParameter(5);
      minu.FixParameter(6);
      minu.FixParameter(7);
      
      for(int iens=0;iens<nens;iens++)
        {
          dM2Pi_fit[iens]=dM2Pi.data[iens].data[iboot];
          err_dM2Pi_fit[iens]=dM2Pi.data[iens].err();
          dM2K_fit[iens]=dM2K.data[iens].data[iboot];
          err_dM2K_fit[iens]=dM2K.data[iens].err();
        }
      
      //minimize
      minu.Migrad();
      
      //get back parameters
      double dum;
      minu.GetParameter(0,A.data[iboot],dum);
      minu.GetParameter(1,B.data[iboot],dum);
      minu.GetParameter(2,C.data[iboot],dum);
      minu.GetParameter(3,D.data[iboot],dum);
      
      double lat_med[4]={lat[0].med(),lat[1].med(),lat[2].med(),lat[3].med()};
      if(iboot==nboot) C2=chi2(A.data[iboot],B[iboot],C[iboot],D[iboot],lat_med,true);
    }
  
  //calculate the chi2
  cout<<"A=("<<A<<"), B=("<<B<<"), C=("<<C<<"), D=("<<D<<")"<<endl;
  cout<<"Chi2 = "<<C2<<" / "<<2*nens-npars<<" = "<<C2/(2*nens-npars)<<endl;
  
  delete[] ml_fit;
  delete[] dM2Pi_fit;
  delete[] err_dM2Pi_fit;
  delete[] dM2K_fit;
  delete[] err_dM2K_fit;
}


void plot_funz_a2(const char *out_path,const char *title,const char *xlab,const char *ylab,double *X,bvec &Y,bvec &par,boot &chiral_extrap_cont,int plot_PiK)
{
  //setup the plot
  grace out(out_path);
  out.plot_size(800,600);
  out.plot_title(combine("Continuum extrapolation of %s",title).c_str());
  out.axis_label(xlab,ylab);
  
  //plot the function with error
  double X_pol[100],X2_pol[100];
  bvec Y_pol(100,nboot,njack);
  for(int iens=0;iens<100;iens++)
    {
      X_pol[iens]=0.1/99*iens;
      X2_pol[iens]=X_pol[iens]*X_pol[iens];
	for(int iboot=plot_iboot=0;iboot<nboot+1;plot_iboot=iboot++)
	  {
	    double W;
	    if(plot_PiK) W=ms_phys.med();
	    else W=ml_phys.med();
	    Y_pol.data[iens].data[iboot]=fun_fit(par[0][iboot],par[1][iboot],par[2][iboot],par[3][iboot],ml_phys[iboot],W,X_pol[iens]/hc);
	  }
    }
  out.set(1,"yellow");
  out.polygon(X2_pol,Y_pol);
  out.new_set();
  out.set(1,"red");
  out.set_line_size(2);
  out.ave_line(X2_pol,Y_pol);
  out.new_set();

  //plot the data with error
  for(int ib=0;ib<nbeta;ib++)
    {
      out.set(4,"none",set_symbol[ib],set_color[ib],"filled");
      out.set_legend(set_legend_fm[ib]);
      out.set_line_size(2);
      out.fout<<"@type xydy"<<endl;
      out.fout<<X[ib]*X[ib]<<" "<<Y.data[ib]<<endl;
      out.new_set();
    }
  
  //plot the extrapolated point with error
  out.set(4,"none","circle","indigo","filled");
  out.set_line_size(3);
  out.set_legend("Physical point");
  out.print_graph(0,chiral_extrap_cont);
  out.new_set();
}

void load_iboot(int *iboot_jack,char *ens_name)
{
  char path[1024];
  sprintf(path,"../../../RADIATIVE_CHARMONIUM/DATA1/%s/iboot",ens_name);
  
  FILE *fiboot=fopen(path,"r");
  if(fiboot==NULL)
    {
      perror(combine("Error opening file iboot for ensamble %s",path).c_str());
      exit(1);
    }
  int nr=fread(iboot_jack,sizeof(int),100,fiboot);
  if(nr!=100)
    {
      perror(combine("Error loading iboot data for ensamble %s",path).c_str());
      exit(1);
    }
  fclose(fiboot);
}

int main(int narg,char **arg)
{
  init_latpars();
  
  //read ensemble list, meson masses and meson name
  FILE *an_input_file=open_file("analysis_pars","r");
  read_formatted_from_file_expecting((char*)&approx_Pi_FSE,an_input_file,"%d","approx_Pi_FSE");
  read_formatted_from_file_expecting((char*)&approx_K_FSE,an_input_file,"%d","approx_K_FSE");
  read_formatted_from_file_expecting((char*)&nens,an_input_file,"%d","nens");
  lmass=new double[nens];
  L=new int[nens];
  ibeta=new int[nens];
  FSE_K_corr=new double[nens];
  H1_corr=new double[nens];
  H2_corr=new double[nens];
  dM2Pi=bvec(nens,nboot,njack);
  dM2K=bvec(nens,nboot,njack);
  MPi=bvec(nens,nboot,njack);
  MK=bvec(nens,nboot,njack);
  eps_gam=bvec(nens,nboot,njack);
  
  ofstream bare_data_table("bare_data_table");
  ofstream dim_data_table("dim_data_table");
  
  for(int iens=0;iens<nens;iens++)
    {
      char path[1024];
      read_formatted_from_file((char*)&(ibeta[iens]),an_input_file,"%d","ibeta");
      read_formatted_from_file((char*)&(lmass[iens]),an_input_file,"%lg","lmass");
      read_formatted_from_file((char*)&(L[iens]),an_input_file,"%d","L");
      read_formatted_from_file((char*)&(H1_corr[iens]),an_input_file,"%lg","H1_corr");
      read_formatted_from_file((char*)&(H2_corr[iens]),an_input_file,"%lg","H2_corr");
      read_formatted_from_file((char*)&(FSE_K_corr[iens]),an_input_file,"%lg","FSE_K_corr");
      read_formatted_from_file(path,an_input_file,"%s","path");
      
      //load bare data
      jack aMPi(njack),aMK(njack);
      jack a2dM2Pi(njack),a2dM2K(njack);
      aMPi.load(combine("../%s/MPi",path).c_str());
      a2dM2Pi.load(combine("../%s/DeltaM2Pi",path).c_str());
      aMK.load(combine("../%s/MK",path).c_str());
      a2dM2K.load(combine("../%s/DeltaM2K",path).c_str());
      
      //write the bare data table
      bare_data_table<<iens<<" "<<smart_print(aMPi)<<" "<<smart_print(aMK)<<" "<<smart_print(a2dM2Pi)<<" "<<smart_print(a2dM2K)<<endl;
      
      //load iboot
      int iboot_jack[100];
      load_iboot(iboot_jack,path);
      
      //convert jack from boot
      boot_from_jack(MPi.data[iens],aMPi,iboot_jack);
      boot_from_jack(MK.data[iens],aMK,iboot_jack);
      boot_from_jack(dM2Pi.data[iens],a2dM2Pi,iboot_jack);
      boot_from_jack(dM2K.data[iens],a2dM2K,iboot_jack);
      MPi[iens]/=lat[ibeta[iens]];
      MK[iens]/=lat[ibeta[iens]];
      dM2Pi[iens]/=-sqr(lat[ibeta[iens]]);
      dM2K[iens]/=sqr(lat[ibeta[iens]]);
      eps_gam[iens]=dM2K[iens]/dM2Pi[iens]-1;
      
      //write the dimensional data table
      dim_data_table<<iens<<" "<<smart_print(MPi[iens])<<" "<<smart_print(MK[iens])<<" "
		    <<smart_print(dM2Pi[iens])<<" "<<smart_print(dM2K[iens])<<endl;
    }
  fclose(an_input_file);
  
  //define ml and ref ml
  ml=bvec(nens,nboot,njack);
  for(int iens=0;iens<nens;iens++)
    {
      int b=ibeta[iens],r=ref_ml_beta[b];
      //define ml
      cout<<iens<<" "<<b<<" "<<lmass[iens]<<endl;
      ml[iens]=lmass[iens]/lat[b]/Zp[b];
      //set the lighter mass
      if(r==-1||fabs(ml[r].med()-0.050)>fabs(ml[iens].med()-0.050)) ref_ml_beta[b]=iens;
    }
  cout<<"---"<<endl;
  for(int ib=0;ib<nbeta;ib++)
    if(ref_ml_beta[ib]!=-1)
      cout<<"Ref "<<ib<<" = "<<ref_ml_beta[ib]<<", "<<smart_print(ml[ref_ml_beta[ib]]*1e3)<<" MeV"<<endl;
  cout<<"---"<<endl;

  //perform the fit
  boot A(nboot,njack),B(nboot,njack),C(nboot,njack),D(nboot,njack);
  fit(A,B,C,D);
  
  //chiral extrapolation
  bvec dM2Pi_chir(nbeta,nboot,njack),dM2K_chir(nbeta,nboot,njack);
  boot dM2Pi_chir_cont(nboot,njack),dM2K_chir_cont(nboot,njack);
  bvec dM2Pi_estr_ml(nbeta,nboot,njack),dM2K_estr_ml(nbeta,nboot,njack);
  boot eps_gam_chir_cont(nboot,njack);
  
  //take extrapolation
  for(int iboot=0;iboot<nboot+1;iboot++)
    {
      dM2Pi_chir_cont.data[iboot]=fun_fit(A[iboot],B[iboot],C[iboot],D[iboot],ml_phys[iboot],ml_phys[iboot],0);
      dM2K_chir_cont.data[iboot]=fun_fit(A[iboot],B[iboot],C[iboot],D[iboot],ml_phys[iboot],ms_phys[iboot],0);
      eps_gam_chir_cont.data[iboot]=dM2K_chir_cont.data[iboot]/dM2Pi_chir_cont.data[iboot]-1;
      for(int ib=0;ib<nbeta;ib++)
	{
	  int r=ref_ml_beta[ib];
	  dM2Pi_chir[ib].data[iboot]=fun_fit(A[iboot],B[iboot],C[iboot],D[iboot],ml_phys[iboot],ml_phys[iboot],lat[ib][iboot]);
	  dM2K_chir[ib].data[iboot]=fun_fit(A[iboot],B[iboot],C[iboot],D[iboot],ml_phys[iboot],ms_phys[iboot],lat[ib][iboot]);
	  if(r!=-1)
	    {
	      dM2Pi_estr_ml.data[ib].data[iboot]=dM2Pi[r][iboot]*
		fun_fit(A[nboot],B[nboot],C[nboot],D[nboot],ml_phys[nboot],ml_phys[nboot],0)/
		fun_fit(A[nboot],B[nboot],C[nboot],D[nboot],ml[r][nboot],ml[r][nboot],0);
	      dM2K_estr_ml.data[ib].data[iboot]=dM2Pi[r][iboot]*
		fun_fit(A[nboot],B[nboot],C[nboot],D[nboot],ml_phys[nboot],ms_phys[nboot],0)/
		fun_fit(A[nboot],B[nboot],C[nboot],D[nboot],ml[r][nboot],ms_phys[nboot],0);
	    }
	}
    }
  
  //chiral and continuum
  cout<<"dM2Pi = "<<smart_print(dM2Pi_chir_cont*1e6)<<" MeV^2, physical: 1261 MeV^2"<<endl;
  cout<<"dM2K = "<<smart_print(dM2K_chir_cont*1e6)<<" MeV^2, \"physical(eps=0.7)\": 2000 MeV^2"<<endl;
  cout<<endl;
  
  //plots
  par_res_fit=bvec(4,nboot,njack);
  
  par_res_fit.data[0]=A;
  par_res_fit.data[1]=B;
  par_res_fit.data[2]=C;
  par_res_fit.data[3]=D;
  
  const char tag_ml[1024]="m\\sl\\N\\S\\oMS\\O,2GeV\\N (GeV)";
  const char tag_a2[1024]="a\\S2\\N (fm)";
  double lat_med_fm[4]={lat[0].med()*hc,lat[1].med()*hc,lat[2].med()*hc,lat[3].med()*hc};

  char dM2Pi_label[]={"\\xD\\0M\\S2\\N\\s\\xp\\N\\0 (GeV\\S2\\N)"};
  char dM2K_label[]={"\\xD\\0M\\S2\\N\\sK\\N (GeV\\S2\\N)"};
  char eps_gam_label[]={"\\xe\\0\\s\\xg"};
  
  //compute the corrections
  bvec corr_dM2Pi(nens,nboot,njack),corr_dM2K(nens,nboot,njack);
  for(int iens=0;iens<nens;iens++)
    {
      for(int iboot=0;iboot<nboot+1;iboot++)
	{
	  corr_dM2Pi[iens].data[iboot]=fin_inf_diff(A[iboot],iens,MPi[iens][iboot],fun_fit(A[iboot],B[iboot],C[iboot],D[iboot],ml[iens][iboot],ml[iens][iboot],lat[ibeta[iens]][iboot]),0);
	  corr_dM2K[iens].data[iboot]=fin_inf_diff(A[iboot],iens,MK[iens][iboot],fun_fit(A[iboot],B[iboot],C[iboot],D[iboot],ml[iens][iboot],ms_phys[iboot],lat[ibeta[iens]][iboot]),1);
	}
      
      cout<<iens<<" corr_Pi: "<<smart_print(corr_dM2Pi[iens])<<", corr_K: "<<smart_print(corr_dM2K[iens])<<endl;
    }
  
  plot_funz_ml("dM2Pi_funz_ml.xmg",dM2Pi_label,tag_ml,dM2Pi_label,ml,dM2Pi,corr_dM2Pi,par_res_fit,ml_phys.med(),dM2Pi_chir_cont,0);
  plot_funz_a2("dM2Pi_funz_a2.xmg",dM2Pi_label,tag_a2,dM2Pi_label,lat_med_fm,dM2Pi_estr_ml,par_res_fit,dM2Pi_chir_cont,0);
    
  plot_funz_ml("dM2K_funz_ml.xmg",dM2K_label,tag_ml,dM2K_label,ml,dM2K,corr_dM2K,par_res_fit,ml_phys.med(),dM2K_chir_cont,1);
  plot_funz_a2("dM2K_funz_a2.xmg",dM2K_label,tag_a2,dM2K_label,lat_med_fm,dM2K_estr_ml,par_res_fit,dM2K_chir_cont,1);

  plot_funz_ml("eps_gam.xmg",eps_gam_label,tag_ml,"",ml,eps_gam,corr_dM2Pi*0,par_res_fit,ml_phys.med(),eps_gam_chir_cont,2);

  //write results  
  dM2Pi_chir_cont.write_to_binfile("result_dM2Pi");
  dM2K_chir_cont.write_to_binfile("result_dM2K");

  
  return 0;
}
