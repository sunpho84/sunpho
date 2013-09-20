#include "include.h"
#include "../nf2/common_pars.cpp"

int include_bound=1;
int nth,nens,nth_to_use;
const int nbeta=4,nth_max=8,nens_max=13,njacks=16;
int ibeta[nens_max],use[nens_max];
double lmass[nens_max];
int deg_z=2;
int iboot;
int stranged;

const double MDstar=2.01,MD0star=2.32;
const double MDSstar=2.112,MD0Sstar=2.32;
const double MD_ph=1.868,MK_ph=0.4944,MP_ph=0.135;
const double Q2_max_ph[2]={sqr(MD_ph-MP_ph),sqr(MD_ph-MK_ph)};
bvec Q2,FP,F0,FT,MD,MP,EP,Z;

bvec ml;
int ref_ml_beta[4]={-1,-1,-1,-1};

const char set_color[nbeta][1024]={"black","blue","red","green4"};
const char set_fill_color[nbeta][1024]={"grey","turquoise","yellow","green"};
const char set_symbol[nbeta][1024]={"square","circle","triangle up","triangle left"};
const char set_legend[nbeta][1024]={"\\xb\\0=3.80","\\xb\\0=3.90","\\xb\\0=4.05","\\xb\\0=4.20"};
const char set_legend_fm[nbeta][1024]={"a = 0.098 fm","a = 0.085 fm","a = 0.067 fm","a = 0.054 fm"};

int study_T_fr_P=1;

template <class T1,class T2> T1 fun_Z(T1 MP,T1 MD,T2 Q2)
{
  double T0=0;//(MD_ph+MP_ph)*sqr(sqrt(MD_ph)-sqrt(MP_ph));
  T1 TP=sqr(MD+MP);
  T1 Z=sqrt(TP-Q2)-sqrt(TP-T0);
  Z/=sqrt(TP-Q2)+sqrt(TP-T0);
  
  return Z;
}

template <class T1,class T2> T1 fun_Zone(T1 MP,T1 MD,T2 Z)
{
  double T0=0;//(MD_ph+MP_ph)*sqr(sqrt(MD_ph)-sqrt(MP_ph));
  T1 TP=sqr(MD+MP),B=sqrt(TP-T0),A=B*(1+Z)/(1-Z);
  
  return TP-sqr(A);
}

double fun_fit_par(double A,double B,double C,double ml,double a)
{return A + B*ml + C*a*a;}

void plot_funz_ml(const char *out_path,const char *title,const char *xlab,const char *ylab,bvec &X,bvec &Y,bvec &par,double X_phys,double (*fun)(double,double,double,double,double),boot &chiral_extrap_cont)
{
  //setup the plot
  grace out(out_path);
  out.plot_size(800,600);
  out.plot_title(combine("Chiral extrapolation of %s",title).c_str());
  out.axis_label(xlab,ylab);
  
  if(fun!=NULL)
    {
      //plot the function with error
      int npoint=100;
      double X_pol[npoint];
      bvec Y_pol(npoint,nboot,njack);
      for(int ipoint=0;ipoint<npoint;ipoint++) X_pol[ipoint]=0.0599/(npoint-1)*ipoint+0.0001;
      for(int ib=0;ib<nbeta;ib++)
        {
          bvec Y_pol(npoint,nboot,njack);
          for(int ipoint=0;ipoint<npoint;ipoint++)
            for(iboot=0;iboot<=nboot;iboot++)
              Y_pol.data[ipoint].data[iboot]=fun(par[0][iboot],par[1][iboot],par[2][iboot],X_pol[ipoint],lat[ib][iboot]);
          
          out.set(1,set_fill_color[ib]);
          out.polygon(X_pol,Y_pol);
          out.new_set();
          out.set(1,set_color[ib]);
          out.set_line_size(2);
          out.ave_line(X_pol,Y_pol);
          out.new_set();
        }
      //plot continuum curve
      for(int ipoint=0;ipoint<npoint;ipoint++)
        for(iboot=0;iboot<=nboot;iboot++)
          Y_pol.data[ipoint]=fun(par[0][iboot],par[1][iboot],par[2][iboot],X_pol[ipoint],0);
      out.set(1,"magenta");
      out.set_line_size(3);
      out.ave_line(X_pol,Y_pol);
      out.new_set();
    }
  
  //plot the original data with error  
  for(int ib=0;ib<nbeta;ib++)
    {
      out.set(4,"none",set_symbol[ib],set_color[ib],"filled");
      out.set_legend(set_legend_fm[ib]);
      out.set_line_size(2);
      out.fout<<"@type xydy"<<endl;
      for(int iens=0;iens<nens;iens++)
	if(use[iens])
	  if(ibeta[iens]==ib)
	    out.fout<<X[iens].med()<<" "<<Y.data[iens]<<endl;
      out.new_set();
    }
  
  //plot the extrapolated point with error
  out.set(4,"none","circle","indigo","filled");
  out.set_line_size(3);
  out.set_legend("Physical point");
  out.print_graph(X_phys,chiral_extrap_cont);
  out.new_set();
}

void plot_funz_a2(const char *out_path,const char *title,const char *xlab,const char *ylab,double *X,bvec &Y,bvec &par,double (*fun)(double,double,double,double,double),boot &chiral_extrap_cont)
{
  //setup the plot
  grace out(out_path);
  out.plot_size(800,600);
  out.plot_title(combine("Continuum extrapolation of %s",title).c_str());
  out.axis_label(xlab,ylab);
  
  //plot the function with error
  int npoints_plot=100;
  double X_pol[npoints_plot],X2_pol[npoints_plot];
  bvec Y_pol(npoints_plot,nboot,njack);
  for(int ipoint=0;ipoint<npoints_plot;ipoint++)
    {
      X_pol[ipoint]=0.1/(npoints_plot-1)*ipoint;
      X2_pol[ipoint]=X_pol[ipoint]*X_pol[ipoint];
      for(iboot=0;iboot<=nboot;iboot++)
	Y_pol[ipoint].data[iboot]=fun(par[0][iboot],par[1][iboot],par[2][iboot],ml_phys[iboot],X_pol[ipoint]/hc);
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
      out.fout<<"@type xydy"<<endl;
      out.set(4,"none",set_symbol[ib],set_color[ib],"filled");
      out.set_legend(set_legend_fm[ib]);
      out.set_line_size(2);
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

double *X_fit;
double *Y_fit,*err_Y_fit;

//calculate the chi square
void chi2_par(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  for(int iens=0;iens<nens;iens++)
    if(use[iens])
      {
	double Y_teo=fun_fit_par(p[0],p[1],p[2],ml[iens][iboot],lat[ibeta[iens]][iboot]);
	double contr=pow((Y_fit[iens]-Y_teo)/err_Y_fit[iens],2);
	ch+=contr;
	if(flag==1) cout<<"contr iens "<<iens<<", ml="<<ml[iens][iboot]<<", a="<<lat[ibeta[iens]][iboot]<<", ("<<Y_fit[iens]<<"-"<<Y_teo<<")/"<<err_Y_fit[iens]<<"="<<contr<<endl;
      }
}

boot fit_par(bvec &X,bvec &Y,const char *title,int include_ml_term=1,int include_a2_term=1)
{
  //copy X
  X_fit=new double[nens];
  for(int iens=0;iens<nens;iens++) if(use[iens]) X_fit[iens]=X[iens].med();
  Y_fit=new double[nens];
  err_Y_fit=new double[nens];
  
  TMinuit minu;
  
  int npars=3;
  minu.DefineParameter(0,"A",1.0,0.1,0,0);
  minu.DefineParameter(1,"B",0,0.1,0,0);
  minu.DefineParameter(2,"C",0.0,0.1,0,0);
  if(!include_ml_term) minu.FixParameter(1);
  if(!include_a2_term) minu.FixParameter(2);
  minu.SetFCN(chi2_par);
  
  boot C2(nboot,njack);
  boot A(nboot,njack),B(nboot,njack),C(nboot,njack);
  for(iboot=0;iboot<=nboot;iboot++)
    {
      if(iboot!=nboot) minu.SetPrintLevel(-1);
      else minu.SetPrintLevel(1);
      
      for(int iens=0;iens<nens;iens++)
	if(use[iens])
	  {
	    Y_fit[iens]=Y.data[iens].data[iboot];
	    err_Y_fit[iens]=Y.data[iens].err();
	  }
      
      //minimize
      minu.Migrad();
      minu.mnimpr();
      minu.Migrad();
      
      //get back parameters
      double dum;
      minu.GetParameter(0,A.data[iboot],dum);
      minu.GetParameter(1,B.data[iboot],dum);
      minu.GetParameter(2,C.data[iboot],dum);
      
      double par_boot[3]={A[iboot],B[iboot],C[iboot]};
      chi2_par(npars,NULL,C2.data[iboot],par_boot,iboot==nboot);
    }
  
  //calculate the chi2
  int ndofs=nens-(1+include_a2_term+include_ml_term);
  cout<<"A=("<<A<<"), B=("<<B<<"), C=("<<C<<")"<<endl;
  cout<<"Chi2 = "<<C2<<" / "<<ndofs<<" = "<<C2/ndofs<<endl;
  
  delete[] X_fit;
  delete[] Y_fit;
  delete[] err_Y_fit;
  
  //plots

  //chiral extrapolation
  bvec Y_chir(nbeta,nboot,njack);
  boot Y_chir_cont(nboot,njack);
  bvec Y_estr_ml(nbeta,nboot,njack);
  
  for(iboot=0;iboot<=nboot;iboot++)
    {
      Y_chir_cont.data[iboot]=fun_fit_par(A[iboot],B[iboot],C[iboot],ml_phys[iboot],0);
      for(int ib=0;ib<nbeta;ib++)
	{
	  int r=ref_ml_beta[ib];
	  Y_chir[ib].data[iboot]=fun_fit_par(A[iboot],B[iboot],C[iboot],ml_phys[iboot],lat[ib][iboot]);
	  if(r!=-1)
	    Y_estr_ml.data[ib].data[iboot]=Y[r][iboot]*
	      fun_fit_par(A[nboot],B[nboot],C[nboot],ml_phys[nboot],lat[ib][iboot])/
	      fun_fit_par(A[nboot],B[nboot],C[nboot],ml[r][nboot],lat[ib][iboot]);
	}
    }
  
  //chiral and continuum
  cout<<title<<" = "<<Y_chir_cont<<endl;
  cout<<endl;
  
  bvec par_res_fit_F(3,nboot,njack);
  
  par_res_fit_F[0]=A;
  par_res_fit_F[1]=B;
  par_res_fit_F[2]=C;
  
  const char tag_ml[1024]="m\\sl\\N\\SMS,2GeV\\N (GeV)";
  const char tag_a2[1024]="a\\S2\\N (fm)";
  double lat_med_fm[4]={lat[0].med()*hc,lat[1].med()*hc,lat[2].med()*hc,lat[3].med()*hc};
  
  plot_funz_ml(combine("%s_funz_ml.xmg",title).c_str(),title,tag_ml,title,ml,Y,par_res_fit_F,
	       ml_phys.med(),fun_fit_par,Y_chir_cont);
  plot_funz_a2(combine("%s_funz_a2.xmg",title).c_str(),title,tag_a2,title,lat_med_fm,
	       Y_estr_ml,par_res_fit_F,fun_fit_par,Y_chir_cont);
  
  return Y_chir_cont;
}

int icombo(int iens,int ith)
{return iens*nth_to_use+ith;}

int ipar_combo(int d,int iff)
{
  if(d<0||d>=deg_z) crash("d");
  if(iff<0||iff>=3) crash("iff");
  
  if(d==0&&iff==0) return deg_z*1;
  else return deg_z*iff+d;
}

template <class T1,class T2> void put_or_remove_pole(T1 &fP,T1 &f0,T1 &fT,T2 q2,int put_remove)
{
  T2 factP,fact0,factT;
  if(!stranged)
    {
      factP=(1-q2/sqr(MDstar));
      factT=(1-q2/sqr(MDstar));
      fact0=(1-q2/sqr(MD0star));
    }
  else
    {
      factP=(1-q2/sqr(MDSstar));
      factT=(1-q2/sqr(MDSstar));
      fact0=(1-q2/sqr(MD0Sstar));
    }
  if(put_remove==1)
    {
      fP*=factP;
      fT*=factT;
      f0*=fact0;
    }
  else
    {
      fP/=factP;
      fT/=factT;
      f0/=fact0;
    }
}

template <class T1,class T2> void put_pole(T1 &fP,T1 &f0,T1 &fT,T2 q2)
{put_or_remove_pole(fP,f0,fT,q2,0);}
template <class T1,class T2> void remove_pole(T1 &fP,T1 &f0,T1 &fT,T2 q2)
{put_or_remove_pole(fP,f0,fT,q2,1);}

void boot_from_jack(boot &out,jack in,int *iboot_jack)
{
  int nboot=out.nboot;
  int njack=in.njack;
  for(iboot=0;iboot<nboot;iboot++) out.data[iboot]=in.data[iboot_jack[iboot]];
  out.data[nboot]=in.data[njack];
}

void load_iboot(int *iboot_jack,char *ens_name)
{
  char path[1024];
  sprintf(path,"/Users/francesco/QCD/LAVORI/RADIATIVE_CHARMONIUM/DATA1/%s/iboot",ens_name);
  
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

boot fun_plot_Z_chir_cont(bvec p,double x)
{
  boot out=p[0];
  for(int d=1;d<deg_z;d++) out+=pow(x,d)*p[d];
  if(include_bound) out+=pow(x,deg_z)*p[deg_z-1]/(deg_z);
  return out;
}

boot fun_plot_Z(bvec *p,double x,int iens)
{
  boot out=p[0][iens];
  for(int d=1;d<deg_z;d++) out+=pow(x,d)*p[d][iens];
  if(include_bound) out+=pow(x,deg_z)*p[deg_z-1][iens]/(deg_z);
  return out;
}

template <class T> T fun_fit_Z(T *p,T x,int iff)
{
  T out=p[ipar_combo(0,iff)];
  for(int d=1;d<deg_z;d++) out+=pow(x,d)*p[ipar_combo(d,iff)];
  if(include_bound) out+=pow(x,deg_z)*p[ipar_combo(deg_z-1,iff)]/(deg_z);
  return out;
}

double *x_fit_Z;
double *y_num_fit_Z[3];
double *y_err_fit_Z[3];
void chi2_z(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  for(int ith=0;ith<nth_to_use;ith++)
    for(int iff=0;iff<3;iff++)
      if(ith||iff==1)
	{
	  double y_teo=fun_fit_Z(p,x_fit_Z[ith],iff);
	  double contr=pow((y_num_fit_Z[iff][ith]-y_teo)/y_err_fit_Z[iff][ith],2);
	  ch+=contr;
	  if(flag==1) cout<<"contr ith "<<ith<<", iff="<<iff<<", z="<<x_fit_Z[ith]<<" ("<<
			y_num_fit_Z[iff][ith]<<"-"<<y_teo<<")/"<<y_err_fit_Z[iff][ith]<<"="<<contr<<endl;
	}
}

void fit_Z(bvec *pars_P,bvec *pars_0,bvec *pars_T)
{
  TMinuit minu;
  minu.SetFCN(chi2_z);
  int npars=deg_z*3;
  //set pars
  for(int ipar=0;ipar<npars;ipar++) minu.DefineParameter(ipar,combine("P%d",ipar).c_str(),1,0.0001,0,0);
  minu.FixParameter(0);
  
  x_fit_Z=new double[nth_to_use];
  for(int iff=0;iff<3;iff++)
    {
      y_err_fit_Z[iff]=new double[nth_to_use];
      y_num_fit_Z[iff]=new double[nth_to_use];
    }
  
  for(int iens=0;iens<nens;iens++)
    if(use[iens])
      {
	for(int ith=0;ith<nth_to_use;ith++)
	  {
	    int ic=icombo(iens,ith);
	    y_err_fit_Z[0][ith]=FP[ic].err();
	    y_err_fit_Z[1][ith]=F0[ic].err();
	    y_err_fit_Z[2][ith]=FT[ic].err();
	  }
	
	boot C2(nboot,njack);
	for(iboot=0;iboot<=nboot;iboot++)
	  {
	    for(int ith=0;ith<nth_to_use;ith++)
	      {
		int ic=icombo(iens,ith);
		y_num_fit_Z[0][ith]=FP[ic][iboot];
		y_num_fit_Z[1][ith]=F0[ic][iboot];
		y_num_fit_Z[2][ith]=FT[ic][iboot];
		x_fit_Z[ith]=Z[ic][iboot];
	      }
	
	    //quiet if not central boot
	    if(iboot!=nboot) minu.SetPrintLevel(-1);
	    else minu.SetPrintLevel(1);
	    //minu.Migrad();
	    //minu.mnimpr();
	    minu.Migrad();
      
	    //get back parameters
	    double dum;
	    double par_boot[3*deg_z];
	    for(int d=0;d<deg_z;d++)
	      {
		minu.GetParameter(ipar_combo(d,0),pars_P[d][iens].data[iboot],dum);
		minu.GetParameter(ipar_combo(d,1),pars_0[d][iens].data[iboot],dum);
		minu.GetParameter(ipar_combo(d,2),pars_T[d][iens].data[iboot],dum);
		par_boot[0*deg_z+d]=pars_P[d][iens][iboot];
		par_boot[1*deg_z+d]=pars_0[d][iens][iboot];
		par_boot[2*deg_z+d]=pars_T[d][iens][iboot];
	      }
	  
	    chi2_z(npars,NULL,C2.data[iboot],par_boot,iboot==nboot);
	  }
	
	int ndofs=nth_to_use*3-2-5;
	cout<<"Chi2 = "<<C2<<" / "<<ndofs<<" = "<<C2/ndofs<<endl;
	
	//max Z
	double Z_min=min(Z[icombo(iens,0)].med(),Z[icombo(iens,nth_to_use-1)].med());
	double Z_max=0;//max(Z[icombo(iens,0)].med(),Z[icombo(iens,nth_to_use-1)].med());
	
	//plot
	int npoints_plot=100;
	double x_plot[npoints_plot],dx=(Z_max-Z_min)/(npoints_plot-1);
	bvec y_plot_P(npoints_plot,nboot,njack),y_plot_0(npoints_plot,nboot,njack),y_plot_T(npoints_plot,nboot,njack);
	for(int ipoint_plot=0;ipoint_plot<npoints_plot;ipoint_plot++)
	  {
	    x_plot[ipoint_plot]=Z_min+dx*ipoint_plot;
	    y_plot_P[ipoint_plot]=fun_plot_Z(pars_P,x_plot[ipoint_plot],iens);
	    y_plot_0[ipoint_plot]=fun_plot_Z(pars_0,x_plot[ipoint_plot],iens);
	    y_plot_T[ipoint_plot]=fun_plot_Z(pars_T,x_plot[ipoint_plot],iens);
	  }
	
	{
	  ofstream outP(combine("fit_in_Z/plots/FP_ens_%d.xmg",iens).c_str());
	  ofstream out0(combine("fit_in_Z/plots/F0_ens_%d.xmg",iens).c_str());
	  ofstream outT(combine("fit_in_Z/plots/FT_ens_%d.xmg",iens).c_str());
	  outP<<"@type xydy"<<endl;
	  out0<<"@type xydy"<<endl;
	  outT<<"@type xydy"<<endl;
	  outP<<write_polygon(x_plot,y_plot_P);
	  out0<<write_polygon(x_plot,y_plot_0);
	  outT<<write_polygon(x_plot,y_plot_T);
	  outP<<"&"<<endl;
	  out0<<"&"<<endl;
	  outT<<"&"<<endl;
	  outP<<write_ave_line(x_plot,y_plot_P);
	  out0<<write_ave_line(x_plot,y_plot_0);
	  outT<<write_ave_line(x_plot,y_plot_T);
	  outP<<"&"<<endl;
	  out0<<"&"<<endl;
	  outT<<"&"<<endl;
	  outP<<"@type xydy"<<endl;
	  out0<<"@type xydy"<<endl;
	  outT<<"@type xydy"<<endl;
	  for(int ith=0;ith<nth_to_use;ith++)
	    {
	      if(ith) outP<<Z[icombo(iens,ith)].med()<<" "<<FP[icombo(iens,ith)]<<endl;
	      out0<<Z[icombo(iens,ith)].med()<<" "<<F0[icombo(iens,ith)]<<endl;
	      if(ith) outT<<Z[icombo(iens,ith)].med()<<" "<<FT[icombo(iens,ith)]<<endl;
	    }
	}
      }
  
}

int main(int narg,char **arg)
{
  init_latpars();
  
  //read ensemble list, meson masses and meson name
  FILE *an_input_file=open_file("analysis_pars","r");
  char base_path_ff[200],base_path_MD[200],base_path_MP[200];
  
  read_formatted_from_file_expecting(base_path_ff,an_input_file,"%s","base_path_ff");
  read_formatted_from_file_expecting(base_path_MD,an_input_file,"%s","base_path_MD");
  read_formatted_from_file_expecting(base_path_MP,an_input_file,"%s","base_path_MP");
  
  int use_also_divita;
  read_formatted_from_file_expecting((char*)(&use_also_divita),an_input_file,"%d","use_also_divita");
  if(use_also_divita) crash("divita is not our friend");
  read_formatted_from_file_expecting((char*)(&stranged),an_input_file,"%d","stranged");

  read_formatted_from_file_expecting((char*)(&nth),an_input_file,"%d","nth");
  read_formatted_from_file_expecting((char*)(&nth_to_use),an_input_file,"%d","nth_to_use");
  read_formatted_from_file_expecting((char*)(&nens),an_input_file,"%d","nens");
  if(nth>nth_max) crash("nth>nth_max");
  if(nth_to_use>nth) crash("nth_to_use>nth");
  if(nens>nens_max) crash("nens>nens_max");
  
  //read data
  ofstream bare_data_table("bare_data_table");
  jvec f_q2(nens*nth_to_use,njacks);
  char ens_path[nens][100];
  FP=F0=FT=EP=Q2=Z=bvec(nens*nth_to_use,nboot,njack);
  MD=MP=bvec(nens,nboot,njacks);
  
  for(int iens=0;iens<nens;iens++)
    {
      read_formatted_from_file((char*)&(use[iens]),an_input_file,"%d","use");
      read_formatted_from_file((char*)&(ibeta[iens]),an_input_file,"%d","ibeta");
      read_formatted_from_file((char*)&(lmass[iens]),an_input_file,"%lg","lmass");
      read_formatted_from_file(ens_path[iens],an_input_file,"%s","ens_path");
      
      //load iboot
      int iboot_jack[100];
      load_iboot(iboot_jack,ens_path[iens]);
      
      //read data
      jvec EP_j(nth,njack),Q2_j(nth,njack),FP_j(nth,njack),F0_j(nth,njack),FT_j(nth,njack),F0S_j(nth,njack);
      EP_j.load(combine(base_path_ff,ens_path[iens]).c_str(),0);
      Q2_j.load(combine(base_path_ff,ens_path[iens]).c_str(),1);
      FP_j.load(combine(base_path_ff,ens_path[iens]).c_str(),2);
      F0_j.load(combine(base_path_ff,ens_path[iens]).c_str(),3);
      FT_j.load(combine(base_path_ff,ens_path[iens]).c_str(),4);
      F0S_j.load(combine(base_path_ff,ens_path[iens]).c_str(),5);
      
      //read the corr
      jvec corr(5,njack);
      int ff_tag[2]={0,1};
      corr.load(combine("/Users/francesco/QCD/LAVORI/B3PTS_DAMIR/%s/CORRECTIVE/corrective_factor",
			ens_path[iens]).c_str(),ff_tag[stranged]);
      
      for(int ith=0;ith<5;ith++)
	{
	  FP_j[ith]/=corr[ith].med();
	  F0_j[ith]/=corr[ith].med();
	  FT_j[ith]/=corr[ith].med();
	  F0S_j[ith]/=corr[ith].med();
	}
	  
      //convert jack to boot and pass to dimensionful quantity
      jack MD_j(njacks),MP_j(njacks);
      MD_j.load(combine(base_path_MD,ens_path[iens]).c_str(),0);
      MP_j.load(combine(base_path_MP,ens_path[iens]).c_str(),0);
      boot_from_jack(MD[iens],MD_j,iboot_jack);
      boot_from_jack(MP[iens],MP_j,iboot_jack);
      MD[iens]/=lat[ibeta[iens]];
      MP[iens]/=lat[ibeta[iens]];
      
      //write the bare data table while converting jack to boot
      for(int ith=0;ith<nth_to_use;ith++)
        {
          int ic=icombo(iens,ith),ib=ibeta[iens];
          boot_from_jack(EP.data[ic],EP_j[ith],iboot_jack);
          boot_from_jack(Q2.data[ic],Q2_j[ith],iboot_jack);
          boot_from_jack(FP.data[ic],FP_j[ith],iboot_jack);
          boot_from_jack(F0.data[ic],F0_j[ith],iboot_jack);
          boot_from_jack(FT.data[ic],FT_j[ith],iboot_jack);
          
          //pass to dimensionful quantity
          Q2[ic]/=sqr(lat[ib]);          
          EP[ic]/=lat[ib];
	  
	  //compute
	  Z[ic]=fun_Z(MP[iens],MD[iens],Q2[ic]);
	  
	  //remove pole
	  remove_pole(FP[ic],F0[ic],FT[ic],Q2[ic]);
	  
	  //divide FT by FP
	  if(study_T_fr_P) FT[ic]/=FP[ic];
	}
    }
  fclose(an_input_file);
  
  //define ml and ref ml
  ml=bvec(nens,nboot,njack);
  for(int iens=0;iens<nens;iens++)
    if(use[iens])
      {      
	int b=ibeta[iens],r=ref_ml_beta[b];
	//define ml
	cout<<iens<<" "<<b<<" "<<lmass[iens]<<endl;
	ml[iens]=lmass[iens]/lat[b]/Zp[b];
	//set the lighter mass
	if(r==-1||fabs(ml[r].med()-0.050)>fabs(ml[iens].med()-0.050)) ref_ml_beta[b]=iens;
      }
  cout<<"---"<<endl;
  for(int ib=0;ib<nbeta;ib++) if(ref_ml_beta[ib]!=-1) cout<<"Ref "<<ib<<" = "<<ref_ml_beta[ib]<<", "<<ml[ref_ml_beta[ib]]<<" MeV"<<endl;
  cout<<"---"<<endl;
  
  //fit the ff as function of z
  bvec pars_P[deg_z],pars_0[deg_z],pars_T[deg_z];
  for(int d=0;d<deg_z;d++) pars_P[d]=pars_0[d]=pars_T[d]=bvec(nens,nboot,njack);
  fit_Z(pars_P,pars_0,pars_T);
  
  //estrapolate pars to the continuum
  bvec pars_P_chir_cont(deg_z,nboot,njack),pars_0_chir_cont(deg_z,nboot,njack),pars_T_chir_cont(deg_z,nboot,njack);
  for(int d=0;d<deg_z;d++)
    {
      int include_ml_term=1,include_a2_term=1;
      //if(d>0) include_ml_term=0;
      pars_P_chir_cont[d]=fit_par(ml,pars_P[d],combine("fit_in_Z/plots/pars_P_deg_%d",d).c_str(),
				  include_ml_term,include_a2_term);
      pars_0_chir_cont[d]=fit_par(ml,pars_0[d],combine("fit_in_Z/plots/pars_0_deg_%d",d).c_str(),
				  include_ml_term,include_a2_term);
      pars_T_chir_cont[d]=fit_par(ml,pars_T[d],combine("fit_in_Z/plots/pars_T_deg_%d",d).c_str(),
				  include_ml_term,include_a2_term);
    }
    
  //max Z
  double Z_min=fun_Z(stranged?MK_ph:MP_ph,MD_ph,Q2_max_ph[stranged]);
  double Z_max=0;
  
  if(include_bound && deg_z>2) crash(""); 
  
  //plot
  int npoints_plot=100;
  double Z_plot[npoints_plot],dZ=(Z_max-Z_min)/(npoints_plot-1);
  double Q2_plot[npoints_plot];
  bvec y_plot_P(npoints_plot,nboot,njack),y_plot_0(npoints_plot,nboot,njack),y_plot_T(npoints_plot,nboot,njack);
  for(int ipoint_plot=0;ipoint_plot<npoints_plot;ipoint_plot++)
    {
      Z_plot[ipoint_plot]=Z_min+dZ*ipoint_plot;
      y_plot_P[ipoint_plot]=fun_plot_Z_chir_cont(pars_P_chir_cont,Z_plot[ipoint_plot]);
      y_plot_0[ipoint_plot]=fun_plot_Z_chir_cont(pars_0_chir_cont,Z_plot[ipoint_plot]);
      y_plot_T[ipoint_plot]=fun_plot_Z_chir_cont(pars_T_chir_cont,Z_plot[ipoint_plot]);
      
      Q2_plot[ipoint_plot]=fun_Zone(stranged?MK_ph:MP_ph,MD_ph,Z_plot[ipoint_plot]);
      
      //return to T
      if(study_T_fr_P) y_plot_T[ipoint_plot]*=y_plot_P[ipoint_plot];

      put_pole(y_plot_P[ipoint_plot],y_plot_0[ipoint_plot],y_plot_T[ipoint_plot],Q2_plot[ipoint_plot]);      
      
      //exp
      if(study_T_fr_P) y_plot_T[ipoint_plot]/=y_plot_P[ipoint_plot];
    }
  
  for(int Z_Q2_flag=0;Z_Q2_flag<2;Z_Q2_flag++)
    {
      char Z_Q2_tag[2][4]={"Z","Q2"};
      ofstream outP(combine("fit_in_Z/plots/FP_chir_cont_funz_%s.xmg",Z_Q2_tag[Z_Q2_flag]).c_str());
      ofstream out0(combine("fit_in_Z/plots/F0_chir_cont_funz_%s.xmg",Z_Q2_tag[Z_Q2_flag]).c_str());
      ofstream outT(combine("fit_in_Z/plots/FT_chir_cont_funz_%s.xmg",Z_Q2_tag[Z_Q2_flag]).c_str());
      
      double *K_plot=Z_Q2_flag?Q2_plot:Z_plot;
      
      outP<<write_polygon(K_plot,y_plot_P);
      out0<<write_polygon(K_plot,y_plot_0);
      outT<<write_polygon(K_plot,y_plot_T);
      outP<<"&"<<endl;
      out0<<"&"<<endl;
      outT<<"&"<<endl;
      outP<<write_ave_line(K_plot,y_plot_P);
      out0<<write_ave_line(K_plot,y_plot_0);
      outT<<write_ave_line(K_plot,y_plot_T);
      outP<<"&"<<endl;
      out0<<"&"<<endl;
      outT<<"&"<<endl;
    }
  
  return 0;  
}  
