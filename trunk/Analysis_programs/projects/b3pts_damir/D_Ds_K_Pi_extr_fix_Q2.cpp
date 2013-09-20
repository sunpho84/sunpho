#include "include.h"
#include "../nf2/common_pars.cpp"

int nth,nens,nth_to_use;
const int nbeta=4,nth_max=8,nens_max=13,njacks=16;
int ibeta[nens_max],use[nens_max];
double lmass[nens_max];

int plot_ijack;
int stranged;
int include_a2_term=1;
int include_ml_term=1;

double ext_ch2;
const double MDstar=2.01,MD0star=2.32;
const double MDSstar=2.112,MD0Sstar=2.32;
const double MD_ph=1.868,MK_ph=0.4944,MP_ph=0.135;
const double Q2_max_ph[2]={sqr(MD_ph-MP_ph),sqr(MD_ph-MK_ph)};
bvec Q2,FP,F0,FT,MD,MP,EP;

bvec ml;
int ref_ml_beta[4]={-1,-1,-1,-1};

bvec par_res_fit_F;

const char set_color[nbeta][1024]={"black","blue","red","green4"};
const char set_fill_color[nbeta][1024]={"grey","turquoise","yellow","green"};
const char set_symbol[nbeta][1024]={"square","circle","triangle up","triangle left"};
const char set_legend[nbeta][1024]={"\\xb\\0=3.80","\\xb\\0=3.90","\\xb\\0=4.05","\\xb\\0=4.20"};
const char set_legend_fm[nbeta][1024]={"a = 0.098 fm","a = 0.085 fm","a = 0.067 fm","a = 0.054 fm"};

int study_T_fr_P=1;
int reciproc_data=1;
int plot_iboot;

double fun_fit_F(double A,double B,double C,double ml,double a)
{return A*(1 + B*ml + C*a*a);}

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
            for(int iboot=plot_iboot=0;iboot<=nboot;plot_iboot=iboot++)
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
        for(int iboot=plot_iboot=0;iboot<nboot+1;plot_iboot=iboot++)
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
  double X_pol[100],X2_pol[100];
  bvec Y_pol(100,nboot,njack);
  for(int iens=0;iens<100;iens++)
    if(use[iens])
      {
	X_pol[iens]=0.1/99*iens;
	X2_pol[iens]=X_pol[iens]*X_pol[iens];
	for(int iboot=plot_iboot=0;iboot<nboot+1;plot_iboot=iboot++)
	  Y_pol.data[iens].data[iboot]=fun(par[0][iboot],par[1][iboot],par[2][iboot],ml_phys[iboot],X_pol[iens]/hc);
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
void chi2_fun(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  for(int iens=0;iens<nens;iens++)
    if(use[iens])
      {
        double Y_teo=fun_fit_F(p[0],p[1],p[2],X_fit[iens],lat[ibeta[iens]][plot_iboot]);
        double contr=pow((Y_fit[iens]-Y_teo)/err_Y_fit[iens],2);
        ch+=contr;
      }
}

boot fit(bvec &X,bvec &Y,const char *title)
{
  //copy X
  X_fit=new double[nens];
  for(int iens=0;iens<nens;iens++) if(use[iens]) X_fit[iens]=X[iens].med();
  Y_fit=new double[nens];
  err_Y_fit=new double[nens];
  
  TMinuit minu;
  
  int npars=3;
  minu.DefineParameter(0,"A",1.0,0.0001,0,0);
  minu.DefineParameter(1,"B",0.0,0.0001,0,0);
  minu.DefineParameter(2,"C",0.0,0.0001,0,0);
  if(!include_ml_term) minu.FixParameter(1);
  if(!include_a2_term) minu.FixParameter(2);
  minu.SetFCN(chi2_fun);
  
  boot C2(nboot,njack);
  boot A(nboot,njack),B(nboot,njack),C(nboot,njack);
  for(int iboot=0;iboot<nboot+1;iboot++)
    {
      plot_iboot=iboot;
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
      chi2_fun(npars,NULL,C2.data[iboot],par_boot,0);
    }
  
  int ninc_ens=0;
  for(int iens=0;iens<nens;iens++)
    if(use[iens]) ninc_ens++;
  
  //calculate the chi2
  cout<<"A=("<<A<<"), B=("<<B<<"), C=("<<C<<")"<<endl;
  cout<<"Chi2 = "<<C2<<" / "<<ninc_ens-npars<<" = "<<C2/(ninc_ens-npars)<<endl;
  
  delete[] X_fit;
  delete[] Y_fit;
  delete[] err_Y_fit;
  
  //plots

  //chiral extrapolation
  bvec Y_chir(nbeta,nboot,njack);
  boot Y_chir_cont(nboot,njack);
  bvec Y_estr_ml(nbeta,nboot,njack);
  for(int iboot=0;iboot<nboot+1;iboot++)
    {
      Y_chir_cont.data[iboot]=fun_fit_F(A[iboot],B[iboot],C[iboot],ml_phys[iboot],0);
      for(int ib=0;ib<nbeta;ib++)
	{
	  int r=ref_ml_beta[ib];
	  Y_chir[ib].data[iboot]=fun_fit_F(A[iboot],B[iboot],C[iboot],ml_phys[iboot],lat[ib][iboot]);
	  if(r!=-1)
	    Y_estr_ml.data[ib].data[iboot]=Y[r][iboot]*fun_fit_F(A[nboot],B[nboot],C[nboot],ml_phys[nboot],0)/fun_fit_F(A[nboot],B[nboot],C[nboot],ml[r][nboot],0);
	}
    }
  
  //chiral and continuum
  cout<<title<<" = "<<Y_chir_cont<<endl;
  cout<<endl;
  
  par_res_fit_F=bvec(3,nboot,njack);
  
  par_res_fit_F.data[0]=A;
  par_res_fit_F.data[1]=B;
  par_res_fit_F.data[2]=C;
  
  const char tag_ml[1024]="m\\sl\\N\\SMS,2GeV\\N (GeV)";
  const char tag_a2[1024]="a\\S2\\N (fm)";
  double lat_med_fm[4]={lat[0].med()*hc,lat[1].med()*hc,lat[2].med()*hc,lat[3].med()*hc};
  
  plot_funz_ml(combine("fixQ2_extr/plots/%s_funz_ml.xmg",title).c_str(),title,tag_ml,title,ml,Y,par_res_fit_F,
	       ml_phys.med(),fun_fit_F,Y_chir_cont);
  plot_funz_a2(combine("fixQ2_extr/plots/%s_funz_a2.xmg",title).c_str(),title,tag_a2,title,lat_med_fm,
	       Y_estr_ml,par_res_fit_F,fun_fit_F,Y_chir_cont);
  
  return Y_chir_cont;
}

int icombo(int iens,int ith)
{return iens*nth_to_use+ith;}

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
  for(int iboot=0;iboot<nboot;iboot++) out.data[iboot]=in.data[iboot_jack[iboot]];
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

boot parab_spline(boot *A,double *x,double xint,const char *path=NULL)
{
  boot out(nboot,njack);
  ofstream plot;
  if(path!=NULL) plot.open(path);
  
  boot a(nboot,njack),b(nboot,njack),c(nboot,njack);
  for(int iboot=0;iboot<=nboot;iboot++)
    {
      double y[3];
      for(int i=0;i<3;i++) y[i]=A[i][iboot];
      parabolic_spline(a.data[iboot],b.data[iboot],c.data[iboot],x,y);
      out.data[iboot]=a[iboot]*xint*xint+b[iboot]*xint+c[iboot];
    }
  
  if(path)
    {
      bvec inte(100,nboot,njack);
      double xin=min(xint,min(x[0],x[2])),xfin=max(xint,max(x[2],x[0])),dx=(xfin-xin)/99;
      double xn[100];
      for(int ix=0;ix<100;ix++)
	{
	  xn[ix]=xin+dx*ix;
	  inte[ix]=a*xn[ix]*xn[ix]+b*xn[ix]+c;
	}
      plot<<write_polygon(xn,inte);
      plot<<"&"<<endl;
      plot<<"@type xydy"<<endl;
      for(int i=0;i<3;i++) plot<<x[i]<<" "<<A[i]<<endl;
      plot<<"&"<<endl;
      plot<<xint<<" "<<out<<endl;
      plot<<"&"<<endl;
    }
  
  return out;
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
  FP=F0=FT=EP=Q2=bvec(nens*nth_to_use,nboot,njack);
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

	  remove_pole(FP[ic],F0[ic],FT[ic],Q2[ic]);
          
	  //pass to study FT/FP
	  if(study_T_fr_P) if(ith) FT[ic]/=FP[ic];
	  
	  if(reciproc_data)
            {
              FP[ic]=1/FP[ic];
              F0[ic]=1/F0[ic];
              FT[ic]=1/FT[ic];
            }
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
  
  //find Q2_max_ref
  double Q2_min_ref=0;
  double Q2_max_ref=Q2_max_ph[stranged];
  for(int iens=0;iens<nens;iens++)
    if(use[iens])
      {
	Q2_max_ref=min(Q2_max_ref,Q2[icombo(iens,0)].med());
	Q2_min_ref=max(Q2_min_ref,Q2[icombo(iens,nth_to_use-1)].med());
      }
  //Q2_min_ref*=0;
  //Q2_max_ref*=1.35;
  cout<<"Q2_min_ref: "<<Q2_min_ref<<endl;
  cout<<"Q2_max_ref: "<<Q2_max_ref<<endl;
  
  //prepare value of Q2 ref
  int nth_to_int=8;
  double Q2_ref[nth_to_int];
  bvec F0_chir_cont(nth_to_int,nboot,njack),FP_chir_cont(nth_to_int,nboot,njack),FT_chir_cont(nth_to_int,nboot,njack);
  for(int ith_to_int=0;ith_to_int<nth_to_int;ith_to_int++)
    {
      Q2_ref[ith_to_int]=Q2_max_ref-(Q2_max_ref-Q2_min_ref)/(nth_to_int-1)*ith_to_int;
      
      bvec FP_int(nens,nboot,njack),FT_int(nens,nboot,njack),F0_int(nens,nboot,njack);
      for(int iens=0;iens<nens;iens++)
	if(use[iens])
	  {
	    //find central
	    int ice_0=0,ice_P=1;
	    do ice_0++;
	    while(ice_0<3 && Q2_ref[ith_to_int]<(Q2[icombo(iens,ice_0)].med()+Q2[icombo(iens,ice_0+1)].med())/2);
	    do ice_P++;
	    while(ice_P<3 && Q2_ref[ith_to_int]<(Q2[icombo(iens,ice_P)].med()+Q2[icombo(iens,ice_P+1)].med())/2);
	    
	    double x_0[3]={Q2[icombo(iens,ice_0-1)].med(),Q2[icombo(iens,ice_0)].med(),Q2[icombo(iens,ice_0+1)].med()};
	    double x_P[3]={Q2[icombo(iens,ice_P-1)].med(),Q2[icombo(iens,ice_P)].med(),Q2[icombo(iens,ice_P+1)].med()};
	    FP_int[iens]=parab_spline(FP.data+icombo(iens,ice_P-1),x_P,Q2_ref[ith_to_int],
				      combine("fixQ2_extr/plots/spline_fP_iens_%d_ith_%d.xmg",iens,ith_to_int).c_str());
	    FT_int[iens]=parab_spline(FT.data+icombo(iens,ice_P-1),x_P,Q2_ref[ith_to_int],
				      combine("fixQ2_extr/plots/spline_fT_iens_%d_ith_%d.xmg",iens,ith_to_int).c_str());
	    F0_int[iens]=parab_spline(F0.data+icombo(iens,ice_0-1),x_0,Q2_ref[ith_to_int],
				      combine("fixQ2_extr/plots/spline_f0_iens_%d_ith_%d.xmg",iens,ith_to_int).c_str());
	    
	    //put back kinematical bound
	    //if(ith_to_int==(nth_to_int-1)) F0_int[iens]=FP_int[iens]=(FP_int[iens]+F0_int[iens])/2;
	  }
      
      //perform the fit
      F0_chir_cont[ith_to_int]=fit(ml,F0_int,combine("F0_th_%d",ith_to_int).c_str());
      FP_chir_cont[ith_to_int]=fit(ml,FP_int,combine("FP_th_%d",ith_to_int).c_str());
      FT_chir_cont[ith_to_int]=fit(ml,FT_int,combine("FT_th_%d",ith_to_int).c_str());
      
      if(reciproc_data)
	{
	  FP_chir_cont[ith_to_int]=1/FP_chir_cont[ith_to_int];
	  F0_chir_cont[ith_to_int]=1/F0_chir_cont[ith_to_int];
	  FT_chir_cont[ith_to_int]=1/FT_chir_cont[ith_to_int];
	}

      //convert back FT from FT/FP
      if(study_T_fr_P) FT_chir_cont[ith_to_int]*=FP_chir_cont[ith_to_int];

      put_pole(FP_chir_cont[ith_to_int],F0_chir_cont[ith_to_int],FT_chir_cont[ith_to_int],Q2_ref[ith_to_int]);
    }
  
  //plot data
  {
    ofstream out0("fixQ2_extr/plots/F0_chir_cont.xmg");
    ofstream outP("fixQ2_extr/plots/FP_chir_cont.xmg");
    ofstream outT("fixQ2_extr/plots/FT_chir_cont.xmg");
    out0<<"@type xydy"<<endl;
    outP<<"@type xydy"<<endl;
    outT<<"@type xydy"<<endl;
    for(int ith_to_int=0;ith_to_int<nth_to_int;ith_to_int++)
      {
	out0<<Q2_ref[ith_to_int]<<" "<<F0_chir_cont[ith_to_int]<<endl;
	outP<<Q2_ref[ith_to_int]<<" "<<FP_chir_cont[ith_to_int]<<endl;
	outT<<Q2_ref[ith_to_int]<<" "<<FT_chir_cont[ith_to_int]/FP_chir_cont[ith_to_int]<<endl;
      }
  }
  
  //save data
  char result_path[100]="fixQ2_extr/results_Q2_F0_FP_FT";
  {
    bvec temp(nth_to_int,nboot,njack);
    for(int ith=0;ith<nth_to_int;ith++) temp[ith]=Q2_ref[ith];
    temp.write_to_binfile(result_path);
  }
  F0_chir_cont.append_to_binfile(result_path);
  FP_chir_cont.append_to_binfile(result_path);
  FT_chir_cont.append_to_binfile(result_path);
  
  //vdm test
  boot ghat(nboot,njack),fdstar(nboot,njack);
  ghat.fill_gauss(15.8,0.8,198);
  fdstar.fill_gauss(0.278,0.016,575);
  {
    double xout[100];
    bvec yout(100,nboot,njack);
    double dx=Q2_max_ph[0]/99;
    for(int iq=0;iq<100;iq++)
      {
	xout[iq]=dx*iq;
	yout[iq]=0.5*MDstar*ghat*fdstar/(MDstar*MDstar-xout[iq]);
      }
    ofstream plot("fixQ2_extr/plots/vdm.xmg");
    plot<<write_polygon(xout,yout);
    plot<<"&"<<endl;
    plot<<write_ave_line(xout,yout)<<endl;
  }
  cout<<0.5*MDstar*ghat*fdstar<<" "<<MDstar*MDstar<<endl;
  
  return 0;
}
