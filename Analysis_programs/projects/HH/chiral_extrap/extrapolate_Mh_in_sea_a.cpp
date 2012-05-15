#include "../HH_common.cpp"

const int nbeta=4;
char **base_corrs_path,**ens_name;
int nens,*T,*ibeta,*nmass;
double **mass,*sea_mass;

//read data
bvec *Mh;

bvec ml;
int ref_ml_beta[4]={-1,-1,-1,-1};

bvec par_res_fit_Mh;

const char set_color[nbeta][1024]={"black","blue","red","green4"};
const char set_fill_color[nbeta][1024]={"grey","turquoise","yellow","green"};
const char set_symbol[nbeta][1024]={"square","circle","triangle up","triangle left"};
const char set_legend[nbeta][1024]={"\\xb\\0=3.80","\\xb\\0=3.90","\\xb\\0=4.05","\\xb\\0=4.20"};
const char set_legend_fm[nbeta][1024]={"a = 0.098 fm","a = 0.085 fm","a = 0.067 fm","a = 0.054 fm"};

int plot_iboot;
int include_a4;
int include_380;

double fun_fit_Mh(double A,double B,double C,double D,double ml,double a)
{
  double a1=a/lat[1][plot_iboot];
  return A*(1 + B*ml + C*a1*a1 + D*a1*a1*a1*a1);
}

void plot_funz_ml(const char *out_path,const char *title,const char *xlab,const char *ylab,bvec &X,bvec &Y,bvec &par,double X_phys,double (*fun)(double,double,double,double,double,double),boot &chiral_extrap_cont)
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
	    for(int iboot=plot_iboot=0;iboot<nboot+1;plot_iboot=iboot++)
	      Y_pol.data[ipoint].data[iboot]=fun(par[0][iboot],par[1][iboot],par[2][iboot],par[3][iboot],X_pol[ipoint],lat[ib][iboot]);
	  
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
	  Y_pol.data[ipoint]=fun(par[0][iboot],par[1][iboot],par[2][iboot],par[3][iboot],X_pol[ipoint],0);
      //out.set(1,"magenta");
      //out.polygon(X_pol,Y_pol);
      //out.new_set();
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
      for(int iens=0;iens<nens;iens++) if(ibeta[iens]==ib) out.fout<<X[iens].med()<<" "<<Y.data[iens]<<endl;
      out.new_set();
    }
  
  //plot the extrapolated point with error
  out.set(4,"none","circle","indigo","filled");
  out.set_line_size(3);
  out.set_legend("Physical point");
  out.print_graph(X_phys,chiral_extrap_cont);
  out.new_set();
}

bool fitting_fK;
double *X_fit;
double *Y_fit,*err_Y_fit;

//calculate the chi square
double chi2(double A,double B,double C,double D,double *a)
{
  double ch2=0;
  
  for(int iens=0;iens<nens;iens++)
    if(include_380 || ibeta[iens]!=0)
      ch2+=pow((Y_fit[iens]-fun_fit_Mh(A,B,C,D,X_fit[iens],a[ibeta[iens]]))/err_Y_fit[iens],2);
  
  return ch2;
}

//wrapper for the calculation of the chi2
void chi2_wr(int &npar,double *fuf,double &ch,double *p,int flag)
{
  double A=p[0],B=p[1],C=p[2],D=p[3];
  double *a=p+4;
  ch=chi2(A,B,C,D,a);
}

void fit(boot &A,boot &B,boot &C,boot &D,bvec &X,bvec &Y)
{
  //copy X
  X_fit=new double[nens];
  for(int iens=0;iens<nens;iens++) X_fit[iens]=X[iens].med();
  Y_fit=new double[nens];
  err_Y_fit=new double[nens];
  
  TMinuit minu;
  minu.SetPrintLevel(-1);
  
  int npars=4;
  minu.DefineParameter(0,"A",0.0,0.0001,0,0);
  minu.DefineParameter(1,"B",0.0,0.0001,0,0);
  minu.DefineParameter(2,"C",0.0,0.0001,0,0);
  minu.DefineParameter(3,"D",0.0,0.0001,0,0);
  if(!include_a4)
    {
      minu.FixParameter(3);
      npars--;
    }
  minu.SetFCN(chi2_wr);
  
  double C2;
  for(int iboot=0;iboot<nboot+1;iboot++)
    {
      if(iboot>0)
        minu.SetPrintLevel(-1);
      
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
          Y_fit[iens]=Y.data[iens].data[iboot];
          err_Y_fit[iens]=Y.data[iens].err();
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
      if(iboot==0) C2=chi2(A.data[iboot],B[iboot],C[iboot],D[iboot],lat_med);
    }
  
  int ninc_ens=0;
  for(int iens=0;iens<nens;iens++)
    if(ibeta[iens]!=0 || include_380) ninc_ens++;
  
  //calculate the chi2
  cout<<"A=("<<A<<"), B=("<<B<<"), C=("<<C<<"), D=("<<D<<")"<<endl;
  cout<<"Chi2 = "<<C2<<" / "<<ninc_ens-npars<<" = "<<C2/(ninc_ens-npars)<<endl;
  
  delete[] X_fit;
  delete[] Y_fit;
  delete[] err_Y_fit;
}


void plot_funz_a2(const char *out_path,const char *title,const char *xlab,const char *ylab,double *X,bvec &Y,bvec &par,double (*fun)(double,double,double,double,double,double),boot &chiral_extrap_cont)
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
	  Y_pol.data[iens].data[iboot]=fun(par[0][iboot],par[1][iboot],par[2][iboot],par[3][iboot],ml_phys[iboot],X_pol[iens]/hc);
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

int main(int narg,char **arg)
{
  int set_color[4]={1,4,2,15};
  int set_color2[4]={7,9,11,3};
  int set_symbol[4]={2,1,4,5};
  int set_linewitdth[4]={1,2,2,2};
  char set_name[4][10]={"3.80","3.90","4.05","4.20"};

  init_latpars();
  
  //read ensemble list, meson masses and meson name
  FILE *an_input_file=open_file("analysis_pars","r");
  char ens_list_path[1024],meson_mass_file[1024],meson_name[1024];
  read_formatted_from_file_expecting(ens_list_path,an_input_file,"%s","ens_list_path");
  read_formatted_from_file_expecting(meson_mass_file,an_input_file,"%s","meson_mass_file");
  read_formatted_from_file_expecting(meson_name,an_input_file,"%s","meson_name");
  read_formatted_from_file_expecting((char*)&include_a4,an_input_file,"%d","include_a4");
  read_formatted_from_file_expecting((char*)&include_380,an_input_file,"%d","include_380");
  fclose(an_input_file);
  
  //load ensembles list and parameters
  load_ensembles_list(base_corrs_path,ens_name,nens,T,ibeta,nmass,mass,sea_mass,ens_list_path);
  
  //define ml and ref ml
  ml=bvec(nens,nboot,njack);
  for(int iens=0;iens<nens;iens++)
    {
      int b=ibeta[iens],r=ref_ml_beta[b];
      //define ml
      cout<<iens<<" "<<b<<" "<<sea_mass[iens]<<endl;
      ml[iens]=sea_mass[iens]/lat[b]/Zp[b];
      //set the lighter mass
      if(r==-1||fabs(ml[r].med()-0.050)>fabs(ml[iens].med()-0.050)) ref_ml_beta[b]=iens;
    }
  cout<<"---"<<endl;
  for(int ib=0;ib<nbeta;ib++) if(ref_ml_beta[ib]!=-1) cout<<"Ref "<<ib<<" = "<<ref_ml_beta[ib]<<", "<<ml[ref_ml_beta[ib]]<<" MeV"<<endl;
  cout<<"---"<<endl;

  //load data
  int nm=11;
  Mh=new bvec[nm];
  bvec Mh_chir_cont(nm,nboot,njack);
  for(int im=0;im<nm;im++)
    {
      Mh[im]=bvec(nens,nboot,njack);
      Mh[im].load(meson_mass_file,im);
      
      //perform the fit
      boot A(nboot,njack),B(nboot,njack),C(nboot,njack),D(nboot,njack);
      fit(A,B,C,D,ml,Mh[im]);
      
      //chiral extrapolation
      boot Mh_chir[nbeta];
      bvec Mh_estr_ml(nbeta,nboot,njack);
      for(int ib=0;ib<nbeta;ib++)
	{
	  Mh_chir[ib]=boot(nboot,njack);
	  for(int iboot=0;iboot <nboot+1;iboot++)
	    {
	      int r=ref_ml_beta[ib];
	      Mh_chir_cont[im].data[iboot]=fun_fit_Mh(A[iboot],B[iboot],C[iboot],D[iboot],ml_phys[iboot],0);
	      Mh_chir[ib].data[iboot]=fun_fit_Mh(A[iboot],B[iboot],C[iboot],D[iboot],ml_phys[iboot],lat[ib][iboot]);
	      
	      if(r!=-1)
		Mh_estr_ml.data[ib].data[iboot]=Mh[im][r][iboot]*fun_fit_Mh(A[nboot],B[nboot],C[nboot],D[nboot],ml_phys[nboot],0)/fun_fit_Mh(A[nboot],B[nboot],C[nboot],D[nboot],ml[r][nboot],0);
	    }
	}
      
      //chiral and continuum
      cout<<"Mh = ("<<Mh_chir_cont[im]*1000<<") MeV"<<endl;
      cout<<endl;

      par_res_fit_Mh=bvec(4,nboot,njack);
      
      par_res_fit_Mh.data[0]=A;
      par_res_fit_Mh.data[1]=B;
      par_res_fit_Mh.data[2]=C;
      par_res_fit_Mh.data[3]=D;
      
      const char tag_ml[1024]="m\\sl\\N\\SMS,2GeV\\N (GeV)";
      const char tag_a2[1024]="a\\S2\\N (fm)";
      double lat_med_fm[4]={lat[0].med()*hc,lat[1].med()*hc,lat[2].med()*hc,lat[3].med()*hc};

      plot_funz_ml(combine("Mh_funz_ml_%d.xmg",im).c_str(),meson_name,tag_ml,meson_name,ml,Mh[im],par_res_fit_Mh,ml_phys.med(),fun_fit_Mh,Mh_chir_cont[im]);
      plot_funz_a2(combine("Mh_funz_a2_%d.xmg",im).c_str(),meson_name,tag_a2,meson_name,lat_med_fm,Mh_estr_ml,par_res_fit_Mh,fun_fit_Mh,Mh_chir_cont[im]);
    }
  
  Mh_chir_cont.write_to_binfile("results");
  
  /////////////////////////////////////////////////////////////////////////////////////////////////
  
  double min_rep_hmass=0.7,min_rep_hmass_bare[4];
  
  //write the charm at beta 3.90
  ofstream charm_beta_390("charm_beta_390.xmg");
  charm_beta_390<<"@type xydy"<<endl;
  for(int iens=0;iens<nens;iens++)
    if(ibeta[iens]==1)
      charm_beta_390<<sea_mass[iens]<<" "<<Mh[1][iens]<<endl;
  
  //prepare reciprocal mass and substets
  double rep_hmass[nm+1];
  int nmin_rep_hmass=0;
  for(int im=0;im<nm;im++)
    {
      rep_hmass[im]=1/ref_hmass[im].med();
      if(rep_hmass[im]>=min_rep_hmass) nmin_rep_hmass++;
    }
  double rep_min_hmass[nmin_rep_hmass];
  double rep_hmass_bare[4][nref_hmass];
  double rep_min_hmass_constrained[nmin_rep_hmass+1];
  int itemp=0;
  bvec Mh_min_chir_cont_constrained(nmin_rep_hmass+1,nboot,njack);
  for(int im=0;im<nm;im++)
    {
      if(rep_hmass[im]>=min_rep_hmass)
	{
	  rep_min_hmass_constrained[itemp]=rep_min_hmass[itemp]=rep_hmass[im];
	  Mh_min_chir_cont_constrained[itemp]=Mh_chir_cont[im];
	  itemp++;
	}
      for(int ib=0;ib<4;ib++)
	rep_hmass_bare[ib][im]=rep_hmass[im]/(lat[ib].med()*Zp[ib].med());
    }
  
  for(int ib=0;ib<4;ib++)
    min_rep_hmass_bare[ib]=min_rep_hmass/(lat[ib].med()*Zp[ib].med());
  
  Mh_min_chir_cont_constrained[nmin_rep_hmass].fill_gauss(0,0.000000001,22332132);
  rep_min_hmass_constrained[nmin_rep_hmass]=0;
  
  //perform static limit extrapolation at fixed lattice spacing
  bvec pars_linear_fit_static_fixlatt[4];
  bvec pars_linear_fit_static_fixlatt_constrained[4];
  bvec pars_parab_fit_static_fixlatt[4];
  bvec pars_parab_fit_static_fixlatt_constrained[4];
  bvec pars_parab_fit_static_fixlatt_bare[4];
  for(int ib=0;ib<4;ib++)
    {
      //copy Mh for the ref ensemble
      bvec y(nm,nboot,njack),yconst(nmin_rep_hmass+1,nboot,njack);
      int itemp=0;
      for(int im=0;im<nm;im++)
	{
	  y[itemp]=Mh[im][ref_ml_beta[ib]];
	  if(rep_hmass[im]>=min_rep_hmass)
	    {
	      yconst[itemp]=Mh[im][ref_ml_beta[ib]];
	      itemp++;
	    }
	}
      yconst[nmin_rep_hmass].fill_gauss(0,0.000000001,22332132);
      pars_linear_fit_static_fixlatt[ib]=poly_fit(rep_hmass,y,1,min_rep_hmass,1.1);
      pars_parab_fit_static_fixlatt[ib]=poly_fit(rep_hmass,y,2,min_rep_hmass,1.1);
      pars_parab_fit_static_fixlatt_bare[ib]=poly_fit(rep_hmass_bare[ib],y,2,min_rep_hmass_bare[ib],rep_hmass_bare[ib][0]);
      pars_linear_fit_static_fixlatt_constrained[ib]=poly_fit(rep_min_hmass_constrained,yconst,1,0,1.1);
      pars_parab_fit_static_fixlatt_constrained[ib]=poly_fit(rep_min_hmass_constrained,yconst,2,0,1.1);
      cout<<"Mh in the static limit at fixed lattice spacing "<<lat[ib]<<" = "<<pars_linear_fit_static_fixlatt[ib][0]<<endl;
    }
  
  //perform static limit in the continuum
  bvec pars_linear_fit_static=poly_fit(rep_hmass,Mh_chir_cont,1,min_rep_hmass,1.1);
  bvec pars_linear_fit_static_constrained=poly_fit(rep_min_hmass_constrained,Mh_min_chir_cont_constrained,1,0,1.1);
  cout<<"Static limit linear extrapolation in the continuum: "<<pars_linear_fit_static[0]<<endl;
  bvec pars_parab_fit_static=poly_fit(rep_hmass,Mh_chir_cont,2,min_rep_hmass,1.1);
  bvec pars_parab_fit_static_constrained=poly_fit(rep_min_hmass_constrained,Mh_min_chir_cont_constrained,2,0,1.1);
  cout<<"Static limit parab extrapolation in the continuum: "<<pars_parab_fit_static[0]<<endl;  
  
  //scaling in bare units
  {
    int iset=0;
    
    ofstream static_limit_mh_bare("scaling_static_limit_mh_bare.xmg");
    
    //x-axis
    static_limit_mh_bare<<"@    xaxis  label \"1/m\\sh\\N\\S\\oMS\\O,2 GeV\\N [GeV\\S-1\\N]\""<<endl<<
      "@    xaxis  label char size 1.390000"<<endl<<
      "@    xaxis  tick major 0.2"<<endl;
    
    //y-axis
    static_limit_mh_bare<<"@    yaxis  label \"M\\SVect\\N\\shh\\N-M\\SPseudo\\N\\shh\\N [GeV]\""<<endl<<
      "@    yaxis  label char size 1.390000"<<endl<<
      "@    yaxis  tick major 0.05"<<endl;
    
    //plot the lines at fixed lattice spacing
    for(int ib=0;ib<4;ib++)
      {
	static_limit_mh_bare<<
	  "@    s"<<iset<<" fill type 1"<<endl<<
	  "@    s"<<iset<<" fill color "<<set_color2[ib]<<endl<<
	  "@    s"<<iset<<" line color "<<set_color2[ib]<<endl<<
	  "@    s"<<iset+1<<" line color "<<set_color[ib]<<endl;
	static_limit_mh_bare<<write_poly_with_error(pars_parab_fit_static_fixlatt_bare[ib],0,rep_hmass_bare[ib][0]*1.1);
	iset+=2;
      }
    
    //write representative ensembles
    for(int ib=0;ib<4;ib++)
      for(int istep=0;istep<2;istep++)
	{
	  static_limit_mh_bare<<"@type xydy"<<endl<<
	    "@ s"<<iset<<" symbol color "<<set_color[ib]<<endl<<
	    "@ s"<<iset<<" symbol "<<set_symbol[ib]<<endl<<
	    "@ s"<<iset<<" errorbar color "<<set_color[ib]<<endl<<
	    "@ s"<<iset<<" line type 0"<<endl;
	  if(istep==1)
	    static_limit_mh_bare<<"@ s"<<iset<<" symbol linewidth 2"<<endl<<
	      "@ s"<<iset<<" errorbar linewidth 2"<<endl<<
	      "@ s"<<iset<<" legend \"\\xb="<<set_name[ib]<<"\\N\""<<endl<<
	      "@ s"<<iset<<" errorbar riser linewidth 2"<<endl;
	  
	  int lim=(istep==0) ? nref_hmass : nmin_rep_hmass;
	  
	  for(int im=0;im<lim;im++) static_limit_mh_bare<<rep_hmass_bare[ib][im]<<" "<<Mh[im][ref_ml_beta[ib]]<<endl;
	  static_limit_mh_bare<<"&"<<endl;
	  iset++;
	}
  }
	  
  //plot the four representative ensemble as a function of 1/m_h with various possibilities
  for(int istep=0;istep<16;istep++)
    {
      int iset=0;
      ifstream header("header.xmg");
      ofstream static_limit(combine("scaling_static_limit_%02d.xmg",istep).c_str());
      
      char a[10000];
      while(header.getline(a,10000)) static_limit<<a<<endl;
      
      //title
      static_limit<<"@    title \"Vector-Pseudoscalar quarkonium mass difference\""<<endl;
      
      //zoom
      static_limit<<
	"@    world 0, 0, 1.07, 0.24"<<endl<<
	"@    view 0.150000, 0.150000, 1.150000, 0.850000"<<endl;
      
      //legend position
      static_limit<<"@    legend 0.975, 0.4"<<endl;
      
      //x-axis label of charm and bottom
      static_limit<<"@with string"<<endl<<
	"@    string on"<<endl<<
	"@    string loctype view"<<endl<<
	"@    string 0.943627450981, 0.115196078431"<<endl<<
	"@    string color 13"<<endl<<
	"@    string rot 0"<<endl<<
	"@    string font 0"<<endl<<
	"@    string just 0"<<endl<<
	"@    string char size 1.210000"<<endl<<
	"@    string def \"Charm\""<<endl<<
	"@with string"<<endl<<
	"@    string on"<<endl<<
	"@    string loctype view"<<endl<<
	"@    string 0.359068627451, 0.118872549019"<<endl<<
	"@    string color 13"<<endl<<
	"@    string rot 0"<<endl<<
	"@    string font 0"<<endl<<
	"@    string just 0"<<endl<<
	"@    string char size 1.210000"<<endl<<
	"@    string def \"Bottom\""<<endl;
      
      //x-axis
      static_limit<<"@    xaxis  label \"1/m\\sh\\N\\S\\oMS\\O,2 GeV\\N [GeV\\S-1\\N]\""<<endl<<
	"@    xaxis  label char size 1.390000"<<endl<<
	"@    xaxis  tick major 0.2"<<endl;
      
      //y-axis
      static_limit<<"@    yaxis  label \"M\\SVect\\N\\shh\\N-M\\SPseudo\\N\\shh\\N [GeV]\""<<endl<<
	"@    yaxis  label char size 1.390000"<<endl<<
	"@    yaxis  tick major 0.05"<<endl;
      
      //linear static extrapolation at fixed lattice spacing
      if(istep>=4)
	{
	  bvec *pars,par;
	  switch(istep)
	    {
	    case 4:
	      pars=pars_linear_fit_static_fixlatt;
	      par=pars_linear_fit_static;
	      break;
	    case 5:
	      pars=pars_linear_fit_static_fixlatt;
	      par=pars_linear_fit_static_constrained;
	      break;
	    case 6:
	      pars=pars_parab_fit_static_fixlatt;
	      par=pars_parab_fit_static;
	      break;
	    default:
	      pars=pars_parab_fit_static_fixlatt;
	      par=pars_parab_fit_static_constrained;
	      break;
	    }
	  
	  //plot the lines at fixed lattice spacing
	  for(int ib=0;ib<4;ib++)
	    {
	      static_limit<<
		"@    s"<<iset<<" fill type 1"<<endl<<
		"@    s"<<iset<<" fill color "<<set_color2[ib]<<endl<<
		"@    s"<<iset<<" line color "<<set_color2[ib]<<endl<<
		"@    s"<<iset+1<<" line color "<<set_color[ib]<<endl;
	      static_limit<<write_poly_with_error(pars[ib],0,1.1);
	      iset+=2;
	    }
	  
	  //constrained fit in the continuum
	  static_limit<<
	    "@    s"<<iset<<" fill type 1"<<endl<<
	    "@    s"<<iset<<" fill color 6"<<endl<<
	    "@    s"<<iset<<" line color 6"<<endl<<
	    "@    s"<<iset+1<<" line color 8"<<endl;
	  static_limit<<write_poly_with_error(par,0,1.1);
	  iset+=2;
	}
      
      //write representative ensembles data used for fit (or, only charm)
      if(istep>=1)
	for(int ib=0;ib<4;ib++)
	  {
	    static_limit<<"@type xydy"<<endl<<
	      "@ s"<<iset<<" symbol color "<<set_color[ib]<<endl<<
	      "@ s"<<iset<<" symbol "<<set_symbol[ib]<<endl<<
	      "@ s"<<iset<<" errorbar color "<<set_color[ib]<<endl<<
	      "@ s"<<iset<<" legend \"\\xb="<<set_name[ib]<<"\\N\""<<endl<<
	      "@ s"<<iset<<" line type 0"<<endl<<
	      "@ s"<<iset<<" symbol linewidth "<<set_linewitdth[ib]<<endl<<
	      "@ s"<<iset<<" errorbar linewidth "<<set_linewitdth[ib]<<endl<<
	      "@ s"<<iset<<" errorbar riser linewidth "<<set_linewitdth[ib]<<endl;
	    //if istep==1 or 2 plot only charm
	    if(istep==1 || istep==2) static_limit<<rep_hmass[1]<<" "<<Mh[1][ref_ml_beta[ib]]<<endl;
	    else for(int im=0;im<nmin_rep_hmass;im++) static_limit<<rep_min_hmass[im]<<" "<<Mh[im][ref_ml_beta[ib]]<<endl;
	    static_limit<<"&"<<endl;
	    iset++;
	  }
      
      //plot remaining data
      if(istep>=8)
	for(int ib=0;ib<4;ib++)
	  {
	    static_limit<<"@type xydy"<<endl<<
	      "@ s"<<iset<<" symbol color "<<set_color[ib]<<endl<<
	      "@ s"<<iset<<" symbol "<<set_symbol[ib]<<endl<<
	      "@ s"<<iset<<" errorbar color "<<set_color[ib]<<endl<<
	      "@ s"<<iset<<" line type 0"<<endl<<
	      "@ s"<<iset<<" symbol linewidth "<<min((double)set_linewitdth[ib],1.5)<<endl<<
	      "@ s"<<iset<<" errorbar linewidth "<<min((double)set_linewitdth[ib],1.5)<<endl<<
	      "@ s"<<iset<<" errorbar riser linewidth "<<min((double)set_linewitdth[ib],1.5)<<endl;
	    for(int im=nmin_rep_hmass;im<min(nm,istep-4);im++) static_limit<<rep_hmass[im]<<" "<<Mh[im][ref_ml_beta[ib]]<<endl;
	    static_limit<<"&"<<endl;
	    iset++;
	  }
      
      //draw charm vertical line
      static_limit<<"@type xy"<<endl<<"@    s"<<iset<<" line linestyle 3"<<endl
		  <<"@    s"<<iset<<" line linewidth 2.0"<<endl
		  <<"@    s"<<iset<<" line color 6"<<endl;
      static_limit<<1.005/ref_hmass[1].med()<<" "<<0<<endl;
      static_limit<<1.005/ref_hmass[1].med()<<" "<<0.3<<endl;
      static_limit<<"&"<<endl;
      iset++;
      
      //draw bottom vertical line
      static_limit<<"@    s"<<iset<<" line linestyle 3"<<endl
		  <<"@    s"<<iset<<" line linewidth 2.0"<<endl
		  <<"@    s"<<iset<<" line color 6"<<endl;
      static_limit<<1.015/ref_hmass[nm-1].med()<<" "<<0<<endl;
      static_limit<<1.015/ref_hmass[nm-1].med()<<" "<<0.3<<endl;
      static_limit<<"&"<<endl;
      iset++;
      
      //physical values
      static_limit<<"@    s"<<iset<<" line type 0"<<endl<<
	"@    s"<<iset<<" symbol 10"<<endl<<
	"@    s"<<iset<<" legend \"Experimental\""<<endl<<
	"@    s"<<iset<<" symbol color 10"<<endl<<
	"@    s"<<iset<<" errorbar color 10"<<endl<<
	"@    s"<<iset<<" symbol linewidth 2"<<endl<<
	"@    s"<<iset<<" errorbar linewidth 2"<<endl<<
	"@    s"<<iset<<" errorbar riser linewidth 2"<<endl<<
	"@type xydy"<<endl<<
	"0.867 0.1164 0"<<endl<<
	"0.20245 0.065 0.008"<<endl<<
	"&"<<endl;
      iset++;
      
      //continuum limit mass per mass on fitted data
      if(istep>=2)
	{
	  static_limit<<"@type xydy"<<endl;
	  static_limit<<
	    "@    s"<<iset<<" symbol 3"<<endl<<
	    "@    s"<<iset<<" symbol color 8"<<endl<<
	    "@    s"<<iset<<" symbol fill 1 "<<endl<<
	    "@    s"<<iset<<" symbol fill color 8 "<<endl<<
	    "@    s"<<iset<<" legend \"Extrapol.\""<<endl<<
	    "@    s"<<iset<<" symbol linewidth 2"<<endl<<
	    "@    s"<<iset<<" line type 0"<<endl<<
	    "@    s"<<iset<<" errorbar linewidth 2"<<endl<<
            "@    s"<<iset<<" errorbar color 8"<<endl<<
	    "@    s"<<iset<<" errorbar riser linewidth 2"<<endl;
	  if(istep==2) static_limit<<rep_hmass[1]<<" "<<Mh_chir_cont[1]<<endl;
	  else for(int im=0;im<nmin_rep_hmass;im++) static_limit<<rep_hmass[im]<<" "<<Mh_chir_cont[im]<<endl;
	  static_limit<<"&"<<endl;
	  iset++;
	}
      //remaining data
      if(istep>=8)
	{
	  static_limit<<"@type xydy"<<endl;
	  static_limit<<
	    "@    s"<<iset<<" symbol 3"<<endl<<
	    "@    s"<<iset<<" symbol color 8"<<endl<<
	    "@    s"<<iset<<" symbol linewidth 1"<<endl<<
	    "@    s"<<iset<<" line type 0"<<endl<<
	    "@    s"<<iset<<" errorbar linewidth 1"<<endl<<
	    "@    s"<<iset<<" errorbar color 8"<<endl<<
	    "@    s"<<iset<<" errorbar riser linewidth 1"<<endl;
	  for(int im=nmin_rep_hmass;im<min(nm,istep-4);im++) static_limit<<rep_hmass[im]<<" "<<Mh_chir_cont[im]<<endl;
	  static_limit<<"&"<<endl;
	  iset++;
	}
      
      //constrained parabolic result
      if(istep>=7)
	{      
	  boot ex=pars_parab_fit_static_constrained[0]+1/ref_hmass[nm-1]*(pars_parab_fit_static_constrained[1]+1/ref_hmass[nm-1]*pars_parab_fit_static_constrained[2]);
	  static_limit<<
	    "@with line"<<endl<<
	    "@    line on"<<endl<<
	    "@    line loctype view"<<endl<<
	    "@    line 0.348039215686, 0.240196078431, 0.416666666667, 0.226715686275"<<endl<<
	    "@    line linewidth 2.0"<<endl<<
	    "@    line linestyle 1"<<endl<<
	    "@    line color 1"<<endl<<
	    "@    line arrow 1"<<endl<<
	    "@    line arrow type 0"<<endl<<
	    "@    line arrow length 1.000000"<<endl<<
	    "@    line arrow layout 1.000000, 1.000000"<<endl<<
	    "@line def@with string"<<endl<<
	    "@    string on"<<endl<<
	    "@    string loctype view"<<endl<<
	    "@    string 0.421568627451, 0.218137254902"<<endl<<
	    "@    string color 1"<<endl<<
	    "@    string rot 0"<<endl<<
	    "@    string font 0"<<endl<<
	    "@    string just 0"<<endl<<
	    "@    string char size 1.210000"<<endl<<
	    "@    string def \""<<int(ex.med()*1000)<<" ("<<int(ex.err()*1000)<<") MeV\""<<endl;
	}
    }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  return 0;
}
