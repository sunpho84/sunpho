#include "include.h"

double nucleonp=938.27;
double nucleon0=939.57;

double xi0=1314.86;
double xim=1321.31;

double sigmap=1189.37;
double sigmam=1197.449;
double sigma0=1192.642;

double a[3]={0.0885/197,0.0815/197,0.0619/197};
int ib;

int tmin_2pts;
int tmax_2pts;

int tmin_mass_slope;
int tmax_mass_slope;

int tmin_elec_slope=4;
int tmax_elec_slope=10;

int tmin_sigma_elec_slope=4;
int tmax_sigma_elec_slope=10;

static const int ext_islsrc=1,ext_islsink=0;
struct mega_header_t
{
  int twist,nf,nsrc;
  int T,LX,LY,LZ;
  int nk;
  int nslsrc,nslsink;
  double ksm,beta,ksea,musea,csw;
  double k[12];
  int sstepsrc[3];
  int size;
  int icorr(int im1,int im2,int im3,int parity,int islsrc=ext_islsrc,int islsink=ext_islsink)
  {return islsrc+nslsrc*(islsink+nslsink*(parity+2*(im1+nk*(im2+nk*im3))));}
};

int T,TH,njacks=15;
mega_header_t mega_header;
const double eu=2.0/3,ed=-1.0/3;
const double e2=4*M_PI/137;

void separator()
{cout<<"------------------------------"<<endl;}

void read_analysis_pars()
{
  FILE *fin=open_file("analysis_pars","r");
  read_formatted_from_file_expecting((char*)&ib,fin,"%d","ib");
  read_formatted_from_file_expecting((char*)&TH,fin,"%d","L");
  T=2*TH;
  read_formatted_from_file_expecting((char*)&tmin_2pts,fin,"%d","tmin_2pts");
  read_formatted_from_file_expecting((char*)&tmax_2pts,fin,"%d","tmax_2pts");
  read_formatted_from_file_expecting((char*)&tmin_mass_slope,fin,"%d","tmin_mass_slope");
  read_formatted_from_file_expecting((char*)&tmax_mass_slope,fin,"%d","tmax_mass_slope");
  fclose(fin);
}

////////////////////////////////////////// hadrons //////////////////////////////////////

struct hadr_t
{
  int m1,m2,m3;
  int n1,n2,n3;
  int a(){return mega_header.icorr(m1,m2,m3,0);}
  int b(){return mega_header.icorr(n1,n2,n3,0);}
  int c(){return mega_header.icorr(m1,m2,m3,1);}
  int d(){return mega_header.icorr(n1,n2,n3,1);}
  hadr_t(int m1,int m2,int m3) : m1(m1),m2(m2),m3(m3)
    {
      n1=(m1/2)*2+!(m1%2);
      n2=(m2/2)*2+!(m2%2);
      n3=(m3/2)*2+!(m3%2);
    }
 private:
  hadr_t() {}
};

hadr_t nuc_d_000(0,0,0);
hadr_t nuc_d_001(0,0,1);
hadr_t nuc_d_010(0,1,0);
hadr_t nuc_d_011(0,1,1);
hadr_t nuc_d_100(1,0,0);
hadr_t nuc_d_101(1,0,1);
hadr_t nuc_d_110(1,1,0);
hadr_t nuc_d_111(1,1,1);

hadr_t nuc_d(0,1,1);
hadr_t sig_d(0,3,0); //or 031?
hadr_t xi_d(2,1,2); //or 213?

/////////////////////////////// loading from file //////////////////////////////

struct record_t
{
  int ntot_combo;
  int size;
  double *vect;
  record_t(mega_header_t head)
  {
    int nparity_combo=2,nsme_combo=head.nslsrc*head.nslsink;
    int nk_combo=head.nk*head.nk*head.nk;
    ntot_combo=nparity_combo*nk_combo*nsme_combo;
    size=ntot_combo*head.T*2*(njacks+1);
    vect=new double[size];
  }
  
  //read jacknives
  void read_vect(FILE *fin)
  {
    int iconf;
    int rc1=fread(&iconf,sizeof(int),1,fin);
    if(rc1!=1) crash("otbained %d while reading iconf expecting %d",rc1,1);
    cout<<iconf<<endl;
    int rc2=fread(vect,sizeof(double),size,fin);
    if(rc2!=size) crash("otbained %d while reading vect expecting %d",rc2,size);
  }
  ~record_t() {delete [] vect;}
private:
  record_t(){}
};

ostream& operator<<(ostream &out,mega_header_t &in)
{
  out<<"Twist: "<<in.twist<<endl;
  out<<"Nf: "<<in.nf<<endl;
  out<<"Nsrc: "<<in.nsrc<<endl;
  out<<"T: "<<in.T<<endl;
  out<<"LX: "<<in.LX<<endl;
  out<<"LY: "<<in.LY<<endl;
  out<<"LZ: "<<in.LZ<<endl;
  out<<"Nk: "<<in.nk<<endl;
  out<<"Nslsrc: "<<in.nslsrc<<endl;
  out<<"Nslsink: "<<in.nslsink<<endl;
  out<<"Ksm: "<<in.ksm<<endl;
  out<<"Beta: "<<in.beta<<endl;
  out<<"Ksea: "<<in.ksea<<endl;
  out<<"Musea: "<<in.musea<<endl;
  out<<"Csw: "<<in.csw<<endl;
  for(int i=0;i<2*in.nk;i++) cout<<"K["<<i<<"]: "<<in.k[i]<<endl;
  for(int i=0;i<in.nslsrc+in.nslsink;i++) cout<<"Sstepsrc["<<i<<"]: "<<in.sstepsrc[i]<<endl;
  out<<"Size: "<<in.size<<endl;

  return out;
}

////////////////////////////////////// holding various contractions ////////////////////////

struct corr_t
{
  jvec *data;
  corr_t(const char *path,double coeff=+1)
  {
    FILE *fin=open_file(path,"r");

    //load header
    int rc=fread(&mega_header,sizeof(mega_header_t),1,fin);
    if(rc!=1) crash("obtained %d while reading from path: %s",rc,path);
    cout<<mega_header<<endl;
    
    //read the record
    record_t *record=new record_t(mega_header);
    data=(jvec*)malloc(record->ntot_combo*sizeof(jvec));
    if(data==NULL) crash("allocating: %d",record->ntot_combo*sizeof(jvec));
  
    //reading vector
    record->read_vect(fin);
    for(int icombo=0;icombo<record->ntot_combo;icombo++) data[icombo].create(T,njacks);
    for(int ijack=0;ijack<=njacks;ijack++)
      {
	for(int icombo=0;icombo<record->ntot_combo;icombo++)
	  for(int t=0;t<T;t++)
	    {
	      data[icombo][t][ijack]=record->vect[ijack+(njacks+1)*(0+2*(t+T*icombo))];
	      //cout<<record->vect[ijack+(njacks+1)*(0+2*(t+T*icombo))]<<endl;
	    }
      }
    
    //add the sign
    for(int icombo=0;icombo<record->ntot_combo;icombo++) data[icombo]*=coeff;
    
    fclose(fin);
    delete record;
  }
  ~corr_t(){free(data);}

  jvec operator[](hadr_t &h)
  {
    //cout<<h.a()<<" "<<h.b()<<" "<<h.c()<<" "<<h.d()<<" : ";
    
    return (data[h.a()]+data[h.b()]-(data[h.c()]+data[h.d()]).simmetric()).subset(0,TH)/4;
  }
private:
  corr_t() {};
};

//  jvec All,Eel,Ele,Lal,Lee,Lla,Lll,Llp,Lls,Llt,Lpl,Lsl,Ltl,Pll,Sll,Tll;
////////////////////////////////////////////////////////////////////////////////////////////

int main(int narg,char **arg)
{
  read_analysis_pars();
  
  corr_t direct_Lll("data/direct_oC5C5o-Lll_conf.1.dat");
  corr_t exchange_Lll("data/exchange_oC5C5o-Lll_conf.1.dat",-1);
  
  //comparing
  cout<<"nucleon_mass_000: "<<constant_fit(aperiodic_effective_mass(direct_Lll[nuc_d_000]-exchange_Lll[nuc_d_000]),tmin_2pts,tmax_2pts,"plots/comparison_rcombo/nucleon_mass_000.xmg")<<endl;
  cout<<"nucleon_mass_001: "<<constant_fit(aperiodic_effective_mass(direct_Lll[nuc_d_001]-exchange_Lll[nuc_d_001]),tmin_2pts,tmax_2pts,"plots/comparison_rcombo/nucleon_mass_001.xmg")<<endl;
  cout<<"nucleon_mass_010: "<<constant_fit(aperiodic_effective_mass(direct_Lll[nuc_d_010]-exchange_Lll[nuc_d_010]),tmin_2pts,tmax_2pts,"plots/comparison_rcombo/nucleon_mass_010.xmg")<<endl;
  cout<<"nucleon_mass_011: "<<constant_fit(aperiodic_effective_mass(direct_Lll[nuc_d_011]-exchange_Lll[nuc_d_011]),tmin_2pts,tmax_2pts,"plots/comparison_rcombo/nucleon_mass_011.xmg")<<endl;
  cout<<"nucleon_mass_100: "<<constant_fit(aperiodic_effective_mass(direct_Lll[nuc_d_100]-exchange_Lll[nuc_d_100]),tmin_2pts,tmax_2pts,"plots/comparison_rcombo/nucleon_mass_100.xmg")<<endl;
  cout<<"nucleon_mass_101: "<<constant_fit(aperiodic_effective_mass(direct_Lll[nuc_d_101]-exchange_Lll[nuc_d_101]),tmin_2pts,tmax_2pts,"plots/comparison_rcombo/nucleon_mass_101.xmg")<<endl;
  cout<<"nucleon_mass_110: "<<constant_fit(aperiodic_effective_mass(direct_Lll[nuc_d_110]-exchange_Lll[nuc_d_110]),tmin_2pts,tmax_2pts,"plots/comparison_rcombo/nucleon_mass_110.xmg")<<endl;
  cout<<"nucleon_mass_111: "<<constant_fit(aperiodic_effective_mass(direct_Lll[nuc_d_111]-exchange_Lll[nuc_d_111]),tmin_2pts,tmax_2pts,"plots/comparison_rcombo/nucleon_mass_111.xmg")<<endl;
  
  {ofstream naz("nazario_tables/direct_Lll_000");naz<<aperiodic_effective_mass(direct_Lll[nuc_d_000]);}
  {ofstream naz("nazario_tables/direct_Lll_001");naz<<direct_Lll[nuc_d_001];}
  {ofstream naz("nazario_tables/direct_Lll_010");naz<<direct_Lll[nuc_d_010];}
  {ofstream naz("nazario_tables/direct_Lll_011");naz<<direct_Lll[nuc_d_011];}
  {ofstream naz("nazario_tables/direct_Lll_100");naz<<direct_Lll[nuc_d_100];}
  {ofstream naz("nazario_tables/direct_Lll_101");naz<<direct_Lll[nuc_d_101];}
  {ofstream naz("nazario_tables/direct_Lll_110");naz<<direct_Lll[nuc_d_110];}
  {ofstream naz("nazario_tables/direct_Lll_111");naz<<direct_Lll[nuc_d_111];}
  {ofstream naz("nazario_tables/exchange_Lll_000");naz<<exchange_Lll[nuc_d_000];}
  {ofstream naz("nazario_tables/exchange_Lll_001");naz<<exchange_Lll[nuc_d_001];}
  {ofstream naz("nazario_tables/exchange_Lll_010");naz<<exchange_Lll[nuc_d_010];}
  {ofstream naz("nazario_tables/exchange_Lll_011");naz<<exchange_Lll[nuc_d_011];}
  {ofstream naz("nazario_tables/exchange_Lll_100");naz<<exchange_Lll[nuc_d_100];}
  {ofstream naz("nazario_tables/exchange_Lll_101");naz<<exchange_Lll[nuc_d_101];}
  {ofstream naz("nazario_tables/exchange_Lll_110");naz<<exchange_Lll[nuc_d_110];}
  {ofstream naz("nazario_tables/exchange_Lll_111");naz<<exchange_Lll[nuc_d_111];}
  
  separator(); ////////////// hadron mass /////////////
  
  //nucleon
  jvec nucleon=direct_Lll[nuc_d]-exchange_Lll[nuc_d];
  jack nucleon_mass=constant_fit(aperiodic_effective_mass(nucleon),tmin_2pts,tmax_2pts,"plots/nucleon_mass.xmg");
  cout<<"Nucleon mass: "<<smart_print(nucleon_mass)<<" = "<<smart_print(nucleon_mass/a[ib])<<" MeV (Nucleon0 exp: "
      <<nucleon0<<" MeV)"<<endl;
  
  //sigma
  jvec sigma=direct_Lll[sig_d]-exchange_Lll[sig_d];
  jack sigma_mass=constant_fit(aperiodic_effective_mass(sigma),tmin_2pts,tmax_2pts,"plots/sigma_mass.xmg");
  cout<<"Sigma mass: "<<smart_print(sigma_mass)<<" = "<<smart_print(sigma_mass/a[ib])<<" MeV (Sigma0 exp: "
      <<sigma0<<" MeV)"<<endl;
  
  //xi
  jvec xi=direct_Lll[xi_d]-exchange_Lll[xi_d];
  jack xi_mass=constant_fit(aperiodic_effective_mass(xi),tmin_2pts,tmax_2pts,"plots/xi_mass.xmg");
  cout<<"Xi mass: "<<smart_print(xi_mass)<<" = "<<smart_print(xi_mass/a[ib])<<" MeV (Xi0 exp: "
      <<xi0<<" MeV)"<<endl;
  
  separator(); //////////// mass correction ////////////
  
  corr_t direct_Sll("data/direct_oC5C5o-Sll_conf.1.dat");
  corr_t exchange_Sll("data/exchange_oC5C5o-Sll_conf.1.dat",-1);
  
  corr_t direct_Lsl("data/direct_oC5C5o-Lsl_conf.1.dat");
  corr_t exchange_Lsl("data/exchange_oC5C5o-Lsl_conf.1.dat",-1);
  
  corr_t direct_Lls("data/direct_oC5C5o-Lls_conf.1.dat");
  corr_t exchange_Lls("data/exchange_oC5C5o-Lls_conf.1.dat",-1);
  
  jvec nucleon_dm=direct_Lls[nuc_d]-exchange_Lls[nuc_d]+exchange_Lsl[nuc_d]-exchange_Sll[nuc_d];
  jvec sigma_dm=direct_Lls[sig_d]+direct_Sll[sig_d]-exchange_Lls[sig_d]-exchange_Sll[sig_d];
  jvec xi_dm=direct_Lsl[xi_d]-exchange_Lsl[xi_d];
  
  jack mslope_nuc=constant_fit(numerical_derivative(nucleon_dm/nucleon),tmin_mass_slope,tmax_mass_slope,
			       "plots/nucleon_mass_slope.xmg");
  jack mslope_sig=constant_fit(numerical_derivative(sigma_dm/sigma),tmin_mass_slope,tmax_mass_slope,
			       "plots/sigma_mass_slope.xmg");
  jack mslope_xi=constant_fit(numerical_derivative(xi_dm/xi),tmin_mass_slope,tmax_mass_slope,
			      "plots/xi_mass_slope.xmg");
  
  cout<<"Mslope_nuc: "<<smart_print(mslope_nuc)<<endl;
  cout<<"Mslope_sig: "<<smart_print(mslope_sig)<<endl;
  cout<<"Mslope_xi: "<<smart_print(mslope_xi)<<endl;
  
  separator(); //////////// kappa correction ////////////
  
  corr_t direct_Pll("data/direct_oC5C5o-Pll_conf.1.dat");
  corr_t exchange_Pll("data/exchange_oC5C5o-Pll_conf.1.dat",-1);

  corr_t direct_Lpl("data/direct_oC5C5o-Lpl_conf.1.dat");
  corr_t exchange_Lpl("data/exchange_oC5C5o-Lpl_conf.1.dat",-1);
  
  corr_t direct_Llp("data/direct_oC5C5o-Llp_conf.1.dat");
  corr_t exchange_Llp("data/exchange_oC5C5o-Llp_conf.1.dat",-1);
  
  jvec nucleon_dk=direct_Llp[nuc_d]-direct_Lpl[nuc_d]+direct_Pll[nuc_d]+
    -exchange_Llp[nuc_d]+exchange_Lpl[nuc_d]-exchange_Pll[nuc_d];
  jvec sigma_dk=direct_Llp[sig_d]+direct_Pll[sig_d]-exchange_Llp[sig_d]-exchange_Pll[sig_d];
  jvec xi_dk=direct_Lpl[xi_d]-exchange_Lpl[xi_d];
  
  cout<<"Kslope_nuc: "<<smart_print(constant_fit(numerical_derivative(nucleon_dk/nucleon),tmin_mass_slope,tmax_mass_slope,
						 "plots/nucleon_kappa_slope.xmg"))<<endl;
  cout<<"Kslope_sig: "<<smart_print(constant_fit(numerical_derivative(sigma_dk/sigma),tmin_mass_slope,tmax_mass_slope,
						 "plots/sigma_kappa_slope.xmg"))<<endl;
  cout<<"Kslope_xi: "<<smart_print(constant_fit(numerical_derivative(xi_dk/xi),tmin_mass_slope,tmax_mass_slope,
						"plots/xi_kappa_slope.xmg"))<<endl;
  
  separator(); ///////////// elettromagnetic slope //////////////
  
  corr_t direct_Eel("data/direct_oC5C5o-Eel_conf.1.dat");
  corr_t exchange_Eel("data/exchange_oC5C5o-Eel_conf.1.dat",-1);
  
  corr_t direct_Ele("data/direct_oC5C5o-Ele_conf.1.dat");
  corr_t exchange_Ele("data/exchange_oC5C5o-Ele_conf.1.dat",-1);

  corr_t direct_Lee("data/direct_oC5C5o-Lee_conf.1.dat");
  corr_t exchange_Lee("data/exchange_oC5C5o-Lee_conf.1.dat",-1);
  
  corr_t direct_All("data/direct_oC5C5o-All_conf.1.dat");
  corr_t exchange_All("data/exchange_oC5C5o-All_conf.1.dat",-1);
  
  corr_t direct_Lal("data/direct_oC5C5o-Lal_conf.1.dat");
  corr_t exchange_Lal("data/exchange_oC5C5o-Lal_conf.1.dat",-1);
  
  corr_t direct_Lla("data/direct_oC5C5o-Lla_conf.1.dat");
  corr_t exchange_Lla("data/exchange_oC5C5o-Lla_conf.1.dat",-1);
  
  corr_t direct_Tll("data/direct_oC5C5o-Tll_conf.1.dat");
  corr_t exchange_Tll("data/exchange_oC5C5o-Tll_conf.1.dat",-1);
  
  corr_t direct_Ltl("data/direct_oC5C5o-Ltl_conf.1.dat");
  corr_t exchange_Ltl("data/exchange_oC5C5o-Ltl_conf.1.dat",-1);
  
  corr_t direct_Llt("data/direct_oC5C5o-Llt_conf.1.dat");
  corr_t exchange_Llt("data/exchange_oC5C5o-Llt_conf.1.dat",-1);
  
  jvec nucleon_de=direct_Eel[nuc_d]+direct_All[nuc_d]+direct_Tll[nuc_d]-
    (exchange_Lla[nuc_d]+exchange_Llt[nuc_d]+exchange_Ele[nuc_d]+exchange_All[nuc_d]+
     exchange_Tll[nuc_d]-exchange_Lal[nuc_d]-exchange_Ltl[nuc_d]);
  
  jack nucleon_de_slope=constant_fit(numerical_derivative(nucleon_de/nucleon),tmin_elec_slope,
				     tmax_elec_slope,"plots/nucleon_elec_slope.xmg");
  cout<<"Eslope_nuc: "<<smart_print(nucleon_de_slope)<<endl;
  
  jvec sigma_de=direct_Ele[sig_d]-exchange_Ele[sig_d];
  jack sigma_de_slope=constant_fit(-numerical_derivative(sigma_de/sigma),tmin_sigma_elec_slope,tmax_sigma_elec_slope,
				   "plots/sigma_elec_slope.xmg");
  
  cout<<"Eslope_sig: "<<smart_print(sigma_de_slope)<<endl;
  cout<<"Sigma_e_corr: "<<smart_print(sigma_de_slope*sqr(eu-ed)*e2/a[ib])<<" MeV"<<endl;
  cout<<"Expected: "<<(sigmap+sigmam-2*sigma0)<<" MeV"<<endl;
  
  return 0;
}
