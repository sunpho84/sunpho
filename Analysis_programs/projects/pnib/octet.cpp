#include "include.h"

double sigmap=1189.37;
double sigmam=1197.449;
double sigma0=1192.642;
double a=0.086/197;

struct mega_header_t
{
  int twist,nf,nsrc;
  int T,LX,LY,LZ;
  int nk;
  int nslsrc,nslsink;
  double ksm,beta,ksea,musea,csw;
  double k[8];
  int sstepsrc[3];
  static const int ext_islsrc=1,ext_islsink=0;
  int icorr(int im1,int im2,int im3,int parity,int islsrc=ext_islsrc,int islsink=ext_islsink)
  {return islsrc+nslsrc*(islsink+nslsink*(parity+2*(im1+nk*(im2+nk*im3))));}
};

int T=48,TH=24,njacks=26,clust_size=4,nconfs=njacks*clust_size;
mega_header_t mega_header;
double eu=2.0/3,ed=-1.0/3;
double e2=4*M_PI/137;

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

hadr_t nuc_d(0,1,0);
hadr_t sig_d(0,3,0);
hadr_t xi_d(2,1,2);

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
    size=ntot_combo*head.T*2;
    vect=new double[size];
    //cout<<"Ntot_combo: "<<ntot_combo<<endl;
    //cout<<"Size: "<<size<<endl;
  }
  
  //reading vect
  void read_vect(FILE *fin)
  {
    //read iconf
    int iconf;
    int rc=fread(&iconf,sizeof(int),1,fin);
    if(rc!=1) crash("otbained %d while reading iconf",rc);
    //printf("Conf: %d\n",iconf);
    
    //read the vector properly
    rc=fread(vect,sizeof(double),size,fin);
    if(rc!=size) crash("otbained %d while reading vect",rc);
    //for(int i=0;i<size;i++) cout<<i<<" "<<vect[i]<<endl;
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
    if(rc!=1) crash("obtained %d while reading",rc);
    //cout<<mega_header<<endl;
    
    //read the record
    record_t *record=new record_t(mega_header);
    data=(jvec*)malloc(record->ntot_combo*sizeof(jvec));
    for(int icombo=0;icombo<record->ntot_combo;icombo++) data[icombo].create(T,njacks);
    //reading vector
    for(int iconf=0;iconf<nconfs;iconf++)
      {
	record->read_vect(fin);
	int ijack=iconf/clust_size;
	for(int icombo=0;icombo<record->ntot_combo;icombo++)
	  for(int t=0;t<T;t++)
	    data[icombo][t][ijack]+=record->vect[0+2*(t+T*icombo)]/clust_size;
      }
    //clusterize and add the sign
    for(int icombo=0;icombo<record->ntot_combo;icombo++)
      {
	data[icombo].clusterize();
	data[icombo]*=coeff;
      }
    
    fclose(fin);
    delete record;
  }
  ~corr_t(){free(data);}

  jvec operator[](hadr_t &h)
  {cout<<h.a()<<" "<<h.b()<<" "<<h.c()<<" "<<h.d()<<" : ";
    //return data[h.a()];}
    //return data[h.a()]-data[h.c()].simmetric();}
    return (data[h.a()]+data[h.b()]-(data[h.c()]+data[h.d()]).simmetric()).subset(0,TH)/4;}
private:
  corr_t() {};
};

//  jvec All,Eel,Ele,Lal,Lee,Lla,Lll,Llp,Lls,Llt,Lpl,Lsl,Ltl,Pll,Sll,Tll;
////////////////////////////////////////////////////////////////////////////////////////////

int main(int narg,char **arg)
{
  corr_t direct_Lll("data/direct_oC5C5o-Lll_conf.1.dat");
  corr_t exchange_Lll("data/exchange_oC5C5o-Lll_conf.1.dat",-1);
  
  //comparing
  cout<<"nucleon_mass_000: "<<constant_fit(aperiodic_effective_mass(direct_Lll[nuc_d_000]-exchange_Lll[nuc_d_000]),11,15,"plots/comparison_rcombo/nucleon_mass_000.xmg")<<endl;
  cout<<"nucleon_mass_001: "<<constant_fit(aperiodic_effective_mass(direct_Lll[nuc_d_001]-exchange_Lll[nuc_d_001]),11,15,"plots/comparison_rcombo/nucleon_mass_001.xmg")<<endl;
  cout<<"nucleon_mass_010: "<<constant_fit(aperiodic_effective_mass(direct_Lll[nuc_d_010]-exchange_Lll[nuc_d_010]),11,15,"plots/comparison_rcombo/nucleon_mass_010.xmg")<<endl;
  cout<<"nucleon_mass_011: "<<constant_fit(aperiodic_effective_mass(direct_Lll[nuc_d_011]-exchange_Lll[nuc_d_011]),11,15,"plots/comparison_rcombo/nucleon_mass_011.xmg")<<endl;
  cout<<"nucleon_mass_100: "<<constant_fit(aperiodic_effective_mass(direct_Lll[nuc_d_100]-exchange_Lll[nuc_d_100]),11,15,"plots/comparison_rcombo/nucleon_mass_100.xmg")<<endl;
  cout<<"nucleon_mass_101: "<<constant_fit(aperiodic_effective_mass(direct_Lll[nuc_d_101]-exchange_Lll[nuc_d_101]),11,15,"plots/comparison_rcombo/nucleon_mass_101.xmg")<<endl;
  cout<<"nucleon_mass_110: "<<constant_fit(aperiodic_effective_mass(direct_Lll[nuc_d_110]-exchange_Lll[nuc_d_110]),11,15,"plots/comparison_rcombo/nucleon_mass_110.xmg")<<endl;
  cout<<"nucleon_mass_111: "<<constant_fit(aperiodic_effective_mass(direct_Lll[nuc_d_111]-exchange_Lll[nuc_d_111]),11,15,"plots/comparison_rcombo/nucleon_mass_111.xmg")<<endl;
  
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
  
  //nucleon
  jvec nucleon=direct_Lll[nuc_d]-exchange_Lll[nuc_d];
  jack nucleon_mass=constant_fit(aperiodic_effective_mass(nucleon),11,15,"plots/nucleon_mass.xmg");
  cout<<"Nucleon mass: "<<smart_print(nucleon_mass)<<" = "<<smart_print(nucleon_mass/a)<<" MeV"<<endl;
  
  //sigma
  jvec sigma=direct_Lll[sig_d]-exchange_Lll[sig_d];
  jack sigma_mass=constant_fit(aperiodic_effective_mass(sigma),11,15,"plots/sigma_mass.xmg");
  cout<<"Sigma mass: "<<smart_print(sigma_mass)<<" = "<<smart_print(sigma_mass/a)<<" MeV"<<endl;
  
  //xi
  jvec xi=direct_Lll[xi_d]-exchange_Lll[xi_d];
  jack xi_mass=constant_fit(aperiodic_effective_mass(xi),11,15,"plots/xi_mass.xmg");
  cout<<"Xi mass: "<<smart_print(xi_mass)<<" = "<<smart_print(xi_mass/a)<<" MeV"<<endl;  
  
  corr_t direct_Sll("data/direct_oC5C5o-Sll_conf.1.dat");
  corr_t exchange_Sll("data/exchange_oC5C5o-Sll_conf.1.dat",-1);
  
  corr_t direct_Lsl("data/direct_oC5C5o-Lsl_conf.1.dat");
  corr_t exchange_Lsl("data/exchange_oC5C5o-Lsl_conf.1.dat",-1);
  
  corr_t direct_Lls("data/direct_oC5C5o-Lls_conf.1.dat");
  corr_t exchange_Lls("data/exchange_oC5C5o-Lls_conf.1.dat",-1);
  
  jvec nucleon_dm=direct_Lls[nuc_d]-exchange_Lls[nuc_d]+exchange_Lsl[nuc_d]-exchange_Sll[nuc_d];
  jvec sigma_dm=direct_Lls[sig_d]+direct_Sll[sig_d]-exchange_Lls[sig_d]-exchange_Sll[sig_d];
  jvec xi_dm=direct_Lsl[xi_d]-exchange_Lsl[xi_d];
  
  cout<<"Mslope_nuc: "<<constant_fit(numerical_derivative(nucleon_dm/nucleon),6,15,"plots/nucleon_mass_slope.xmg")<<endl;
  cout<<"Mslope_sig: "<<constant_fit(numerical_derivative(sigma_dm/sigma),6,15,"plots/sigma_mass_slope.xmg")<<endl;
  cout<<"Mslope_xi: "<<constant_fit(numerical_derivative(xi_dm/xi),6,15,"plots/xi_mass_slope.xmg")<<endl;
  
  ///////////// elettromagnetic slope //////////////
  
  corr_t direct_Ele("data/direct_oC5C5o-Ele_conf.1.dat");
  corr_t exchange_Ele("data/exchange_oC5C5o-Ele_conf.1.dat",-1);
  
  jvec sigma_de=direct_Ele[sig_d]-exchange_Ele[sig_d];
  jack sigma_de_slope=constant_fit(-numerical_derivative(sigma_de/sigma),3,9,"plots/sigma_elec_slope.xmg");
  cout<<"Eslope_sig: "<<smart_print(sigma_de_slope)<<endl;
  cout<<"Sigma_e_corr: "<<smart_print(sigma_de_slope*sqr(eu-ed)*e2/a)<<" MeV"<<endl;
  cout<<"Expected: "<<(sigmap+sigmam-2*sigma0)<<" MeV"<<endl;
  
  return 0;
}
