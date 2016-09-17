#include "include.h"

double sigmap=1189.37;
double sigmam=1197.449;
double sigma0=1192.642;
double a=0.086/197;

//#define DEBUG
#define LOAD_GIUSTI

const int RE=0,IM=1;
const int EVN=0,ODD=1;
const int NO_SIGN=1,CHANGE_SIGN=-1;
#ifdef DEBUG
const int ext_islsrc=0,ext_islsink=0;
#else
const int ext_islsrc=1,ext_islsink=0;
#endif

//const int ntot_mass=2,nsmear=1;
//const int istr=0;
const int ntot_mass=6,nsmear=2;
const int istr=1;

struct  __attribute__ ((packed)) mega_header_t
{
  int twist,nf,nsrc;
  int T,LX,LY,LZ;
  int nk;
  int nslsrc,nslsink;
  double ksm,beta,ksea,musea,csw;
  double k[ntot_mass],m[ntot_mass];
  int sstepsrc[nsmear];
  int sstepsnk[1];
  int size;
  
  int icorr(int im1,int im2,int im3,int parity,int islsrc=ext_islsrc,int islsink=ext_islsink)
  {return islsrc+nslsrc*(islsink+nslsink*(parity+2*(im1+nk*(im2+nk*im3))));}
};

int T,TH,njacks=15,clust_size,nconfs;
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

#ifdef DEBUG
int ra=0,rb=0,rc=0;
#else
//int ra=0,rb=1,rc=1;
int ra=1,rb=1,rc=1;
#endif
hadr_t nuc_d(ra,rb,rc);
hadr_t sig_d(ra,istr+rb,rc);
hadr_t xi_d(istr+ra,rb,istr+rc);

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
    if(rc!=size) crash("otbained %d instead of %d while reading vect",rc,size);
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
  for(int i=0;i<in.nk;i++) cout<<"K["<<i<<"]: "<<in.k[i]<<", M["<<i<<"]: "<<in.m[i]<<endl;
  for(int i=0;i<in.nslsrc;i++) cout<<"Sstepsrc["<<i<<"]: "<<in.sstepsrc[i]<<endl;
  for(int i=0;i<in.nslsink;i++) cout<<"Sstepsnk["<<i<<"]: "<<in.sstepsnk[i]<<endl;
  
  return out;
}

////////////////////////////////////// holding various contractions ////////////////////////

struct corr_t
{
  jvec *data;
  int rpar;
  corr_t(const char *path,int rpar,int ri,int sign) : rpar(rpar)
  {
    FILE *fin=open_file(path,"r");
    
    //load header
    int rc=fread(&mega_header,sizeof(mega_header_t),1,fin);
    if(rc!=1) crash("obtained %d while reading",rc);
    //cout<<mega_header<<endl;
    //read the record
    record_t *record=new record_t(mega_header);
    data=(jvec*)malloc(record->ntot_combo*sizeof(jvec));
    for(int icombo=0;icombo<record->ntot_combo;icombo++)
      {
	data[icombo].create(T,njacks);
	data[icombo]=0;
      }
    //reading vector
    for(int iconf=0;iconf<nconfs;iconf++)
      {
	record->read_vect(fin);
	int ijack=iconf/clust_size;
	for(int icombo=0;icombo<record->ntot_combo;icombo++)
	  for(int t=0;t<T;t++)
	    data[icombo][t][ijack]+=record->vect[ri+2*(t+T*icombo)]*mega_header.nsrc/(mega_header.nsrc+1);
	//cout.precision(16);
	//cout<<record->vect[ri+2*(0+T*0)]<<endl;
      }
    //clusterize and put sign
    for(int icombo=0;icombo<record->ntot_combo;icombo++)
      {
	data[icombo].clusterize(clust_size);
	data[icombo]*=sign;
      }
    //cout<<" "<<data[0][0]<<endl;
    
    fclose(fin);
    delete record;
  }
  ~corr_t(){free(data);}
  
  jvec operator[](hadr_t &h)
  {
    //cout<<h.a()<<" "<<h.b()<<" "<<h.c()<<" "<<h.d()<<endl;
#ifdef DEBUG
    return data[h.a()];}
#else
  return (data[h.a()]+(1-2*rpar)*data[h.b()]-(data[h.c()]+(1-2*rpar)*data[h.d()]).simmetric()).subset(0,TH+1)/4;}
#endif
  
private:
  corr_t() {};
};

//  jvec All,Eel,Ele,Lal,Lee,Lla,Lll,Llp,Lls,Llt,Lpl,Lsl,Ltl,Pll,Sll,Tll;
////////////////////////////////////////////////////////////////////////////////////////////

int main(int narg,char **arg)
{
  FILE *input=open_file("input","r");
  read_formatted_from_file_expecting((char*)(&T),input,"%d","T");
  TH=T/2;
  read_formatted_from_file_expecting((char*)(&nconfs),input,"%d","NConfs");
  clust_size=nconfs/njacks;
  nconfs=clust_size*njacks;
  
  cout<<"Size of the header: "<<sizeof(mega_header_t)<<endl;
  cout<<"Combination of r: "<<ra<<" "<<rb<<" "<<rc<<endl;
  cout<<"Source: "<<(ext_islsrc?"SME":"LOC")<<endl;
#ifdef DEBUG
  cout<<"WARNING avoided loading the four combos"<<endl;
#endif
  
  corr_t direct_Lll("data/direct_oC5C5o-Lll_conf.1.dat",EVN,RE,NO_SIGN);
  direct_Lll.data[nuc_d.a()].print_to_file("plots/direct_Lll_a.xmg");
  direct_Lll.data[nuc_d.b()].print_to_file("plots/direct_Lll_b.xmg");
  direct_Lll.data[nuc_d.c()].print_to_file("plots/direct_Lll_c.xmg");
  direct_Lll.data[nuc_d.d()].print_to_file("plots/direct_Lll_d.xmg");
  corr_t exchange_Lll("data/exchange_oC5C5o-Lll_conf.1.dat",EVN,RE,CHANGE_SIGN);
  
  cout<<mega_header<<endl;
  
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
  jack nucleon_mass=constant_fit(aperiodic_effective_mass(nucleon),10,15,"plots/nucleon_mass.xmg");
  cout<<"Nucleon mass: "<<smart_print(nucleon_mass)<<" = "<<smart_print(nucleon_mass/a)<<" MeV"<<endl;
  
  //sigma
  jvec sigma=direct_Lll[sig_d]-exchange_Lll[sig_d];
  jack sigma_mass=constant_fit(aperiodic_effective_mass(sigma),11,15,"plots/sigma_mass.xmg");
  cout<<"Sigma mass: "<<smart_print(sigma_mass)<<" = "<<smart_print(sigma_mass/a)<<" MeV"<<endl;
  
  //xi
  jvec xi=direct_Lll[xi_d]-exchange_Lll[xi_d];
  jack xi_mass=constant_fit(aperiodic_effective_mass(xi),11,15,"plots/xi_mass.xmg");
  cout<<"Xi mass: "<<smart_print(xi_mass)<<" = "<<smart_print(xi_mass/a)<<" MeV"<<endl;
  
  //S has a minus to be added
  
  corr_t direct_Sll("data/direct_oC5C5o-Sll_conf.1.dat",EVN,RE,NO_SIGN          *CHANGE_SIGN);
  corr_t exchange_Sll("data/exchange_oC5C5o-Sll_conf.1.dat",EVN,RE,CHANGE_SIGN  *CHANGE_SIGN);
  
  corr_t direct_Lsl("data/direct_oC5C5o-Lsl_conf.1.dat",EVN,RE,NO_SIGN          *CHANGE_SIGN);
  corr_t exchange_Lsl("data/exchange_oC5C5o-Lsl_conf.1.dat",EVN,RE,CHANGE_SIGN  *CHANGE_SIGN);
  
  corr_t direct_Lls("data/direct_oC5C5o-Lls_conf.1.dat",EVN,RE,NO_SIGN          *CHANGE_SIGN);
  corr_t exchange_Lls("data/exchange_oC5C5o-Lls_conf.1.dat",EVN,RE,CHANGE_SIGN  *CHANGE_SIGN);

  {
  corr_t direct_Llp("data/direct_oC5C5o-Llp_conf.1.dat",ODD,IM,NO_SIGN);
  direct_Llp.data[nuc_d.a()].print_to_file("plots/direct_Llp_a.xmg");
  direct_Llp.data[nuc_d.b()].print_to_file("plots/direct_Llp_b.xmg");
  direct_Llp.data[nuc_d.c()].print_to_file("plots/direct_Llp_c.xmg");
  direct_Llp.data[nuc_d.d()].print_to_file("plots/direct_Llp_d.xmg");
  }
#ifdef LOAD_GIUSTI
  const int nins=16;
  char ins[nins][4]={"Lll","Lls","Lsl","Sll","Llp","Lpl","Pll","Tll","Ltl","Llt","Eel","Ele","Lee","All","Lal","Lla"};
  int rpar[nins]=   {    0,    0,    0,    0,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0};
  int ri[nins]=     {    0,    0,    0,    0,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0};
  int si[nins]=     {   +1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   +1,   +1,   +1,   +1,   +1,   +1};
  char direxch[2][15]={"direct","exchange"};
  int sign[2]={NO_SIGN,CHANGE_SIGN};
  for(int iexch=0;iexch<2;iexch++)
    for(int i=0;i<nins/*only 4 */;i++)
      {
	cout<<" giustifying "<<ins[i]<<" "<<direxch[iexch]<<endl;
	corr_t corr(combine("data/%s_oC5C5o-%s_conf.1.dat",direxch[iexch],ins[i]).c_str(),rpar[i],ri[i],sign[iexch]*si[i]);
	ofstream out(combine("ascii/%s_%s.txt",direxch[iexch],ins[i]));
	if(!out.good()) crash("check ascii folder");
	out.precision(16);
	jvec c=corr[nuc_d];
	c.print_to_file("/tmp/%s_%s.xmg",ins[i],direxch[iexch]);
	for(int t=0;t<=TH;t++)
	  for(int ijack=0;ijack<=njacks;ijack++)
	    out<<t<<" "<<ijack<<" "<<c[t][ijack]<<endl;
      }
  
  return 0;
#endif
  
  jvec nucleon_dm
    =direct_Sll[nuc_d]-direct_Lsl[nuc_d]+direct_Lls[nuc_d]-
    (exchange_Sll[nuc_d]-exchange_Lsl[nuc_d]+exchange_Lls[nuc_d]);
  jvec sigma_dm=direct_Lls[sig_d]+direct_Sll[sig_d]-exchange_Lls[sig_d]-exchange_Sll[sig_d]; //tofix
  jvec xi_dm=direct_Lsl[xi_d]-exchange_Lsl[xi_d]; //idem
  
  jvec nucleon_dm_use_symm
    =direct_Lls[nuc_d]-
    (exchange_Sll[nuc_d]-exchange_Lsl[nuc_d]+exchange_Lls[nuc_d]);
  
  direct_Lls[nuc_d].print_to_file("plots/direct_Lls.xmg");
  direct_Lsl[nuc_d].print_to_file("plots/direct_Lsl.xmg");
  direct_Sll[nuc_d].print_to_file("plots/direct_Sll.xmg");
  exchange_Lls[nuc_d].print_to_file("plots/exchange_Lls.xmg");
  exchange_Lsl[nuc_d].print_to_file("plots/exchange_Lsl.xmg");
  exchange_Sll[nuc_d].print_to_file("plots/exchange_Sll.xmg");
  (direct_Lls[nuc_d]-exchange_Lls[nuc_d]).print_to_file("plots/Lls.xmg");
  (direct_Lsl[nuc_d]-exchange_Lsl[nuc_d]).print_to_file("plots/Lsl.xmg");
  (direct_Sll[nuc_d]-exchange_Sll[nuc_d]).print_to_file("plots/Sll.xmg");
  
  corr_t direct_All("data/direct_oC5C5o-All_conf.1.dat",EVN,RE,NO_SIGN);
  corr_t exchange_All("data/exchange_oC5C5o-All_conf.1.dat",EVN,RE,CHANGE_SIGN);
  corr_t direct_Lal("data/direct_oC5C5o-Lal_conf.1.dat",EVN,RE,NO_SIGN);
  corr_t exchange_Lal("data/exchange_oC5C5o-Lal_conf.1.dat",EVN,RE,CHANGE_SIGN);
  corr_t direct_Lla("data/direct_oC5C5o-Lla_conf.1.dat",EVN,RE,NO_SIGN);
  direct_Lla.data[nuc_d.a()].print_to_file("plots/direct_Lla_a.xmg");
  direct_Lla.data[nuc_d.b()].print_to_file("plots/direct_Lla_b.xmg");
  direct_Lla.data[nuc_d.c()].print_to_file("plots/direct_Lla_c.xmg");
  direct_Lla.data[nuc_d.d()].print_to_file("plots/direct_Lla_d.xmg");
  corr_t exchange_Lla("data/exchange_oC5C5o-Lla_conf.1.dat",EVN,RE,CHANGE_SIGN);
  direct_Lla[nuc_d].print_to_file("plots/direct_Lla.xmg");
  direct_Lal[nuc_d].print_to_file("plots/direct_Lal.xmg");
  direct_All[nuc_d].print_to_file("plots/direct_All.xmg");
  exchange_Lla[nuc_d].print_to_file("plots/exchange_Lla.xmg");
  exchange_Lal[nuc_d].print_to_file("plots/exchange_Lal.xmg");
  exchange_All[nuc_d].print_to_file("plots/exchange_All.xmg");
  (direct_Lla[nuc_d]-exchange_Lla[nuc_d]).print_to_file("plots/Lla.xmg");
  (direct_Lal[nuc_d]-exchange_Lal[nuc_d]).print_to_file("plots/Lal.xmg");
  (direct_All[nuc_d]-exchange_All[nuc_d]).print_to_file("plots/All.xmg");
  direct_All.data[nuc_d.a()].print_to_file("plots/direct_All_a.xmg");
  direct_All.data[nuc_d.b()].print_to_file("plots/direct_All_b.xmg");
  direct_All.data[nuc_d.c()].print_to_file("plots/direct_All_c.xmg");
  direct_All.data[nuc_d.d()].print_to_file("plots/direct_All_d.xmg");

  ((direct_Lal[nuc_d]-exchange_Lal[nuc_d])/(direct_Lll[nuc_d]-exchange_Lll[nuc_d])).print_to_file("/tmp/Lal_slope.xmg");

  corr_t direct_Tll("data/direct_oC5C5o-Tll_conf.1.dat",EVN,RE,NO_SIGN);
  corr_t exchange_Tll("data/exchange_oC5C5o-Tll_conf.1.dat",EVN,RE,CHANGE_SIGN);
  corr_t direct_Ltl("data/direct_oC5C5o-Ltl_conf.1.dat",EVN,RE,NO_SIGN);
  corr_t exchange_Ltl("data/exchange_oC5C5o-Ltl_conf.1.dat",EVN,RE,CHANGE_SIGN);
  corr_t direct_Llt("data/direct_oC5C5o-Llt_conf.1.dat",EVN,RE,NO_SIGN);
  corr_t exchange_Llt("data/exchange_oC5C5o-Llt_conf.1.dat",EVN,RE,CHANGE_SIGN);
  direct_Llt[nuc_d].print_to_file("plots/direct_Llt.xmg");
  direct_Ltl[nuc_d].print_to_file("plots/direct_Ltl.xmg");
  direct_Tll[nuc_d].print_to_file("plots/direct_Tll.xmg");
  exchange_Llt[nuc_d].print_to_file("plots/exchange_Llt.xmg");
  exchange_Ltl[nuc_d].print_to_file("plots/exchange_Ltl.xmg");
  exchange_Tll[nuc_d].print_to_file("plots/exchange_Tll.xmg");
  (direct_Llt[nuc_d]-exchange_Llt[nuc_d]).print_to_file("plots/Llt.xmg");
  (direct_Ltl[nuc_d]-exchange_Ltl[nuc_d]).print_to_file("plots/Ltl.xmg");
  (direct_Tll[nuc_d]-exchange_Tll[nuc_d]).print_to_file("plots/Tll.xmg");
  
  corr_t direct_Lee("data/direct_oC5C5o-Lee_conf.1.dat",EVN,RE,NO_SIGN);
  corr_t exchange_Lee("data/exchange_oC5C5o-Lee_conf.1.dat",EVN,RE,CHANGE_SIGN);
  corr_t direct_Ele("data/direct_oC5C5o-Ele_conf.1.dat",EVN,RE,NO_SIGN);
  corr_t exchange_Ele("data/exchange_oC5C5o-Ele_conf.1.dat",EVN,RE,CHANGE_SIGN);
  corr_t direct_Eel("data/direct_oC5C5o-Eel_conf.1.dat",EVN,RE,NO_SIGN);
  corr_t exchange_Eel("data/exchange_oC5C5o-Eel_conf.1.dat",EVN,RE,CHANGE_SIGN);
  direct_Eel[nuc_d].print_to_file("plots/direct_Eel.xmg");
  direct_Ele[nuc_d].print_to_file("plots/direct_Ele.xmg");
  direct_Lee[nuc_d].print_to_file("plots/direct_Lee.xmg");
  exchange_Eel[nuc_d].print_to_file("plots/exchange_Eel.xmg");
  exchange_Ele[nuc_d].print_to_file("plots/exchange_Ele.xmg");
  exchange_Lee[nuc_d].print_to_file("plots/exchange_Lee.xmg");
  (direct_Eel[nuc_d]-exchange_Eel[nuc_d]).print_to_file("plots/Eel.xmg");
  (direct_Ele[nuc_d]-exchange_Ele[nuc_d]).print_to_file("plots/Ele.xmg");
  (direct_Lee[nuc_d]-exchange_Lee[nuc_d]).print_to_file("plots/Lee.xmg");
  
  corr_t direct_Pll("data/direct_oC5C5o-Pll_conf.1.dat",ODD,IM,NO_SIGN);
  corr_t exchange_Pll("data/exchange_oC5C5o-Pll_conf.1.dat",ODD,IM,CHANGE_SIGN);
  corr_t direct_Lpl("data/direct_oC5C5o-Lpl_conf.1.dat",ODD,IM,NO_SIGN);
  corr_t exchange_Lpl("data/exchange_oC5C5o-Lpl_conf.1.dat",ODD,IM,CHANGE_SIGN);
  corr_t direct_Llp("data/direct_oC5C5o-Llp_conf.1.dat",EVN,IM,NO_SIGN);
  corr_t exchange_Llp("data/exchange_oC5C5o-Llp_conf.1.dat",ODD,IM,CHANGE_SIGN);
  direct_Llp[nuc_d].print_to_file("plots/direct_Llp.xmg");
  direct_Lpl[nuc_d].print_to_file("plots/direct_Lpl.xmg");
  direct_Pll[nuc_d].print_to_file("plots/direct_Pll.xmg");
  exchange_Llp[nuc_d].print_to_file("plots/exchange_Llp.xmg");
  exchange_Lpl[nuc_d].print_to_file("plots/exchange_Lpl.xmg");
  exchange_Pll[nuc_d].print_to_file("plots/exchange_Pll.xmg");
  (direct_Llp[nuc_d]-exchange_Llp[nuc_d]).print_to_file("plots/Llp.xmg");
  (direct_Lpl[nuc_d]-exchange_Lpl[nuc_d]).print_to_file("plots/Lpl.xmg");
  (direct_Pll[nuc_d]-exchange_Pll[nuc_d]).print_to_file("plots/Pll.xmg");
  
  jvec nucleon_dp=direct_Llp[nuc_d]-exchange_Llp[nuc_d];
  nucleon.print_to_file("plots/nucleon.xmg");
  direct_Lll[nuc_d].print_to_file("plots/direct_Lll.xmg");
  exchange_Lll[nuc_d].print_to_file("plots/exchange_Lll.xmg");
  nucleon_dm.print_to_file("plots/nucleon_dm.xmg");
  nucleon_dp.print_to_file("plots/nucleon_dp.xmg");
  
  (nucleon_dm/nucleon).print_to_file("plots/nucleon_mass_slope_corr.xmg");
  cout<<"Mslope_nuc: "<<constant_fit(numerical_derivative(nucleon_dm/nucleon),6,15,"plots/nucleon_mass_slope.xmg")<<endl;
  cout<<"Mslope_nuc_use_symm: "<<constant_fit(numerical_derivative(nucleon_dm_use_symm/nucleon),6,15,"plots/nucleon_mass_slope_use_symm.xmg")<<endl;
  cout<<"Mslope_sig: "<<constant_fit(numerical_derivative(sigma_dm/sigma),6,15,"plots/sigma_mass_slope.xmg")<<endl;
  cout<<"Mslope_xi: "<<constant_fit(numerical_derivative(xi_dm/xi),6,15,"plots/xi_mass_slope.xmg")<<endl;
  
  ///////////// elettromagnetic slope //////////////
    
  jvec sigma_de=direct_Ele[sig_d]-exchange_Ele[sig_d];
  jack sigma_de_slope=constant_fit(-numerical_derivative(sigma_de/sigma),3,9,"plots/sigma_elec_slope.xmg");
  cout<<"Eslope_sig: "<<smart_print(sigma_de_slope)<<endl;
  cout<<"Sigma_e_corr: "<<smart_print(sigma_de_slope*sqr(eu-ed)*e2/a)<<" MeV"<<endl;
  cout<<"Expected: "<<(sigmap+sigmam-2*sigma0)<<" MeV"<<endl;
  
  return 0;
}
