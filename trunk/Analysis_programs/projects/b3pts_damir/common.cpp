#include "include.h"

//0.713 taken from 1201.5025
const double Zp_fr_Zs[4]={0.580,0.626,0.686,713};
const double Zv[4]={0.5816,0.6103,0.6451,0.686};
const double Zt[4]={0.73,0.74,0.78,0.82};

enum PARITY{EVN,ODD,UNC};
enum REIM{RE,IM};

enum S_2pts_mass{LIGHT_S_2PTS,STRANGE_S_2PTS,H0_S_2PTS,H1_S_2PTS,H2_S_2PTS,H3_S_2PTS,H4_S_2PTS,H5_S_2PTS,H6_S_2PTS,H7_S_2PTS,H8_S_2PTS,H9_S_2PTS};
S_2pts_mass CHARM_S_2PTS=H0_S_2PTS;
S_2pts_mass L_S_2PTS[3]={LIGHT_S_2PTS,STRANGE_S_2PTS,H0_S_2PTS};
S_2pts_mass H_S_2PTS[10]={H0_S_2PTS,H1_S_2PTS,H2_S_2PTS,H3_S_2PTS,H4_S_2PTS,H5_S_2PTS,H6_S_2PTS,H7_S_2PTS,H8_S_2PTS,H9_S_2PTS};

enum Spec_mass{LIGHT_SPEC,STRANGE_SPEC};
Spec_mass SPEC[2]={LIGHT_SPEC,STRANGE_SPEC};

enum S0_3pts_mass{LIGHT_S0_3PTS,STRANGE_S0_3PTS,CHARM_S0_3PTS};
S0_3pts_mass L_S0_3PTS[3]={LIGHT_S0_3PTS,STRANGE_S0_3PTS,CHARM_S0_3PTS};

enum S1_3pts_mass{H0_S1_3PTS,H1_S1_3PTS,H2_S1_3PTS,H3_S1_3PTS,H4_S1_3PTS,H5_S1_3PTS,H6_S1_3PTS,H7_S1_3PTS,H8_S1_3PTS,H9_S1_3PTS};
S1_3pts_mass CHARM_S1_3PTS=H0_S1_3PTS;
S1_3pts_mass H_S1_3PTS[10]={H0_S1_3PTS,H1_S1_3PTS,H2_S1_3PTS,H3_S1_3PTS,H4_S1_3PTS,H5_S1_3PTS,H6_S1_3PTS,H7_S1_3PTS,H8_S1_3PTS,H9_S1_3PTS};

int T,TH,tsep;
const int njack=16;

int ibeta;
int nh=10,nl=2,nm=nh+nl,nm_S0=3;
double *mass;

int nth_S0;
double *th_S0;

char ori_meson[2][3][15]={{"Pi","K(lspec)","D"},{"K(sspec)","eta","Ds"}};
char prod_meson[2][2][15]={{"D","H"},{"Ds","Hs"}};

struct info_3pts
{
  char name[5];
  PARITY pa;
  REIM ri;
public:
  info_3pts(const char *ext_name,PARITY pa,REIM ri) : pa(pa),ri(ri) {strcpy(name,ext_name);};
private:
  info_3pts();
};

int map_pa[2]={1,-1};

info_3pts S0P5("S0P5",EVN,RE);
info_3pts VKP5("VKP5",EVN,IM);
info_3pts V0P5("V0P5",ODD,RE);
info_3pts TKP5("TKP5",ODD,IM);
info_3pts VJVK("VJVK",ODD,IM);
info_3pts VKVJ("VKVJ",ODD,IM);
info_3pts P5VK("P5VK",EVN,IM);
info_3pts AKVK("AKVK",EVN,RE);
info_3pts AJVK("AJVK",EVN,RE); //Horr
info_3pts AKVJ("AKVJ",EVN,RE); //Horr
info_3pts A0VK("A0VK",ODD,IM);
info_3pts AKV0("AKV0",ODD,IM);
info_3pts A0V0("A0V0",ODD,RE); //Horr
info_3pts TJVK("TJVK",ODD,IM); //Horr
info_3pts TKVJ("TKVJ",ODD,IM); //Horr
info_3pts BKVK("BKVK",ODD,RE);
info_3pts BJVK("BJVK",ODD,RE);
info_3pts BKVJ("BKVJ",ODD,RE);
info_3pts BKV0("BKV0",EVN,IM);

int n_ME=19;
info_3pts info_ME[19]={S0P5,VKP5,V0P5,TKP5,VJVK,VKVJ,P5VK,AKVK,AJVK,AKVJ,A0VK,AKV0,A0V0,TJVK,TKVJ,BKVK,BJVK,BKVJ,BKV0};

struct combo_3pts
{
  Spec_mass im_spec;
  S0_3pts_mass im_S0;
  S1_3pts_mass im_S1;
  int ith_S0;
public:
  combo_3pts(Spec_mass im_spec,S0_3pts_mass im_S0,S1_3pts_mass im_S1,int ith_S0) : im_spec(im_spec),im_S0(im_S0),im_S1(im_S1),ith_S0(ith_S0) {};
private:
  combo_3pts();
};

double m_S0(int im)
{
  if(im<0||im>=nm_S0) crash("S0 mass asked %d outside the available interval",im);
  return mass[im];
}

double m_S1(int im)
{
  if(im<0||im>=nh) crash("S1 mass asked %d outside the available interval",im);
  return mass[im+2];
}

double m_spec(int im)
{
  if(im<0||im>=2) crash("spec mass asked %d outside the available interval",im);
  return mass[im];
}

double momentum(double th)
{return M_PI/TH*th;}

//relativistic lattice energy: lattice...
jack latt_en(jack M,double th)
{return 2*asinh(sqrt(3*sqr(sin(momentum(th)/2))+sqr(sinh(M/2))));}
//...and continuum
jack cont_en(jack M,double th)
{return sqrt(M*M+3*sqr(momentum(th)));}

void compute_momentums(jack &P2,jack &P0,double &PK,jack &Q2,jack &Q0,double &QK,jack &m,double th,jack &M)
{
  jack e=cont_en(m,th);
  
  P0=M+e;
  Q0=M-e;
  
  //we are describing the process B->D
  double PB=0;
  double PD=momentum(th);
  
  PK=PB+PD;
  QK=PB-PD;

  P2=P0*P0-3*PK*PK;
  Q2=Q0*Q0-3*QK*QK;
}

void read_set_pars(const char *path)
{
  FILE *set_pars_file=open_file(path,"r");
  
  read_formatted_from_file_expecting((char*)&T,set_pars_file,"%d","T");
  TH=tsep=T/2;
  
  //beta
  read_formatted_from_file_expecting((char*)&ibeta,set_pars_file,"%d","ibeta");
  
  //read the masses
  int nmass; //not really used, only read for consistency
  read_formatted_from_file_expecting((char*)&nmass,set_pars_file,"%d","nmass");
  if(nmass!=nm) crash("recheck nm");
  expect_string_from_file(set_pars_file,"mass_list");
  mass=(double*)malloc(sizeof(double)*nmass);
  for(int imass=0;imass<nmass;imass++) read_formatted_from_file((char*)&(mass[imass]),set_pars_file,"%lg","mass");
  
  //read theta
  read_formatted_from_file_expecting((char*)&nth_S0,set_pars_file,"%d","ntheta");
  expect_string_from_file(set_pars_file,"theta_list");
  th_S0=(double*)malloc(sizeof(double)*nth_S0);
  for(int ith_S0=0;ith_S0<nth_S0;ith_S0++) read_formatted_from_file((char*)&(th_S0[ith_S0]),set_pars_file,"%lg","theta");
  
  fclose(set_pars_file);
}

