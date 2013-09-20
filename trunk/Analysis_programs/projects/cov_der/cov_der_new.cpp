#include "include.h"
#include <map>

char *base_path;
int T=48,TH=24;
int njacks=16;
int ngam=11;
map<int,int> gam_map;
int nsme_lev=2,SME=1;

const int SPIN0=0,SPIN1=1;
const int MAX_COMPONENTS=6;

//indices of operators
struct indices_t
{
  int n;
  int pol[MAX_COMPONENTS];
  int gam[MAX_COMPONENTS];
  int der1[MAX_COMPONENTS];
  int der2[MAX_COMPONENTS];
  double sign[MAX_COMPONENTS];
  void reset() {n=0;}
  indices_t() {reset();}
  void add(int p,int g,int d1,int d2,double s)
  {
    if(n>=MAX_COMPONENTS) crash("asked to define an index holder with %d components",n);
    pol[n]=p;
    gam[n]=g;
    der1[n]=d1;
    der2[n]=d2;
    sign[n]=s;
    n++;
  }
};

//struct containing operator: name (used to load), indices, smearing and anti-hermitianity
struct oper_t
{
  indices_t *ind;
  char name[4];
  int antiherm;
  int spin;
  int sme;
  oper_t(int ext_spin,indices_t &ext_ind,const char *ext_name,int ext_sme,int ext_antiherm) : spin(ext_spin), sme(ext_sme)
  {
    ind=&ext_ind;
    strcpy(name,ext_name);
    antiherm=ext_antiherm;
  }
};

int icombo(int ig_so,int ig_si,int ism_lev_so,int imu1_so,int imu2_so,int ism_lev_si,int imu1_si,int imu2_si)
{
  int ider=imu2_si+4*(imu1_si+4*(ism_lev_si+nsme_lev*(imu2_so+4*(imu1_so+4*ism_lev_so))));
  int ig=gam_map[ig_si]+ngam*gam_map[ig_so];
  
  return ig+ngam*ngam*ider;
}

//compute the average number of standard deviation
double ndev(jvec &temp)
{
  const int dt=19;
  double n=0;
  for(int t=dt;t<T-dt;t++) n+=fabs(temp[t].med()/temp[t].err());
  n/=T-2*dt;
  
  return n;
}

jvec load_corr_all(oper_t &so,oper_t &si)
{
  jvec out(T,njacks);
  out=0;
  int nave=0;
  
  //check spin
  if(so.spin!=si.spin) crash("spin of source %d does not agree wih spin of the sink %d",so.spin,si.spin);
  
  for(int iso=0;iso<so.ind->n;iso++)
    for(int isi=0;isi<si.ind->n;isi++)
      {
	jvec corr=so.ind->sign[iso]*si.ind->sign[isi]*
	  jvec_load(combine("%s2pts_corr",base_path).c_str(),
		    T,njacks,icombo(so.ind->gam[iso],so.sme,so.ind->der1[iso],so.ind->der2[iso],
				    si.ind->gam[isi],si.sme,si.ind->der1[isi],si.ind->der2[isi]));
	out+=corr;
	cout<<"corr "<<nave<<": "<<corr[1]<<endl;
	
	double n=ndev(corr);
	int iszero=(fabs(n)<3);
	if(iszero)
	  {
	    cout<<corr<<endl;
	    crash("seems zero: %d, %lg",iszero,n);
	  }
	nave++;
      }
  
  return out;
}

void prepare_table(const char* out_path,int nop,oper_t *ops,const char *plot_path)
{
  //load all the data
  jvec data[nop*nop];
  int set=0;
  ofstream out(plot_path);
  out<<"@type xydy"<<endl;
  for(int iop_so=0;iop_so<nop;iop_so++)
    for(int iop_si=0;iop_si<nop;iop_si++)
      {
	data[iop_so*nop+iop_si]=load_corr_all(ops[iop_so],ops[iop_si]);
	out<<effective_mass(data[iop_so*nop+iop_si].subset(0,TH))<<"@s"<<set++<<" legend \""
	   <<ops[iop_so].name[0]<<ops[iop_so].sme<<"_"<<ops[iop_si].name[0]<<ops[iop_si].sme<<"\"\n&"<<endl;
	
	cout<<data[iop_so*nop+iop_si][1]<<endl;
	
	if(iop_so*nop+iop_si==0) data[iop_so*nop+iop_si].write_to_binfile(out_path);
	else                   data[iop_so*nop+iop_si].append_to_binfile(out_path);
      }
}

int main(int narg,char **arg)
{
  if(narg<2) crash("Use %s path",arg[0]);
  base_path=arg[1];
  
  gam_map[ 1]=0;
  gam_map[ 2]=1;
  gam_map[ 3]=2;
  gam_map[ 5]=3;
  gam_map[ 9]=4;
  gam_map[10]=5;
  gam_map[11]=6;
  gam_map[12]=7;
  gam_map[16]=8;
  gam_map[17]=9;
  gam_map[18]=10;
  
  //define indices mapping
  indices_t PLAIN,GI_DI,GI_BI,GI,EPSIJK_GJ_DK;
  PLAIN.add(0,0,0,0,+1);

  for(int i=0;i<3;i++)
    {
      GI.add(i+1,i+1,0,0,+1);
      GI_DI.add(i+1,i+1,i+1,0,+1);
      
      GI_BI.add(i+1,i+1,1+(i+1)%3,1+(i+2)%3,+1);
      GI_BI.add(i+1,i+1,1+(i+2)%3,1+(i+1)%3,-1);
      
      EPSIJK_GJ_DK.add(i+1,1+(i+1)%3,1+(i+2)%3,0,+1);
      EPSIJK_GJ_DK.add(i+1,1+(i+2)%3,1+(i+1)%3,0,-1);
    }

  ///////////////////////////////// 0-+ (pseudoscalar) //////////////////////////

  oper_t P5_LOC(SPIN0,PLAIN,"P5",0 ,0);
  oper_t P5_SME(SPIN0,PLAIN,"P5",SME,0);
  oper_t A0_LOC(SPIN0,PLAIN,"A0",0 ,0);
  oper_t A0_SME(SPIN0,PLAIN,"A0",SME,0);
  oper_t B1XD_A1_LOC(SPIN0,GI_DI,"C%d",0 ,0);
  oper_t B1XD_A1_SME(SPIN0,GI_DI,"C%d",SME,0);
  oper_t RHOXB_A1_LOC(SPIN0,GI_BI,"V%d",0 ,0);
  oper_t RHOXB_A1_SME(SPIN0,GI_BI,"V%d",SME ,0);
  oper_t RHO2XB_A1_LOC(SPIN0,GI_BI,"A%d",0 ,0);
  oper_t RHO2XB_A1_SME(SPIN0,GI_BI,"A%d",SME ,0);
  
  const int nop_pseudoscalar=10;
  oper_t ops_pseudoscalar[nop_pseudoscalar]={P5_LOC,A0_LOC,B1XD_A1_LOC,RHOXB_A1_LOC,RHO2XB_A1_LOC,
					     P5_SME,A0_SME,B1XD_A1_SME,RHOXB_A1_SME,RHO2XB_A1_SME};
  
  prepare_table("ops_pseudoscalar",nop_pseudoscalar,ops_pseudoscalar,combine("%s/out_pseudoscalar.xmg",base_path).c_str());
  
  /////////////////////////////////// 1-- (vector) ///////////////////////////////
  
  cout<<"prepare vector"<<endl;
  oper_t VK_LOC(SPIN1,GI,"V%d",0,0);
  oper_t VK_SME(SPIN1,GI,"V%d",SME,0);
  oper_t TK_LOC(SPIN1,GI,"T%d",0,0);
  oper_t TK_SME(SPIN1,GI,"T%d",SME,0);
  oper_t A0XD_T1_LOC(SPIN1,GI_DI,"S0",0,0);
  oper_t A0XD_T1_SME(SPIN1,GI_DI,"S0",SME,0);
  oper_t A1XD_T1_LOC(SPIN1,EPSIJK_GJ_DK,"A%d",0 ,0);
  oper_t A1XD_T1_SME(SPIN1,EPSIJK_GJ_DK,"A%d",SME,0);
  const int nop_vectorial=8;
  //oper_t ops_vectorial[nop_vectorial]={A1XD_T1_LOC,A1XD_T1_SME};
  oper_t ops_vectorial[nop_vectorial]={VK_LOC,TK_LOC,A0XD_T1_LOC,A1XD_T1_LOC,VK_SME,TK_SME,A0XD_T1_SME,A1XD_T1_SME};
  prepare_table("ops_vectorial",nop_vectorial,ops_vectorial,combine("%s/out_vectorial.xmg",base_path).c_str());
  
  return 0;
}
