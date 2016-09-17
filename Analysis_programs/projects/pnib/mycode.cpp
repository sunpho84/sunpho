#include "include.h"

//#define DEBUG

int T;
int njacks;
int nqmass=1;
int nsm_sink=2;
int nr=2;
enum dir_exc_t{DIR,EXC};
enum reim_t{RE,IM};
enum loc_sme_t{LOC,SME};
enum icombo_t{i000,iS00,i0S0,i00S,iP00,i0P0,i00P,iT00,i0T0,i00T,i200,i020,i002,i011,i101,i110};

int ind_corr(icombo_t icombo,loc_sme_t ism_sink,int ima,int ra,int imb,int rb,int imc,int rc,dir_exc_t dir_exc,reim_t ri)
{return
    (ri+2*
     (dir_exc+2*
      (rc+nr*
       (imc+nqmass*
	(rb+nr*
	 (imb+nqmass*
	  (ra+nr*
	   (ima+nqmass*
	    (ism_sink+nsm_sink*icombo)))))))));
}

jvec pure_load(icombo_t icombo,loc_sme_t ism_sink,int ima,int ra,int imb,int rb,int imc,int rc,dir_exc_t dir_exc,reim_t ri)
{return jvec_load("bar_contr",T,njacks,ind_corr(icombo,ism_sink,ima,ra,imb,rb,imc,rc,dir_exc,ri));}

jvec load_nucleon(icombo_t icombo,loc_sme_t ism_sink,int ra,int rb,int rc,reim_t ri)
{return pure_load(icombo,ism_sink,0,ra,0,rb,0,rc,DIR,ri)-pure_load(icombo,ism_sink,0,ra,0,rb,0,rc,EXC,ri);}

#ifdef DEBUG
//int ra=1,rb=0,rc=1;
int ra=0,rb=0,rc=0;
#else
int ra=1,rb=1,rc=1;
#endif

jvec load_nucleon(icombo_t icombo,loc_sme_t ism_sink,reim_t ri,string path="")
{
#ifdef DEBUG
  jvec out=load_nucleon(icombo,ism_sink, ra, rb, rc,ri);
#else
  jvec out=(load_nucleon(icombo,ism_sink, ra, rb, rc,ri)+
	    load_nucleon(icombo,ism_sink,!ra,!rb,!rc,ri)).simmetrized(1)/2;
#endif
  
  if(path!="") out.print_to_file(path.c_str());
  return out;
}

int main()
{
  string ou=combine("plots_%d%d%d",ra,rb,rc);
  
  FILE *input=open_file("input","r");
  read_formatted_from_file_expecting((char*)(&T),input,"%d","T");
  read_formatted_from_file_expecting((char*)(&njacks),input,"%d","NJacks");
  fclose(input);
  cout<<"--"<<endl;
  jvec nuc_000=load_nucleon(i000,SME,RE,ou+"/000.xmg");
  cout<<"--"<<endl;
  
  jvec nuc_S00=load_nucleon(iS00,SME,RE,ou+"/S00.xmg");
  jvec nuc_0S0=load_nucleon(i0S0,SME,RE,ou+"/0S0.xmg");
  jvec nuc_00S=load_nucleon(i00S,SME,RE,ou+"/00S.xmg");
  
  jvec nuc_mass_slope=(nuc_S00+nuc_00S-nuc_0S0)/nuc_000;
  cout<<"Mnuc: "<<smart_print(constant_fit(aperiodic_effective_mass(nuc_000),13,20,ou+"/effective_mass.xmg"))<<endl;
  numerical_derivative(nuc_mass_slope).print_to_file(ou+"/mass_slope.xmg");
  
  jvec nuc_P00=load_nucleon(iP00,SME,IM,ou+"/P00.xmg");
  jvec nuc_0P0=load_nucleon(i0P0,SME,IM,ou+"/0P0.xmg");
  jvec nuc_00P=load_nucleon(i00P,SME,IM,ou+"/00P.xmg");
  
  jvec nuc_200=load_nucleon(i200,SME,RE,ou+"/200.xmg");
  jvec nuc_020=load_nucleon(i020,SME,RE,ou+"/020.xmg");
  jvec nuc_002=load_nucleon(i002,SME,RE,ou+"/002.xmg");
  
  jvec nuc_T00=load_nucleon(iT00,SME,RE,ou+"/T00.xmg");
  jvec nuc_0T0=load_nucleon(i0T0,SME,RE,ou+"/0T0.xmg");
  jvec nuc_00T=load_nucleon(i00T,SME,RE,ou+"/00T.xmg");
  
  jvec nuc_011=load_nucleon(i011,SME,RE,ou+"/011.xmg");
  jvec nuc_101=load_nucleon(i101,SME,RE,ou+"/101.xmg");
  jvec nuc_110=load_nucleon(i110,SME,RE,ou+"/110.xmg");
  
  pure_load(i000,SME,0,ra,0,rb,0,rc,DIR,RE).print_to_file(ou+"/000_dir.xmg");
  pure_load(i000,SME,0,ra,0,rb,0,rc,EXC,RE).print_to_file(ou+"/000_exc.xmg");
  
  pure_load(iS00,SME,0,ra,0,rb,0,rc,DIR,RE).print_to_file(ou+"/S00_dir.xmg");
  pure_load(i0S0,SME,0,ra,0,rb,0,rc,DIR,RE).print_to_file(ou+"/0S0_dir.xmg");
  pure_load(i00S,SME,0,ra,0,rb,0,rc,DIR,RE).print_to_file(ou+"/00S_dir.xmg");
  pure_load(iS00,SME,0,ra,0,rb,0,rc,EXC,RE).print_to_file(ou+"/S00_exc.xmg");
  pure_load(i0S0,SME,0,ra,0,rb,0,rc,EXC,RE).print_to_file(ou+"/0S0_exc.xmg");
  pure_load(i00S,SME,0,ra,0,rb,0,rc,EXC,RE).print_to_file(ou+"/00S_exc.xmg");
  
  pure_load(iT00,SME,0,ra,0,rb,0,rc,DIR,RE).print_to_file(ou+"/T00_dir.xmg");
  pure_load(i0T0,SME,0,ra,0,rb,0,rc,DIR,RE).print_to_file(ou+"/0T0_dir.xmg");
  pure_load(i00T,SME,0,ra,0,rb,0,rc,DIR,RE).print_to_file(ou+"/00T_dir.xmg");
  pure_load(iT00,SME,0,ra,0,rb,0,rc,EXC,RE).print_to_file(ou+"/T00_exc.xmg");
  pure_load(i0T0,SME,0,ra,0,rb,0,rc,EXC,RE).print_to_file(ou+"/0T0_exc.xmg");
  pure_load(i00T,SME,0,ra,0,rb,0,rc,EXC,RE).print_to_file(ou+"/00T_exc.xmg");
  
  pure_load(iP00,SME,0,ra,0,rb,0,rc,DIR,IM).print_to_file(ou+"/P00_dir.xmg");
  pure_load(i0P0,SME,0,ra,0,rb,0,rc,DIR,IM).print_to_file(ou+"/0P0_dir.xmg");
  pure_load(i00P,SME,0,ra,0,rb,0,rc,DIR,IM).print_to_file(ou+"/00P_dir.xmg");
  pure_load(iP00,SME,0,ra,0,rb,0,rc,EXC,IM).print_to_file(ou+"/P00_exc.xmg");
  pure_load(i0P0,SME,0,ra,0,rb,0,rc,EXC,IM).print_to_file(ou+"/0P0_exc.xmg");
  pure_load(i00P,SME,0,ra,0,rb,0,rc,EXC,IM).print_to_file(ou+"/00P_exc.xmg");
  
  pure_load(i200,SME,0,ra,0,rb,0,rc,DIR,RE).print_to_file(ou+"/200_dir.xmg");
  pure_load(i020,SME,0,ra,0,rb,0,rc,DIR,RE).print_to_file(ou+"/020_dir.xmg");
  pure_load(i002,SME,0,ra,0,rb,0,rc,DIR,RE).print_to_file(ou+"/002_dir.xmg");
  pure_load(i200,SME,0,ra,0,rb,0,rc,EXC,RE).print_to_file(ou+"/200_exc.xmg");
  pure_load(i020,SME,0,ra,0,rb,0,rc,EXC,RE).print_to_file(ou+"/020_exc.xmg");
  pure_load(i002,SME,0,ra,0,rb,0,rc,EXC,RE).print_to_file(ou+"/002_exc.xmg");
 
  pure_load(i011,SME,0,ra,0,rb,0,rc,DIR,RE).print_to_file(ou+"/011_dir.xmg");
  pure_load(i101,SME,0,ra,0,rb,0,rc,DIR,RE).print_to_file(ou+"/101_dir.xmg");
  pure_load(i110,SME,0,ra,0,rb,0,rc,DIR,RE).print_to_file(ou+"/110_dir.xmg");
  pure_load(i011,SME,0,ra,0,rb,0,rc,EXC,RE).print_to_file(ou+"/011_exc.xmg");
  pure_load(i101,SME,0,ra,0,rb,0,rc,EXC,RE).print_to_file(ou+"/101_exc.xmg");
  pure_load(i110,SME,0,ra,0,rb,0,rc,EXC,RE).print_to_file(ou+"/110_exc.xmg");
  
  aperiodic_effective_mass(load_nucleon(i000,SME,0,0,0,RE)+load_nucleon(i000,SME,1,1,1,RE)).print_to_file(ou+"/test_000.xmg");
  aperiodic_effective_mass(load_nucleon(i000,SME,0,0,1,RE)+load_nucleon(i000,SME,1,1,0,RE)).print_to_file(ou+"/test_001.xmg");
  aperiodic_effective_mass(load_nucleon(i000,SME,0,1,1,RE)+load_nucleon(i000,SME,1,0,0,RE)).print_to_file(ou+"/test_011.xmg");
  aperiodic_effective_mass(load_nucleon(i000,SME,1,0,1,RE)+load_nucleon(i000,SME,0,1,0,RE)).print_to_file(ou+"/test_101.xmg");
  
  jvec_load("mes_contr_00",T,njacks,0).print_to_file(ou+"/pion_00_r1.xmg");
  jvec_load("mes_contr_LL",T,njacks,0).print_to_file(ou+"/pion_LL_r1.xmg");
  jvec_load("mes_contr_0M",T,njacks,0).print_to_file(ou+"/pion_0M_r1.xmg");
  
  jvec_load("mes_contr_00",T,njacks,1).print_to_file(ou+"/pion_00_s1.xmg");
  jvec_load("mes_contr_LL",T,njacks,1).print_to_file(ou+"/pion_LL_s1.xmg");
  jvec_load("mes_contr_0M",T,njacks,1).print_to_file(ou+"/pion_0M_s1.xmg");
  
  jvec_load("mes_contr_00",T,njacks,2).print_to_file(ou+"/pion_00_r2.xmg");
  jvec_load("mes_contr_LL",T,njacks,2).print_to_file(ou+"/pion_LL_r2.xmg");
  jvec_load("mes_contr_0M",T,njacks,2).print_to_file(ou+"/pion_0M_r2.xmg");
  
  return 0;
}
