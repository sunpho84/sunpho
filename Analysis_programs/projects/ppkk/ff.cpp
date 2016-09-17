#include "../../src/include.h"

const int ntseps=3;
const int tsep[ntseps]={16,20,24};
const int T=128,L=64;
const int njacks=77;
const int nmoms=3;
const int moms[nmoms][3]={{0,0,0},{0,0,1},{0,1,1}};
const int nmels=4;
const int iS=0,iV0=1,iVK=2,iTK=3;
const char tag_mel[nmels][3]={"S","V0","VK","TK"};
const char tag_mom[nmoms][10]={"(0,0,0)","(0,0,1)","(0,1,1)"};

#define LIGHT "ml_0.000678"
#define STRANGE "ms_0.02661"

//return the energy
jack latt_en(jack M,const int *mom)
{
  jack E=sqr(sinh(M/2));
  for(int i=0;i<3;i++) E+=sqr(sin(2*M_PI*mom[i]/L/2));
  return 2*asinh(sqrt(E));
}
jack cont_en(jack M,const int *mom)
{
  jack E=sqr(M);
  for(int i=0;i<3;i++) E+=sqr(2*M_PI*mom[i]/L);
  return sqrt(E);
}

void plot_all_tseps(const char *path,jvec *topl,jack &ref,const char *tag_y)
{
  ofstream fout_eff(path);
  int is=0;
  int symb_type[3]={1,2,3};
  int symb_color[3]={4,2,14};
  double tmin=1,tmax=0;
  //find ymin max
  double ymin=1e300,ymax=-1e300;
  for(int itsep=0;itsep<ntseps;itsep++)
    for(int t=2;t<tsep[itsep]-2;t++)
      {
	ymin=std::min(ymin,topl[itsep][t].med());
	ymax=std::max(ymax,topl[itsep][t].med());
      }
  double ydiff=ymax-ymin;
  double wmin=ref.med()-ydiff/2;
  double wmax=ref.med()+ydiff/2;
  for(int itsep=0;itsep<ntseps;itsep++)
    {
      fout_eff<<"@type xydy"<<endl;
      fout_eff<<"@s"<<is<<" legend \"t\\ssep\\N="<<tsep[itsep]<<"\""<<endl;
      fout_eff<<"@s"<<is<<" line type 0"<<endl;
      fout_eff<<"@s"<<is<<" symbol "<<symb_type[itsep]<<endl;
      fout_eff<<"@s"<<is<<" symbol size 0.5"<<endl;
      fout_eff<<"@s"<<is<<" symbol fill pattern 1"<<endl;
      fout_eff<<"@s"<<is<<" symbol color "<<symb_color[itsep]<<endl;
      fout_eff<<"@s"<<is<<" symbol fill color "<<symb_color[itsep]<<endl;
      fout_eff<<"@s"<<is<<" errorbar color "<<symb_color[itsep]<<endl;
      fout_eff<<"@s"<<is<<" symbol linewidth 2.0"<<endl;
      fout_eff<<"@s"<<is<<" errorbar linewidth 2.0"<<endl;
      fout_eff<<"@s"<<is<<" errorbar riser linewidth 2.0"<<endl;
      
      //print and find tmin/max for the plot
      for(int t=2;t<tsep[itsep]-2;t++)
	{
	  fout_eff<<(double)t/tsep[itsep]<<" "<<topl[itsep][t]<<endl;
	  if(topl[itsep][t].med()>wmin && topl[itsep][t].med()<wmax)
	    {
	      tmin=std::min(tmin,t/(double)tsep[itsep]);
	      tmax=std::max(tmax,t/(double)tsep[itsep]);
	    }
	}
      fout_eff<<"&"<<endl;is++;
    }
  
  for(int itsep=0;itsep<ntseps;itsep++)
    {
      fout_eff<<write_constant_with_error(ref,0,1);
      fout_eff<<"@s"<<is<<" linewidth 2"<<endl;
      fout_eff<<"&"<<endl;is++;
    }
  
  fout_eff<<"@ xaxis label \"t/t\\ssep\\N\""<<endl;
  fout_eff<<"@ xaxis label char size 1.5"<<endl;
  fout_eff<<"@ legend 0.832843137255, 0.406617647059"<<endl;
  fout_eff<<"@ yaxis label \""<<tag_y<<"\""<<endl;
  fout_eff<<"@ yaxis label char size 1.5"<<endl;

  fout_eff<<"@    world "<<tmin-(tmax-tmin)/20<<", "<<wmin<<", "<<tmax+(tmax-tmin)/20<<", "<<wmax<<endl;
  fout_eff<<"&"<<endl;is+=2;
}

int main(int narg,char **arg)
{
  //load 2pts
  jvec cPi(T,njacks);
  jvec cK(T,njacks);
  cPi.load("corr_Z2_l64t128_Ls12_m0.000678_LL_m0.000678_LL_15_15",0);
  cK.load("corr_Z2_l64t128_Ls12_m0.000678_LL_m0.02661_LL_15_15",0);
  cPi=cPi.simmetrized(1)/2;
  cK=cK.simmetrized(1)/2;
  cPi.clusterize();
  cK.clusterize();
  
  //fit masses
  jack MPi(njacks),MK(njacks);
  jack Z2Pi(njacks),Z2K(njacks);
  two_pts_fit(MPi,Z2Pi,cPi,15,64,"Pion.xmg");
  two_pts_fit(MK,Z2K,cK,22,64,"Kaon.xmg");
  
  //compute energies
  jvec EPi(nmoms,njacks),EK(nmoms,njacks);
  for(int imom=0;imom<nmoms;imom++)
    {
      EPi[imom]=latt_en(MPi,moms[imom]);
      EK[imom]=latt_en(MK,moms[imom]);
    }
  
  const int idec=1;
  char rest[3][20]={LIGHT,STRANGE,STRANGE};
  char deca[3][20]={LIGHT,LIGHT,STRANGE};
  char spect[3][20]={LIGHT,LIGHT,LIGHT};
  //load
  jvec pp[ntseps][nmoms][nmels];
  for(int itsep=0;itsep<ntseps;itsep++)
    {
      FILE *fin=open_file(combine("%d/3pt_%s_%s_%s",tsep[itsep],spect[idec],rest[idec],deca[idec]).c_str(),"r");
      //FILE *fin=open_file(combine("%d/3pt_ml_0.000678_ml_0.000678_ml_0.000678",tsep[itsep]).c_str(),"r");
      int i=0;
      for(int imom=0;imom<nmoms;imom++)
	for(int imel=0;imel<nmels;imel++)
	  {
	    pp[itsep][imom][imel]=jvec(tsep[itsep]+1,njacks);
	    pp[itsep][imom][imel].load(fin,i++);
	  }
    }
  
  for(int imel=0;imel<nmels;imel++)
    for(int imom=0;imom<nmoms;imom++)
      {
	//find masses
	jack m1,m2,pre;
	switch (idec)
	  {
	  case 0:
	    m1=MPi;
	    m2=EPi[imom];
	    break;
	  case 1:
	    m1=MPi;
	    m2=EK[imom];
	    break;
	  case 2:
	    m1=MK;
	    m2=EK[imom];
	    break;
	  }
	pre=m1-m2;
	
	//compute aperiodic effective mass and mel corr
	jvec aper[ntseps],mel_corr[ntseps];
	for(int itsep=0;itsep<ntseps;itsep++)
	  {
	    aper[itsep]=aperiodic_effective_mass(pp[itsep][imom][imel]);
	    mel_corr[itsep]=pp[itsep][imom][imel];
	    for(int t=0;t<tsep[itsep];t++)
	      mel_corr[itsep][t]/=exp(-m1*t-m2*(tsep[itsep]-t));
	  }
	
	//find the point closer to the correct behaviour
	int fit_tmin[ntseps];
	jack fit_val[ntseps];
	for(int itsep=0;itsep<ntseps;itsep++)
	  {
	    double dmin=1e300;
	    fit_tmin[itsep]=0;
	    for(int t=1;t<tsep[itsep]-1;t++)
	      {
		jack d=pre-aper[itsep][t];
		double dd=fabs(d.med()/d.err());
		if(dd<dmin)
		  {
		    dmin=dd;
		    fit_tmin[itsep]=t;
		  }
	      }
	    fit_val[itsep]=mel_corr[itsep][fit_tmin[itsep]];
	    cout<<"tmin: "<<fit_tmin[itsep]<<" "<<fit_val[itsep]<<" "<<dmin<<endl;
	  }
	
	//plot the comparison of all eff mass
	plot_all_tseps(combine("plots/eff_mel_%s_mom%d.xmg",tag_mel[imel],imom).c_str(),aper,pre,combine("%s_%s",tag_mel[imel],tag_mom[imom]).c_str());
	//plot the matrix element correlators
	plot_all_tseps(combine("plots/mel_%s_mom%d.xmg",tag_mel[imel],imom).c_str(),mel_corr,fit_val[ntseps-1],combine("%s_%s",tag_mel[imel],tag_mom[imom]).c_str());
      }
  
  return 0;
}
