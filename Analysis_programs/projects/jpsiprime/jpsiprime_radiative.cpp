#include <math.h>
#include "include.h"
#include "eff_bis.cpp"
#include "../nf2/common_pars.cpp"

char basepath[100],plotpath[100];//,outfile[100];
int ntheta_S0,ntheta_S1;
int ith_S0,ith_S1;
double theta;
int T,L,TH,TSEP;
int ibeta;
int njacks;
int nst_P5P5,nst_VKVK,nsm=2;
int fix_smlv;
int fix_ismlv;  
int min_GR_fit,max_GR_fit;
int min_FI_fit,max_FI_fit;
int tm_run;

int case01;

int REAL=0,IMAG=1;
int GR=0,EX=1;
int LC=0,SM=1;

const double Zvm[6]={0.5816,0.6103,0.6451,0.686,1.16,0.992};

void read_input(char *path)
{
  FILE *fin=open_file(path,"r");

  //read basepath
  read_formatted_from_file_expecting(basepath,fin,"%s","basepath");
  read_formatted_from_file_expecting(plotpath,fin,"%s","plotpath");
  //read_formatted_from_file_expecting(outfile,fin,"%s","outfile");
  
  read_formatted_from_file_expecting((char*)&tm_run,fin,"%d","TwistedMass");
  
  //eta'c->jpsi psi'->etac
  read_formatted_from_file_expecting((char*)&case01,fin,"%d","Case01");
  
  //read theta
  read_formatted_from_file_expecting((char*)&ntheta_S0,fin,"%d","nthetaS0");
  read_formatted_from_file_expecting((char*)&ntheta_S1,fin,"%d","nthetaS1");
  read_formatted_from_file_expecting((char*)&ith_S0,fin,"%d","ithS0");
  read_formatted_from_file_expecting((char*)&ith_S1,fin,"%d","ithS1");
  read_formatted_from_file_expecting((char*)&theta,fin,"%lg","theta");
  
  //read T
  read_formatted_from_file_expecting((char*)&T,fin,"%d","T");
  read_formatted_from_file_expecting((char*)&TSEP,fin,"%d","TSep");
  L=TH=T/2;
  
  //read beta
  read_formatted_from_file_expecting((char*)&ibeta,fin,"%d","Beta");
  read_formatted_from_file_expecting((char*)&njacks,fin,"%d","NJacks");
  
  //read n fitted states
  read_formatted_from_file_expecting((char*)&nst_P5P5,fin,"%d","NP5P5");
  read_formatted_from_file_expecting((char*)&nst_VKVK,fin,"%d","NVKVK");
  
  //read fixing
  read_formatted_from_file_expecting((char*)&fix_ismlv,fin,"%d","Fix_ismlv");
  read_formatted_from_file_expecting((char*)&fix_smlv,fin,"%d","Fix_smlv");
  
  //read min and max interval for GR and FI mel
  read_formatted_from_file_expecting((char*)&min_GR_fit,fin,"%d","GR_fit_int");
  read_formatted_from_file((char*)&max_GR_fit,fin,"%d","GR_fit_int");
  read_formatted_from_file_expecting((char*)&min_FI_fit,fin,"%d","FI_fit_int");
  read_formatted_from_file((char*)&max_FI_fit,fin,"%d","FI_fit_int");
  
  fclose(fin);
}

int iZ(int ist,int ism)
{return ist*nsm+ism;}

/////////////////////////// 3pts loading ////////////////////////////

int icombo_3pts(int ith1,int ith2,int reim)
{
  return reim+2*(ith1+ntheta_S0*ith2);
}

jvec load_3pts(const char *name,int ism1,int ism2,int ith1,int ith2,int reim)
{
  return jvec_load(combine("%s/3pts_sp0_%s_%02d_%02d",basepath,name,ism1,ism2).c_str(),T,njacks,icombo_3pts(ith1,ith2,reim))*Zvm[ibeta];
}

jvec load_VIVJ(int ism1,int ism2)
{
  jvec P5thVIVJ_pt1= load_3pts("VIVJ_pt1",ism1,ism2,ith_S0,ith_S1,IMAG);
  jvec P5thVIVJ_pt2=-load_3pts("VIVJ_pt2",ism1,ism2,ith_S0,ith_S1,IMAG);
  jvec P5thVIVJ=(P5thVIVJ_pt1+P5thVIVJ_pt2)/2;

  return (TSEP==TH)?P5thVIVJ.simmetrized(-1):P5thVIVJ.subset(0,TSEP+1);
}

///////////////////////////// M & Z loading /////////////////////////

void MZ_load(const char *name,jvec &M,jvec &Z,int nst)
{
  jvec temp(nst*(nsm+1),njacks);
  temp.load(name,0);
  
  int ient=0;
  //load M
  for(int ist=0;ist<nst;ist++)
    M[ist]=temp[ient++];

  for(int ist=0;ist<nst;ist++)
    for(int ism=0;ism<nsm;ism++)
      Z[iZ(ist,ism)]=temp[ient++];

}

////////////////////////////// creating dt //////////////////////////

jvec dt(jack ZV,jack ZP,jack MV,jack MP)
{
  jvec out(TSEP+1,njacks);
  
  for(int t=0;t<=TSEP;t++)
    out[t]=ZV*ZP*exp(-MV*t)*exp(-MP*(TSEP-t))/(2*MP*2*MV);
  
  return out;
}

/////////////////// print the ratio of 1st and gr state //////////////////////////

void study_operator_content(jack ZV,jack ZP,jack MV,jack MP,int nst)
{
  jvec OP[nsm*nst],OV[nsm*nst];
  for(int ist=0;ist<nst;ist++)
    for(int ism=0;ism<nsm;ism++)
      {
	int ic=iZ(ist,ism);
	OP[ic]=jvec(TH+1,njacks);
	OV[ic]=jvec(TH+1,njacks);
	for(int t=0;t<=TH;t++)
	  {
	    OP[ic][t]=ZP[ic]*exp(-MP[ist]*t)/(2*MP[ist]);
	    OV[ic][t]=ZV[ic]*exp(-MV[ist]*t)/(2*MV[ist]);
	  }
      }
  
  for(int ism=0;ism<nsm;ism++)
    {
      ofstream outP(combine("/tmp/Pcont%d.xmg",ism).c_str());
      outP<<"@type xydy"<<endl;
      outP<<OP[iZ(EX,ism)]/OP[iZ(GR,ism)]<<endl;
      cout<<"P suppression "<<ism<<": "<<ZP[iZ(EX,ism)]/ZP[iZ(GR,ism)]<<endl;
      
      ofstream outV(combine("/tmp/Vcont%d.xmg",ism).c_str());
      outV<<"@type xydy"<<endl;
      outV<<OV[iZ(EX,ism)]/OV[iZ(GR,ism)]<<endl;
      cout<<"V suppression "<<ism<<": "<<ZV[iZ(EX,ism)]/ZV[iZ(GR,ism)]<<endl;
    }
}

int main(int narg,char **arg)
{
  if(narg<2) crash("Use %s input",arg[0]);
  read_input(arg[1]);
  
  //load three points
  jvec P5thVIVJ_00=(case01==0)?load_VIVJ(fix_smlv,00):load_VIVJ(00,fix_smlv);
  jvec P5thVIVJ_20=(case01==0)?load_VIVJ(fix_smlv,20):load_VIVJ(20,fix_smlv);
  
  if(case01==0)
    {
      P5thVIVJ_00.print_to_file("plots_radiative/3pts_smeared_local.xmg");
      P5thVIVJ_20.print_to_file("plots_radiative/3pts_smeared_smeared.xmg");
  
      jvec MP(nst_P5P5,njacks),MV(nst_VKVK,njacks),EV(nst_VKVK,njacks);
      jvec ZP(nsm*nst_P5P5,njacks),ZV(nsm*nst_VKVK,njacks),ZV_moving(nsm*nst_VKVK,njacks);
      MZ_load("MZ_P_00_20",MP,ZP,nst_P5P5);
      MZ_load("MZ_V_00_20",MV,ZV,nst_VKVK);
      MZ_load("EZ_V_00_20",EV,ZV_moving,nst_VKVK);
      jack diff=MV[0]-MP[0];
      cout<<MP<<MV<<EV<<endl;
      
      //takes the two orthogonal combinations
      jvec P5thVIVJ_GR_BEST=ZP[iZ(EX,SM)]*P5thVIVJ_00-ZP[iZ(EX,LC)]*P5thVIVJ_20;
      jvec P5thVIVJ_FI_BEST=ZP[iZ(GR,SM)]*P5thVIVJ_00-ZP[iZ(GR,LC)]*P5thVIVJ_20;
      
      //check relevance of contribution of eta' in smeared correlator
      cout<<"Eta' contribution in smeared 2pts correlator: "<<ZP[iZ(EX,SM)]/ZP[iZ(GR,SM)]<<endl;
      cout<<"Eta' contribution in smeared 3pts correlator: "<<
	ZP[iZ(EX,SM)]*ZP[iZ(GR,LC)]/(ZP[iZ(EX,LC)]*ZP[iZ(GR,SM)])<<endl;
      
      //compute matrixe elements
      jvec C_MEL_GR=P5thVIVJ_GR_BEST/
	dt(ZV_moving[iZ(GR,fix_ismlv)],ZP[iZ(GR,LC)]*ZP[iZ(EX,SM)]-ZP[iZ(GR,SM)]*ZP[iZ(EX,LC)],EV[GR],MP[EX]);
      jvec C_MEL_FI=P5thVIVJ_FI_BEST/
	dt(ZV_moving[iZ(GR,fix_ismlv)],ZP[iZ(EX,LC)]*ZP[iZ(GR,SM)]-ZP[iZ(EX,SM)]*ZP[iZ(GR,LC)],EV[GR],MP[EX]);
      
      //fit over a fixed range
      jack MEL_GR=constant_fit(C_MEL_GR,min_GR_fit,max_GR_fit);
      jack MEL_FI=constant_fit(C_MEL_FI,min_FI_fit,max_FI_fit);
    
      {
	ofstream out(combine("%s/out.xmg",plotpath).c_str());
	int iset=0;
	
	out<<"@type xydy"<<endl;
	out<<"@s"<<iset<<" legend \"00\""<<endl;
	out<<P5thVIVJ_00/dt(ZV_moving[iZ(GR,fix_ismlv)],ZP[iZ(GR,LC)],EV[GR],MP[GR])<<endl;
	out<<"&"<<endl;
	iset++;
	
	out<<"@type xydy"<<endl;
	out<<"@s"<<iset<<" legend \"20\""<<endl;
	out<<P5thVIVJ_20/dt(ZV_moving[iZ(GR,fix_ismlv)],ZP[iZ(GR,SM)],EV[GR],MP[GR])<<endl;
	out<<"&"<<endl;
	iset++;
	
	out<<"@type xydy"<<endl;
	out<<"@s"<<iset<<" legend \"best ground\""<<endl;
	out<<write_constant_fit_plot(C_MEL_GR,MEL_GR,min_GR_fit,max_GR_fit,iset);
	out<<"&"<<endl;
	iset+=3;
	
	out<<"@type xydy"<<endl;
	out<<"@s"<<iset<<" legend \"best first\""<<endl;
	out<<"@s"<<iset<<" line color 7"<<endl;
	out<<"@s"<<iset<<" errorbar color 7"<<endl;
	out<<write_constant_fit_plot(C_MEL_FI,MEL_FI,min_FI_fit,max_FI_fit,iset);
	out<<"&"<<endl;
	iset+=3;
      }
      
      //compute Q2
      jack Q2=sqr(MP[GR]-EV[GR])-3*sqr(theta*M_PI/TH);
      cout<<"Q2 for ground state: "<<smart_print(Q2)<<endl;
      Q2=sqr(MP[EX]-EV[GR])-3*sqr(theta*M_PI/TH);
      cout<<"Q2 for excited state: "<<smart_print(Q2)<<endl;
      
      //testing f.f. for ground state: this is the factor to put
      double qi=M_PI*theta/TH;
      jack f=1/qi*(MV[GR]+MP[GR])/(2*MP[GR]);
      cout<<"f: "<<f<<", ff: "<<f*MEL_GR<<endl;
      //testing f.f. for excited state
      f=1/qi*(MV[GR]+MP[EX])/(2*MP[GR]);
      cout<<"f: "<<f<<", ff: "<<f*MEL_FI<<endl;
    }
  else
    {
      P5thVIVJ_00.print_to_file("plots_radiative/3pts_local_smeared.xmg");
      P5thVIVJ_20.print_to_file("plots_radiative/3pts_smeared_smeared.xmg");
  
      jvec MV(nst_VKVK,njacks),MP(nst_P5P5,njacks),EP(nst_P5P5,njacks);
      jvec ZV(nsm*nst_VKVK,njacks),ZP(nsm*nst_P5P5,njacks),ZP_moving(nsm*nst_P5P5,njacks);
      MZ_load("MZ_V_00_20",MV,ZV,nst_VKVK);
      MZ_load("MZ_P_00_20",MP,ZP,nst_P5P5);
      MZ_load("EZ_P_00_20",EP,ZP_moving,nst_P5P5);
      jack diff=MP[GR]-MV[GR];
      cout<<MV<<MP<<EP<<endl;
      cout<<"MV'/MV: "<<MV[EX]/MV[GR]<<endl;
      
      //takes the two orthogonal combinations
      jvec P5thVIVJ_GR_BEST=ZV[iZ(EX,SM)]*P5thVIVJ_00-ZV[iZ(EX,LC)]*P5thVIVJ_20;
      jvec P5thVIVJ_FI_BEST=ZV[iZ(GR,SM)]*P5thVIVJ_00-ZV[iZ(GR,LC)]*P5thVIVJ_20;
      
      //check relevance of contribution of psi' in smeared correlator
      cout<<"Psi' contribution in smeared 2pts correlator: "<<ZV[iZ(EX,SM)]/ZV[iZ(GR,SM)]<<endl;
      cout<<"Psi' contribution in smeared 3pts correlator: "<<ZV[iZ(EX,SM)]*ZV[iZ(GR,LC)]
	/(ZV[iZ(EX,LC)]*ZV[iZ(GR,SM)])<<endl;
      
      //compute matrixe elements
      jvec C_MEL_GR=P5thVIVJ_GR_BEST/
	dt(ZP_moving[iZ(GR,fix_ismlv)],ZV[iZ(GR,LC)]*ZV[iZ(EX,SM)]-ZV[iZ(GR,SM)]*ZV[iZ(EX,LC)],MV[GR],EP[GR]);
      jvec C_MEL_FI=P5thVIVJ_FI_BEST/							         
	dt(ZP_moving[iZ(GR,fix_ismlv)],ZV[iZ(EX,LC)]*ZV[iZ(GR,SM)]-ZV[iZ(EX,SM)]*ZV[iZ(GR,LC)],MV[EX],EP[GR]);
      
      {
	ofstream out(combine("%s/dtp_fr_dt.xmg",plotpath).c_str());
	out<<"@type xydy"<<endl;
	out<<dt(ZP_moving[iZ(GR,fix_ismlv)],ZV[iZ(EX,LC)]*ZV[iZ(GR,SM)]-ZV[iZ(EX,SM)]*ZV[iZ(GR,LC)],MV[EX],EP[GR])/
	  dt(ZP_moving[iZ(GR,fix_ismlv)],ZV[iZ(GR,LC)]*ZV[iZ(EX,SM)]-ZV[iZ(GR,SM)]*ZV[iZ(EX,LC)],MV[GR],EP[GR]);
      }

      //fit over a fixed range
      jack MEL_GR=constant_fit(C_MEL_GR,min_GR_fit,max_GR_fit);
      jack MEL_FI=constant_fit(C_MEL_FI,min_FI_fit,max_FI_fit);
    
      {
	ofstream out(combine("%s/out.xmg",plotpath).c_str());
	int iset=0;
	
	out<<"@type xydy"<<endl;
	out<<"@s"<<iset<<" legend \"00\""<<endl;
	out<<P5thVIVJ_00/dt(ZP_moving[iZ(GR,fix_ismlv)],ZV[iZ(GR,LC)],MV[GR],EP[GR])<<endl;
	out<<"&"<<endl;
	iset++;
	
	out<<"@type xydy"<<endl;
	out<<"@s"<<iset<<" legend \"20\""<<endl;
	out<<P5thVIVJ_20/dt(ZP_moving[iZ(GR,fix_ismlv)],ZV[iZ(GR,SM)],MV[GR],EP[GR])<<endl;
	out<<"&"<<endl;
	iset++;
	
	out<<"@type xydy"<<endl;
	out<<"@s"<<iset<<" legend \"best ground\""<<endl;
	out<<write_constant_fit_plot(C_MEL_GR,MEL_GR,min_GR_fit,max_GR_fit,iset);
	out<<"&"<<endl;
	iset+=3;
	
	out<<"@type xydy"<<endl;
	out<<"@s"<<iset<<" legend \"best first\""<<endl;
	out<<"@s"<<iset<<" line color 7"<<endl;
	out<<"@s"<<iset<<" errorbar color 7"<<endl;
	out<<write_constant_fit_plot(C_MEL_FI,MEL_FI,min_FI_fit,max_FI_fit,iset);
	out<<"&"<<endl;
	iset+=3;
      }
      
      cout<<"best theta for exc: "<<(MV[EX]*MV[EX]-MP[GR]*MP[GR])/(2*MV[EX])*L/sqrt(3)/M_PI<<endl;
      
      //compute Q2
      jack Q2=sqr(MV[GR]-EP[GR])-3*sqr(theta*M_PI/TH);
      cout<<"Q2 for ground state: "<<smart_print(Q2)<<endl;
      Q2=sqr(MV[EX]-EP[GR])-3*sqr(theta*M_PI/TH);
      cout<<"Q2 for excited state: "<<smart_print(Q2)<<endl;
      
      //testing f.f. for ground state: this is the factor to put
      double qi=M_PI*theta/TH;
      jack f=1/qi*(MP[GR]+MV[GR])/(2*MV[GR]);
      cout<<"f: "<<f<<", ff: "<<f*MEL_GR<<endl;
      //testing f.f. for excited state
      f=1/qi*(MP[GR]+MV[EX])/(2*MV[GR]);
      cout<<"f: "<<f<<", ff: "<<f*MEL_FI<<endl;
    }
  
  return 0;
}
