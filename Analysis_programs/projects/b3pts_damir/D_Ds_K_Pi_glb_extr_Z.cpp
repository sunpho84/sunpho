#include "include.h"
#include "../nf2/common_pars.cpp"

int nth,nens,nth_to_use;
const int nth_max=8,nens_max=13,njacks=16;
int ibeta[nens_max],use[nens_max];
double lmass[nens_max];
int iboot;
int stranged;
int ndof;

//double lat_internal[4]={1/2.0275,1/2.338,1/2.9657,1/3.679};

const double MDstar=2.01,MD0star=2.318;
const double MDSstar=2.1123,MD0Sstar=2.3178;
const double MD_ph=1.8696,MDS_ph=1.9685,MK_ph=0.4937,MP_ph=0.13957;
const double Q2_max_ph[2]={sqr(MD_ph-MP_ph),sqr(MD_ph-MK_ph)};
bvec Q2,FP,F0,FT,MD,MP,EP,Z;

bvec ml;

int study_T_fr_P=1;

template <class T1,class T2> T1 fun_Z(T1 MP,T1 MD,T2 Q2)
{
  double T0=0;//(MD_ph+MP_ph)*sqr(sqrt(MD_ph)-sqrt(MP_ph));
  T1 TP=sqr(MD+MP);
  T1 Z=sqrt(TP-Q2)-sqrt(TP-T0);
  Z/=sqrt(TP-Q2)+sqrt(TP-T0);
  
  return Z;
}

template <class T1,class T2> T1 fun_Zone(T1 MP,T1 MD,T2 Z)
{
  double T0=0;//(MD_ph+MP_ph)*sqr(sqrt(MD_ph)-sqrt(MP_ph));
  T1 TP=sqr(MD+MP),B=sqrt(TP-T0),A=B*(1+Z)/(1-Z);
  
  return TP-sqr(A);
}

int icombo(int iens,int ith)
{return iens*nth_to_use+ith;}

template <class T1,class T2,class T3> void put_or_remove_pole(T1 &fP,T1 &f0,T1 &fT,T2 q2,T3 MD,int put_remove)
{
  int use_diff=1;
  int pole_0=0;
  
  T2 factP,fact0,factT;
  if(!stranged)
    {
      //cout<<MD-MD_ph<<endl;
      factP=(1-q2/sqr(MDstar+(MD-MD_ph)*use_diff));
      factT=(1-q2/sqr(MDstar+(MD-MD_ph)*use_diff));
      fact0=(1-q2/sqr(MD0star+(MD-MD_ph)*use_diff));
    }
  else
    {
      //cout<<MD-MDS_ph<<endl;
      factP=(1-q2/sqr(MDSstar+(MD-MDS_ph))*use_diff);
      factT=(1-q2/sqr(MDSstar+(MD-MDS_ph))*use_diff);
      fact0=(1-q2/sqr(MD0Sstar+(MD-MDS_ph))*use_diff);
    }
  if(put_remove==1)
    {
      fP*=factP;
      fT*=factT;
      if(pole_0) f0*=fact0;
      cout<<" FATTORE CORRETTIVO DA RIMUOVERE F0: "<<fact0<<" Q2: "<<q2<<", MD: "<<
	sqr((stranged?MD0Sstar:MD0star)+(MD-MDS_ph)*use_diff)<<endl;
    }
  else
    {
      fP/=factP;
      fT/=factT;
      if(pole_0) f0/=fact0;
      cout<<" FATTORE CORRETTIVO DA AGGIUNGERE F0: "<<fact0<<" Q2: "<<q2<<", MD: "<<
	sqr((stranged?MD0Sstar:MD0star)+(MD-MDS_ph)*use_diff)<<endl;
    }
}

template <class T1,class T2,class T3> void put_pole(T1 &fP,T1 &f0,T1 &fT,T2 q2,T3 MD)
{
  put_or_remove_pole(fP,f0,fT,q2,MD,0);
}
template <class T1,class T2,class T3> void remove_pole(T1 &fP,T1 &f0,T1 &fT,T2 q2,T3 MD)
{
  put_or_remove_pole(fP,f0,fT,q2,MD,1);
}

void boot_from_jack(boot &out,jack in,int *iboot_jack)
{
  int nboot=out.nboot;
  int njack=in.njack;
  for(iboot=0;iboot<nboot;iboot++) out.data[iboot]=in.data[iboot_jack[iboot]];
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

double *Y_fit,*err_Y_fit;

void chi2(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ndof=-npar;
  ch=0;
  for(int iens=0;iens<nens;iens++)
    if(use[iens])
      for(int ith=0;ith<nth_to_use;ith++)
	{
	  int ic=icombo(iens,ith);
	  
	  double z=Z[ic][iboot],a2=pow(lat[ibeta[iens]][iboot],2),l=ml[iens][iboot];
	  double Y_teo[3]={
	    (p[0]+p[1]*a2+p[2]*l)+z*(p[3]+p[4]*a2+p[5]*l),
	    (p[0]+p[1]*a2+p[2]*l)+z*(p[6]+p[7]*a2+p[8]*l),
	    (p[9]+p[10]*a2+p[11]*l)+z*(p[12]+p[13]*a2+p[14]*l)};
	  
	  for(int iff=0;iff<3;iff++)
	    if(iff==1||ith!=0)
	      {
		double contr=pow((Y_fit[ic*3+iff]-Y_teo[iff])/err_Y_fit[ic*3+iff],2);
		ch+=contr;
		if(flag==1)
		  cout<<"contr ic "<<ic<<", iff "<<iff<<", ml="<<ml[iens][iboot]<<", a="<<
		    lat[ibeta[iens]][iboot]<<", ("<<Y_fit[ic*3+iff]<<"-"<<Y_teo[iff]<<")/"<<
		    err_Y_fit[ic*3+iff]<<"="<<contr<<endl;
		ndof++;
	      }
	}
}

bvec fit(boot &chi2_val)
{
  chi2_val=boot(nboot,njack);
  
  //copy X
  Y_fit=new double[3*nens*nth_to_use];
  err_Y_fit=new double[3*nens*nth_to_use];
  
  TMinuit minu;
  minu.SetFCN(chi2);
  
  minu.DefineParameter(0,"FP0_Q0",0.7,0.2,0,0);
  minu.DefineParameter(1,"FP0_A2",0,0.1,0,0);
  minu.DefineParameter(2,"FP0_ML",0,0.1,0,0);
  minu.DefineParameter(3,"FP_Z_Q0",0,0.1,0,0);
  minu.DefineParameter(4,"FP_Z_A2",0,0.1,0,0);
  minu.DefineParameter(5,"FP_Z_ML",0,0.1,0,0);
  minu.DefineParameter(6,"F0_Z_Q0",0,0.1,0,0);
  minu.DefineParameter(7,"F0_Z_A2",0,0.1,0,0);
  minu.DefineParameter(8,"F0_Z_ML",0,0.1,0,0);
  minu.DefineParameter(9,"FT_Q0",0.8,0.1,0,0);
  minu.DefineParameter(10,"FT_A2",0,0.1,0,0);
  minu.DefineParameter(11,"FT_ML",0,0.1,0,0);
  minu.DefineParameter(12,"FT_Z_Q0",0,0.1,0,0);
  minu.DefineParameter(13,"FT_Z_A2",0,0.1,0,0);
  minu.DefineParameter(14,"FT_Z_ML",0,0.1,0,0);
  //minu.FixParameter(1);
  //minu.FixParameter(2);
  //minu.FixParameter(4);
  //minu.FixParameter(5);
  //minu.FixParameter(7);
  //minu.FixParameter(8);
  minu.SetPrintLevel(-1);
  
  int npar=15;
  bvec pars(npar,nboot,njack);
  for(iboot=0;iboot<=nboot;iboot++)
    {
      for(int iens=0;iens<nens;iens++)
        if(use[iens])
	  for(int ith=0;ith<nth_to_use;ith++)
          {
	    int ic=icombo(iens,ith);

	    Y_fit[ic*3+0]=FP[ic][iboot];
	    Y_fit[ic*3+1]=F0[ic][iboot];
	    Y_fit[ic*3+2]=FT[ic][iboot];
            err_Y_fit[ic*3+0]=FP[ic].err();
            err_Y_fit[ic*3+1]=F0[ic].err();
            err_Y_fit[ic*3+2]=FT[ic].err();
          }

      
      if(iboot==nboot) minu.SetPrintLevel(1);
      
      //minimize
      minu.Migrad();
      
      //get back parameters
      double dum;
      double p[npar];
      for(int ipar=0;ipar<npar;ipar++)
	{
	  minu.GetParameter(ipar,pars[ipar][iboot],dum);
	  p[ipar]=pars[ipar][iboot];
	}
      
      chi2(npar,NULL,chi2_val[iboot],p,iboot==nboot);
    }
  
  delete[] Y_fit;
  delete[] err_Y_fit;
  
  return pars;
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
  FP=F0=FT=EP=Q2=Z=bvec(nens*nth_to_use,nboot,njack);
  MD=MP=bvec(nens,nboot,njacks);
  
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
	  
      //convert jack to boot and pass to dimensionful quantity
      jack MD_j(njacks),MP_j(njacks);
      MD_j.load(combine(base_path_MD,ens_path[iens]).c_str(),0);
      MP_j.load(combine(base_path_MP,ens_path[iens]).c_str(),0);
      boot_from_jack(MD[iens],MD_j,iboot_jack);
      boot_from_jack(MP[iens],MP_j,iboot_jack);
      cout<<"a*MD: "<<smart_print(MD[iens])<<", a*MP: "<<smart_print(MP[iens])<<endl;
      MD[iens]/=lat[ibeta[iens]];
      MP[iens]/=lat[ibeta[iens]];
      
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
	  
	  cout<<Q2[ic]<<endl;
	  
	  //compute
	  Z[ic]=fun_Z(MP[iens],MD[iens],Q2[ic]);
	  
	  //remove pole
	  remove_pole(FP[ic],F0[ic],FT[ic],Q2[ic],MD[iens]);
	  //remove_pole(FP[ic],F0[ic],FT[ic],Q2[ic],stranged?MDS_ph:MD_ph);
	  
	  //divide FT by FP
	  if(study_T_fr_P) FT[ic]/=FP[ic];
	}
    }
  fclose(an_input_file);
  
  cout<<MD<<endl;
  
  //define ml and ref ml
  ml=bvec(nens,nboot,njack);
  for(int iens=0;iens<nens;iens++)
    if(use[iens])
      {      
	int b=ibeta[iens];
	//define ml
	cout<<iens<<" "<<b<<" "<<lmass[iens]<<endl;
	ml[iens]=lmass[iens]/lat[b]/Zp[b];
      }
  
  //compute chi2
  boot chi2_val;
  bvec p=fit(chi2_val);
  
  //max Z
  double Z_min=fun_Z(stranged?MK_ph:MP_ph,MD_ph,Q2_max_ph[stranged]);
  double Z_max=0;
  int npoints_plot=100;
  double Z_plot[npoints_plot],dZ=(Z_max-Z_min)/(npoints_plot-1);
  double Q2_plot[npoints_plot];
  bvec y_plot_P(npoints_plot,nboot,njack),y_plot_0(npoints_plot,nboot,njack),y_plot_T(npoints_plot,nboot,njack);
  
  //plot ensemble per ensemble, in Z
  for(int iens=0;iens<nens;iens++)
    {
      if(use[iens])
	for(int ipoint_plot=0;ipoint_plot<npoints_plot;ipoint_plot++)
	  {
	    double z=Z_plot[ipoint_plot]=Z_min+dZ*ipoint_plot;
	    int b=ibeta[iens];
	    boot a2=pow(lat[b],2),l=ml[iens];
	    y_plot_P[ipoint_plot]=(p[0]+p[1]*a2+p[2]*l)+z*(p[3]+p[4]*a2+p[5]*l);
	    y_plot_0[ipoint_plot]=(p[0]+p[1]*a2+p[2]*l)+z*(p[6]+p[7]*a2+p[8]*l);
	    y_plot_T[ipoint_plot]=(p[9]+p[10]*a2+p[11]*l)+z*(p[12]+p[13]*a2+p[14]*l);
	    
	    //return to T
	    if(study_T_fr_P) y_plot_T[ipoint_plot]*=y_plot_P[ipoint_plot];
	  }

      ofstream outP(combine("fit_in_Z/plots/FP_ens_%d.xmg",iens).c_str());
      ofstream out0(combine("fit_in_Z/plots/F0_ens_%d.xmg",iens).c_str());
      ofstream outT(combine("fit_in_Z/plots/FT_ens_%d.xmg",iens).c_str());
      
      if(use[iens])
	{
	  outP<<write_polygon(Z_plot,y_plot_P);
	  out0<<write_polygon(Z_plot,y_plot_0);
	  outT<<write_polygon(Z_plot,y_plot_T);
	  outP<<"&"<<endl;
	  out0<<"&"<<endl;
	  outT<<"&"<<endl;
	  outP<<write_ave_line(Z_plot,y_plot_P);
	  out0<<write_ave_line(Z_plot,y_plot_0);
	  outT<<write_ave_line(Z_plot,y_plot_T);
	  outP<<"&"<<endl;
	  out0<<"&"<<endl;
	  outT<<"&"<<endl;
	  outP<<"@type xydy"<<endl;
	  out0<<"@type xydy"<<endl;
	  outT<<"@type xydy"<<endl;
	  
	  //data
	  for(int ith=0;ith<nth_to_use;ith++)
	    {
	      int ic=icombo(iens,ith);
	      double X=Z[ic].med();
	      if(ith) outP<<X<<" "<<FP[ic]<<endl;
	      out0<<X<<" "<<F0[ic]<<endl;
	      if(ith) outT<<X<<" "<<FT[ic]<<endl;
	    }
	  outP<<"&"<<endl;
	  out0<<"&"<<endl;
	  outT<<"&"<<endl;
	}
    }
      
  //plot
  for(int ipoint_plot=0;ipoint_plot<npoints_plot;ipoint_plot++)
    {
      double z=Z_plot[ipoint_plot]=Z_min+dZ*ipoint_plot;
      boot a2(nboot,njack),l=ml_phys;
      a2=0;
      y_plot_P[ipoint_plot]=(p[0]+p[1]*a2+p[2]*l)+z*(p[3]+p[4]*a2+p[5]*l);
      y_plot_0[ipoint_plot]=(p[0]+p[1]*a2+p[2]*l)+z*(p[6]+p[7]*a2+p[8]*l);
      y_plot_T[ipoint_plot]=(p[9]+p[10]*a2+p[11]*l)+z*(p[12]+p[13]*a2+p[14]*l);
      
      Q2_plot[ipoint_plot]=fun_Zone(stranged?MK_ph:MP_ph,MD_ph,Z_plot[ipoint_plot]);
      
      //return to T
      if(study_T_fr_P) y_plot_T[ipoint_plot]*=y_plot_P[ipoint_plot];

      put_pole(y_plot_P[ipoint_plot],y_plot_0[ipoint_plot],y_plot_T[ipoint_plot],Q2_plot[ipoint_plot],stranged?MDS_ph:MD_ph);
      
      //exp
      if(study_T_fr_P) y_plot_T[ipoint_plot]/=y_plot_P[ipoint_plot];
    }
  
  //restore pole
  for(int iens=0;iens<nens;iens++)
    if(use[iens])
      for(int ith=0;ith<nth_to_use;ith++)
	{
	  int ic=icombo(iens,ith);
	  put_pole(FP[ic],F0[ic],FT[ic],Q2[ic],MD[iens]);
	  //put_pole(FP[ic],F0[ic],FT[ic],Q2[ic],stranged?MDS_ph:MD_ph);
	}
      
  for(int Z_Q2_flag=0;Z_Q2_flag<2;Z_Q2_flag++)
    {
      char Z_Q2_tag[2][4]={"Z","Q2"};
      ofstream outP(combine("fit_in_Z/plots/FP_chir_cont_funz_%s.xmg",Z_Q2_tag[Z_Q2_flag]).c_str());
      ofstream out0(combine("fit_in_Z/plots/F0_chir_cont_funz_%s.xmg",Z_Q2_tag[Z_Q2_flag]).c_str());
      ofstream outT(combine("fit_in_Z/plots/FT_chir_cont_funz_%s.xmg",Z_Q2_tag[Z_Q2_flag]).c_str());
      
      double *K_plot=Z_Q2_flag?Q2_plot:Z_plot;
      
      outP<<write_polygon(K_plot,y_plot_P);
      out0<<write_polygon(K_plot,y_plot_0);
      outT<<write_polygon(K_plot,y_plot_T);
      outP<<"&"<<endl;
      out0<<"&"<<endl;
      outT<<"&"<<endl;
      outP<<write_ave_line(K_plot,y_plot_P);
      out0<<write_ave_line(K_plot,y_plot_0);
      outT<<write_ave_line(K_plot,y_plot_T);
      outP<<"&"<<endl;
      out0<<"&"<<endl;
      outT<<"&"<<endl;
      outP<<"@type xydy"<<endl;
      out0<<"@type xydy"<<endl;
      outT<<"@type xydy"<<endl;

      if(0)
      for(int iens=0;iens<nens;iens++)
	if(use[iens])
	  {
	    for(int ith=0;ith<nth_to_use;ith++)
	      {
		int ic=icombo(iens,ith);
		double X=(Z_Q2_flag?Q2[ic]:Z[ic]).med();
		if(ith) outP<<X<<" "<<FP[ic]<<endl;
		out0<<X<<" "<<F0[ic]<<endl;
		if(ith) outT<<X<<" "<<FT[ic]<<endl;
	      }
	    outP<<"&"<<endl;
	    out0<<"&"<<endl;
	    outT<<"&"<<endl;
	  }
    }
  
  cout<<p<<endl;
  cout<<"chi2: "<<smart_print(chi2_val)<<"/"<<ndof<<"="<<smart_print(chi2_val/ndof)<<endl;
  
  return 0;  
}  
