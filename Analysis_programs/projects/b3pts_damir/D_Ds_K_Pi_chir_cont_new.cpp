#include "include.h"
#include "../nf2/common_pars.cpp"

int nth_to_use;
int nens,nens_divita=0;
const int nbeta=4,nth_max=8,nens_max=26;
int ibeta[nens_max],use[nens_max];
int stranged;
double ml_ref=0.0050;
bvec ml(nens_max,nboot,njack);
double lmass[nens_max];
int use_also_divita;

const int reciproc_data=1;
const int remove_from_data=1;
const int use_p2=0;

int plot_iboot;

double X_fit[nens_max*nth_max][2];
double Y_fit[nens_max*nth_max][3];
double err_Y_fit[nens_max*nth_max][3];
double ext_ch2_0p,ext_ch2_t;
bvec Q2,FP,F0,FT,MD,MP,EP,P2,X2;
int ref_ml_beta[4]={-1,-1,-1,-1};
const double MD_ph=1.868,MDS_ph=1.9685;
const double MDstar=2.01,MD0star=2.32;
const double MDSstar=2.112,MD0Sstar=2.32;
const double MP_ph=0.135,MP_ref=0.500,MK_ph=0.4944,MK_ref=0.615;
const double Q2_ph_max[2]={sqr(MD_ph-MP_ph),sqr(MD_ph-MK_ph)};
const double P2_ph_at_Q2_zero[2]={sqr((MD_ph*MD_ph-MP_ph*MP_ph)/(2*MD_ph)),sqr((MD_ph*MD_ph-MK_ph*MK_ph)/(2*MD_ph))};
//const double Q2_Pi_ref_max[2]={sqr(MD_ph-MP_ref),sqr(MD_ph-MK_ref)};
//const double P2_Pi_ref_max[2]={sqr(MD_ph-MP_ref),sqr(MD_ph-MK_ref)};

template <class T> T p2_fun(T q2,T m,T M)
{return sqr( (M*M+m*m-q2)/(2*M) )-sqr(m);}

int is_divita_ens(int iens)
{return iens>=nens-nens_divita;}

template <class T1,class T2> void fp0(T1 &fp,T1 &f0,T1 *pars,T2 q2,T2 p2,T2 ml,T2 a2,int add_divita=0)
{
  T2 t2=use_p2?p2:q2;
  T2 mp2=db0.med()*ml;
  
  T1 A_P0=pars[0]+pars[1]*mp2+(pars[2]+add_divita*pars[3])*a2;
  T1 A_P_1=t2*(pars[4]+pars[5]*mp2)+(pars[6]+add_divita*pars[7])*a2*t2;
  T1 A_0_1=t2*(pars[8]+pars[9]*mp2)+(pars[10]+add_divita*pars[11])*a2*t2;
  T1 A_P_2=t2*t2*pars[20];
  T1 A_0_2=t2*t2*pars[21];
  fp=A_P0+A_P_1+A_P_2;
  f0=A_P0+A_0_1+A_0_2;

  if(!remove_from_data)
    {
      T2 factp;
      T2 fact0;
      if(!stranged)
	{
	  factp=(1-q2/sqr(MDstar));
	  fact0=(1-q2/sqr(MD0star));
	}
      else
	{
	  factp=(1-q2/sqr(MDSstar));
	  fact0=(1-q2/sqr(MD0Sstar));
	}
      
      fp/=factp;
      f0/=fact0;
    }
}

template <class T1,class T2> void ft(T1 &fT,T1 *pars,T2 q2,T2 p2,T2 ml,T2 a2,int add_divita=0)
{
  T2 t2=use_p2?p2:q2;
  T2 mp2=db0.med()*ml;
  
  T1 A_T_0=pars[12]+pars[13]*mp2+(pars[14]+add_divita*pars[15])*a2;
  T1 A_T_1=t2*(pars[16]+pars[17]*mp2)+(pars[18]+add_divita*pars[19])*a2*t2;
  T1 A_T_2=t2*t2*pars[22];
  fT=A_T_0+A_T_1+A_T_2;

  if(!remove_from_data)
    {
      T2 factT;
      if(!stranged) factT=(1-q2/sqr(MDstar));
      else factT=(1-q2/sqr(MDSstar));
      
      fT/=factT;
    }
}

template <class T1,class T2> void put_or_remove_pole(T1 &fP,T1 &f0,T1 &fT,T2 q2,int put_remove)
{
  T2 factP,fact0,factT;
  if(!stranged)
    {
      factP=(1-q2/sqr(MDstar));
      factT=(1-q2/sqr(MDstar));
      fact0=(1-q2/sqr(MD0star));
    }
  else
    {
      factP=(1-q2/sqr(MDSstar));
      factT=(1-q2/sqr(MDSstar));
      fact0=(1-q2/sqr(MD0Sstar));
    }
  if(put_remove==1)
    {
      fP*=factP;
      fT*=factT;
      f0*=fact0;
    }
  else
    {
      fP/=factP;
      fT/=factT;
      f0/=fact0;
    }
}

template <class T1,class T2> void put_pole(T1 &fP,T1 &f0,T1 &fT,T2 q2)
{put_or_remove_pole(fP,f0,fT,q2,0);}
template <class T1,class T2> void remove_pole(T1 &fP,T1 &f0,T1 &fT,T2 q2)
{put_or_remove_pole(fP,f0,fT,q2,1);}

int icombo(int iens,int ith)
{return iens*nth_to_use+ith;}

void boot_from_jack(boot &out,jack in,int *iboot_jack)
{
  int nboot=out.nboot;
  int njack=in.njack;
  for(int iboot=0;iboot<nboot;iboot++) out.data[iboot]=in.data[iboot_jack[iboot]];
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

//calculate the chi square
int contr_flag=0;
double chi2(double *p,int npars)
{
  double ch2[3]={0,0,0};
  
  for(int iens=0;iens<nens;iens++)
    if(use[iens])
      for(int ith=0;ith<nth_to_use;ith++)
	{
	  int ic=icombo(iens,ith);
	  
	  double q2=X_fit[ic][0];
	  double p2=P2[ic][plot_iboot];
	  
	  //double MP_max=X_fit[icombo(iens,0)][0];
	  double ml_x=ml[iens][plot_iboot];
	  double a2=sqr(lat[ibeta[iens]][plot_iboot]/lat[1][nboot]);
	  
	  double Y_teo[3];
	  fp0(Y_teo[0],Y_teo[1],p,q2,p2,ml_x,a2,is_divita_ens(iens));
	  ft(Y_teo[2],p,q2,p2,ml_x,a2,is_divita_ens(iens));
	  
	  for(int iff=0;iff<3;iff++)
	    if(iff==1||ith)
	      {
		double contr=pow((Y_fit[ic][iff]-Y_teo[iff])/err_Y_fit[ic][iff],2);
		ch2[iff]+=contr;
		if(contr_flag==1)
		  cout<<"contr ff("<<iff<<") iens("<<iens<<","<<is_divita_ens(iens)<<") ith("<<ith<<"): "<<contr
		      <<" = ("<<Y_fit[ic][iff]<<" - "<<Y_teo[iff]<<") / "<<err_Y_fit[ic][iff]<<endl;
	      }
	}
  
  //export ch2
  ext_ch2_0p=ch2[0]+ch2[1];
  ext_ch2_t=ch2[2];
  
  return ext_ch2_0p+ext_ch2_t;
}

//wrapper for the calculation of the chi2
void chi2_wr(int &npar,double *fuf,double &ch,double *p,int flag)
{ch=chi2(p,npar);}

void fit(bvec &pars)
{
  //copy Y errs
  for(int iens=0;iens<nens;iens++)
    if(use[iens])
      for(int ith=0;ith<nth_to_use;ith++)
        {
          int ic=icombo(iens,ith);
	  if(ith)
	    {
	      err_Y_fit[ic][0]=FP[ith].err();
	      err_Y_fit[ic][2]=FT[ith].err();
	      //cout<<" iens: "<<iens<<" ith: "<<ith<<" "<<Y_fit[ic][2]<<endl;
	    }
	  err_Y_fit[ic][1]=F0[ith].err();
	}
  
  //set minuit function
  TMinuit minu;
  minu.SetFCN(chi2_wr);
  if(plot_iboot!=nboot) minu.SetPrintLevel(-1);
  
  //set pars
  int npars=pars.nel;
  for(int ipar=0;ipar<npars;ipar++) minu.DefineParameter(ipar,combine("P%d",ipar).c_str(),0,0.0001,0,0);
  if(!use_also_divita)
    {
      minu.FixParameter(3);
      minu.FixParameter(7);
      minu.FixParameter(11);
      minu.FixParameter(15);
      minu.FixParameter(19);
    }
  
  //slope in q2 do not show a visible mpi dependance
  if(stranged)
    {
      minu.FixParameter(1);
      minu.FixParameter(13);
    }
  minu.FixParameter(5);
  minu.FixParameter(9);
  minu.FixParameter(17);
  
  //q^4 term not really visible
  minu.FixParameter(20);
  minu.FixParameter(21);
  minu.FixParameter(22);
  
  for(int iboot=0;iboot<=nboot;iboot++)
    {
      plot_iboot=iboot;
      //quiet if not iboot 0
      if(iboot!=nboot) minu.SetPrintLevel(-1);
      else minu.SetPrintLevel(1);
      
      for(int iens=0;iens<nens;iens++)
	if(use[iens])
	  for(int ith=0;ith<nth_to_use;ith++)
	    {
	      int ic=icombo(iens,ith);
	      
	      X_fit[ic][0]=Q2[ic][iboot];
	      X_fit[ic][1]=MP[iens][iboot];
	      
	      Y_fit[ic][0]=FP[ic][iboot];
	      Y_fit[ic][1]=F0[ic][iboot];
	      Y_fit[ic][2]=FT[ic][iboot];
	    }
      
      //minimize
      if(iboot==0)
	{
	  minu.Migrad();
	  minu.mnimpr();
	}
      minu.Migrad();
      
      //get back parameters
      double dum;
      for(int ipar=0;ipar<npars;ipar++)
	minu.GetParameter(ipar,pars.data[ipar].data[iboot],dum);
      
      //make it print separate contributions
      if(iboot==nboot)
	{
	  double t[npars];
	  for(int ipar=0;ipar<npars;ipar++) t[ipar]=pars.data[ipar].data[iboot];
	  contr_flag=1;
	  chi2(t,npars);
	  contr_flag=0;
	}
    }
  
  //print pars
  for(int ipar=0;ipar<npars;ipar++) cout<<"P"<<ipar<<"=("<<pars[ipar]<<")"<<endl;
  cout<<endl;
  
  //count data and print ch2
  int npoints_0p=0,npoints_t=0;
  for(int iens=0;iens<nens;iens++)
    if(use[iens])
      for(int ith=0;ith<nth_to_use;ith++)
	{
	  npoints_0p++;
	  if(ith)
	    {
	      npoints_0p++;
	      npoints_t++;
	    }
	}
  
  cout<<"Chi2_0p = "<<ext_ch2_0p<<" / "<<npoints_0p-2*npars/3<<" = "<<ext_ch2_0p/(npoints_0p-2*npars/3)<<endl;
  cout<<"Chi2_t  = "<<ext_ch2_t<<" / "<<npoints_t-npars/3<<" = "<<ext_ch2_t/(npoints_t-npars/3)<<endl;
}

int main(int narg,char **arg)
{
  init_latpars();
  
  //read ensemble list, meson masses and meson name
  FILE *an_input_file=open_file("analysis_pars","r");
  char base_path_ff[200],base_path_MD[200],base_path_MP[200];
  char base_path_ff_divita[200],base_path_MD_divita[200],base_path_MP_divita[200];
  
  read_formatted_from_file_expecting(base_path_ff,an_input_file,"%s","base_path_ff");
  read_formatted_from_file_expecting(base_path_MD,an_input_file,"%s","base_path_MD");
  read_formatted_from_file_expecting(base_path_MP,an_input_file,"%s","base_path_MP");
  
  read_formatted_from_file_expecting((char*)(&use_also_divita),an_input_file,"%d","use_also_divita");
  if(use_also_divita)
    {
      read_formatted_from_file_expecting(base_path_ff_divita,an_input_file,"%s","base_path_ff_divita");
      read_formatted_from_file_expecting(base_path_MD_divita,an_input_file,"%s","base_path_MD_divita");
      read_formatted_from_file_expecting(base_path_MP_divita,an_input_file,"%s","base_path_MP_divita");
    }
  
  read_formatted_from_file_expecting((char*)(&stranged),an_input_file,"%d","stranged");
  int nth,nth_divita;
  read_formatted_from_file_expecting((char*)(&nth),an_input_file,"%d","nth");
  if(use_also_divita) read_formatted_from_file_expecting((char*)(&nth_divita),an_input_file,"%d","nth_divita");
  read_formatted_from_file_expecting((char*)(&nth_to_use),an_input_file,"%d","nth_to_use");
  read_formatted_from_file_expecting((char*)(&nens),an_input_file,"%d","nens");
  if(use_also_divita)
    {
      read_formatted_from_file_expecting((char*)(&nens_divita),an_input_file,"%d","nens_divita");
      nens+=nens_divita;
    }
  if(nth>nth_max) crash("nth>nth_max");
  if(nth_to_use>nth) crash("nth_to_use>nth");
  if(nens>nens_max) crash("nens>nens_max");
  
  FP=F0=FT=EP=P2=Q2=X2=bvec(nens*nth_to_use,nboot,njack);
  MD=MP=bvec(nens,nboot,njack);
  
  //read data
  ofstream bare_data_table("bare_data_table");
  for(int iens=0;iens<nens;iens++)
    {
      char ens_path[1024];
      read_formatted_from_file((char*)&(use[iens]),an_input_file,"%d","use");
      read_formatted_from_file((char*)&(ibeta[iens]),an_input_file,"%d","ibeta");
      read_formatted_from_file((char*)&(lmass[iens]),an_input_file,"%lg","lmass");
      ml[iens]=lmass[iens]/Zp[ibeta[iens]]/lat[ibeta[iens]];
      read_formatted_from_file(ens_path,an_input_file,"%s","ens_path");
      
      int is_divita=is_divita_ens(iens);
      char *base_path_ff_ptr=is_divita?base_path_ff_divita:base_path_ff;
      
      //read data
      int nth_t=is_divita?nth_divita:nth;
      jvec EP_j(nth_t,njack),Q2_j(nth_t,njack),FP_j(nth_t,njack),F0_j(nth_t,njack),FT_j(nth_t,njack),F0S_j(nth_t,njack);
      EP_j.load(combine(base_path_ff_ptr,ens_path).c_str(),0);
      Q2_j.load(combine(base_path_ff_ptr,ens_path).c_str(),1);
      FP_j.load(combine(base_path_ff_ptr,ens_path).c_str(),2);
      F0_j.load(combine(base_path_ff_ptr,ens_path).c_str(),3);
      FT_j.load(combine(base_path_ff_ptr,ens_path).c_str(),4);
      F0S_j.load(combine(base_path_ff_ptr,ens_path).c_str(),5);

      //set the heavier mass
      int ib=ibeta[iens];
      if(use[iens]) if(ref_ml_beta[ib]==-1||lmass[iens]>lmass[ref_ml_beta[ib]]) ref_ml_beta[ib]=iens;

      //load iboot
      int iboot_jack[100];
      load_iboot(iboot_jack,ens_path);
      
      //convert jack to boot and pass to dimensionful quantity
      jack MD_j(njack),MP_j(njack);
      MD_j.load(combine(is_divita?base_path_MD_divita:base_path_MD,ens_path).c_str(),0);
      MP_j.load(combine(is_divita?base_path_MP_divita:base_path_MP,ens_path).c_str(),0);
      boot_from_jack(MD[iens],MD_j,iboot_jack);
      boot_from_jack(MP[iens],MP_j,iboot_jack);
      MD[iens]/=lat[ib];
      MP[iens]/=lat[ib];
      
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
	  P2[ic]=sqr(EP[ic])-sqr(EP[icombo(iens,0)]);

	  if(remove_from_data) remove_pole(FP[ic],F0[ic],FT[ic],Q2[ic]);
	  if(reciproc_data)
	    {
	      FP[ic]=1/FP[ic];
	      F0[ic]=1/F0[ic];
	      FT[ic]=1/FT[ic];
	    }
	  
	  X2[ic]=use_p2?P2[ic]:Q2[ic];
	  
          bare_data_table<<ens_path<<" "<<EP[ic]<<" "<<Q2[ith]<<" "<<FP[ic]<<" "<<F0[ic]<<" "<<FT[ic]<<" "<<endl;
        }
    }
  
  //perform the fit
  int npars=23;
  bvec pars(npars,nboot,njack);
  fit(pars);
  
  boot fP0_Q20;
  fp0(fP0_Q20,fP0_Q20,pars.data,0.0,P2_ph_at_Q2_zero[stranged],ml_phys_med,0.0);
  if(reciproc_data) fP0_Q20=1/fP0_Q20;
  cout<<"fP,0(Q2=0): "<<fP0_Q20<<endl;
  boot fT_Q20;
  ft(fT_Q20,pars.data,0.0,P2_ph_at_Q2_zero[stranged],ml_phys_med,0.0);
  if(reciproc_data) fT_Q20=1/fT_Q20;
  cout<<"fT(Q2=0): "<<fT_Q20<<endl;
  cout<<"fT/fP(Q2=0): "<<fT_Q20/fP0_Q20<<endl;
  
  //plot of fp and f0 (removed of the pole) ensemble per ensemble 
  for(int iens=0;iens<nens;iens++)
    {
      int is_divita=is_divita_ens(iens);
      ofstream outP0(combine("plots/fP0_ens_%d.xmg",iens).c_str());
      ofstream outT(combine("plots/fT_ens_%d.xmg",iens).c_str());
      ofstream outT_fr_P(combine("plots/fT_fr_P_ens_%d.xmg",iens).c_str());
      ofstream out0_fr_P(combine("plots/f0_fr_P_ens_%d.xmg",iens).c_str());
      double x[100];
      double q2_max=Q2[icombo(iens,0)].med();
      bvec y0(100,nboot,njack),y1(100,nboot,njack),y2(100,nboot,njack),y3(100,nboot,njack),y4(100,nboot,njack);
      
      for(int i=0;i<100;i++)
	{
          double q2=q2_max/99*i;
          double a2=sqr(lat[ibeta[iens]].med()/lat[1].med());
	  double p2=p2_fun(q2,MP[iens].med(),MD[iens].med());
	  fp0(y0[i],y1[i],pars.data,q2,p2,ml[iens].med(),a2,is_divita);
	  ft(y2[i],pars.data,q2,p2,ml[iens].med(),a2,is_divita);
	  x[i]=use_p2?p2:q2;
	  
	  if(!reciproc_data)
	    {
	      y3[i]=y2[i]/y0[i];
	      y4[i]=y1[i]/y0[i];
	    }
	  else
	    {
	      y3[i]=y0[i]/y2[i];
	      y4[i]=y0[i]/y1[i];
	    }
	}
      
      write_polygon(outP0,x,y0,0);outP0<<"&"<<endl<<write_ave_line(x,y0)<<"&"<<endl;
      write_polygon(outP0,x,y1,2);outP0<<"&"<<endl<<write_ave_line(x,y1)<<"&"<<endl;
      write_polygon(outT,x,y2,0);outT<<"&"<<endl<<write_ave_line(x,y2)<<"&"<<endl;
      write_polygon(outT_fr_P,x,y3,0);outT_fr_P<<"&"<<endl<<write_ave_line(x,y3)<<"&"<<endl;
      write_polygon(out0_fr_P,x,y4,0);out0_fr_P<<"&"<<endl<<write_ave_line(x,y4)<<"&"<<endl;
      outP0<<"@type xydy"<<endl;
      outT<<"@type xydy"<<endl;
      outT_fr_P<<"@type xydy"<<endl;
      out0_fr_P<<"@type xydy"<<endl;
      for(int ith=1;ith<nth_to_use;ith++)
	{
	  boot t=FT[icombo(iens,ith)]/FP[icombo(iens,ith)];
	  if(reciproc_data) t=1/t;
	  boot s=F0[icombo(iens,ith)]/FP[icombo(iens,ith)];
	  if(reciproc_data) s=1/s;
	  outP0<<X2[icombo(iens,ith)].med()<<" "<<FP[icombo(iens,ith)]<<endl;
	  outT<<X2[icombo(iens,ith)].med()<<" "<<FT[icombo(iens,ith)]<<endl;
	  outT_fr_P<<X2[icombo(iens,ith)].med()<<" "<<t<<endl;
	  out0_fr_P<<X2[icombo(iens,ith)].med()<<" "<<s<<endl;
	}
      outP0<<"&"<<endl;
      outT<<"&"<<endl;
      outT_fr_P<<"&"<<endl;
      out0_fr_P<<"&"<<endl;
      for(int ith=0;ith<nth_to_use;ith++) outP0<<X2[icombo(iens,ith)].med()<<" "<<F0[icombo(iens,ith)]<<endl;
      outP0<<"&"<<endl;
    }  

  //plot of fp and f0 and fT (removed of the pole) at the ref pion, lattice spacing per lattice spacing
  for(int ib=0;ib<4;ib++)
    if(ref_ml_beta[ib]!=-1)
      {
	int iens=ref_ml_beta[ib];
	
	//interpolate to ref pion
	bvec FP_ref(nth_to_use,nboot,njack);
	bvec F0_ref(nth_to_use,nboot,njack);
	bvec FT_ref(nth_to_use,nboot,njack);
	bvec FT_fr_P_ref(nth_to_use,nboot,njack);
	bvec Q2_ref(nth_to_use,nboot,njack);
	bvec P2_ref(nth_to_use,nboot,njack);
	bvec X2_ref(nth_to_use,nboot,njack);

	for(int ith=0;ith<nth_to_use;ith++)
	  {
	    int ic=icombo(iens,ith);
	    Q2_ref[ith]=Q2[ic]-sqr(MD[iens]-EP[ic]);//this is -Pi_p2
	    Q2_ref[ith]+=sqr(MD[iens]-sqrt(sqr(stranged?MK_ref:MP_ref)-Q2_ref[ith])); //add again Q0^2
	    P2_ref[ith]=p2_fun(Q2_ref[ith],MD[iens]*0+(stranged?MK_ref:MP_ref),MD[iens]);

	    X2_ref[ith]=use_p2?P2_ref[ith].med():Q2_ref[ith].med();
	    
	    boot FP_teo(nboot,njack),F0_teo(nboot,njack),FT_teo(nboot,njack);
	    fp0(FP_teo,F0_teo,pars.data,Q2[ic],P2[ic],Q2[ic]*0+ml_ref,sqr(lat[ib]));
	    ft(FT_teo,pars.data,Q2[ic],P2[ic],Q2[ic]*0+ml_ref,sqr(lat[ib]));
	    fp0(FP_ref[ith],F0_ref[ith],pars.data,Q2_ref[ith],P2_ref[ith],ml[iens],sqr(lat[ib]));
	    ft(FT_ref[ith],pars.data,Q2_ref[ith],P2_ref[ith],ml[iens],sqr(lat[ib]));
	    
	    FP_ref[ith]*=FP[ic]/FP_teo;
	    F0_ref[ith]*=F0[ic]/F0_teo;
	    FT_ref[ith]*=FT[ic]/FT_teo;
	    
	    if(!reciproc_data) FT_fr_P_ref[ith]=FT_ref[ith]/FP_ref[ith];
	    else               FT_fr_P_ref[ith]=FP_ref[ith]/FT_ref[ith];
	  }
	
	ofstream outP0(combine("plots/fP0_Pi_ref_beta_%d.xmg",ib).c_str());
	ofstream outT(combine("plots/fT_Pi_ref_beta_%d.xmg",ib).c_str());
	ofstream outT_fr_P(combine("plots/fT_fr_P_Pi_ref_beta_%d.xmg",ib).c_str());
	double x[100];
	double q2_max=Q2[icombo(iens,0)].med();
	bvec y0(100,nboot,njack),y1(100,nboot,njack),y2(100,nboot,njack),y3(100,nboot,njack);
	for(int i=0;i<100;i++)
	  {
	    double q2=q2_max/99*i;
	    double p2=p2_fun(q2,stranged?MK_ref:MP_ref,MD[iens].med());
	    
	    double a2=sqr(lat[ib].med()/lat[1].med());
	    x[i]=use_p2?p2:q2;

	    fp0(y0[i],y1[i],pars.data,q2,p2,ml_ref,a2);
	    ft(y2[i],pars.data,q2,p2,ml_ref,a2);
	
	    if(!reciproc_data) y3[i]=y2[i]/y0[i];
	    else               y3[i]=y0[i]/y2[i];
	  }
	
	write_polygon(outP0,x,y0,0);outP0<<"&"<<endl<<write_ave_line(x,y0)<<"&"<<endl;
	write_polygon(outP0,x,y1,2);outP0<<"&"<<endl<<write_ave_line(x,y1)<<"&"<<endl;
	write_polygon(outT,x,y2,0);outT<<"&"<<endl<<write_ave_line(x,y2)<<"&"<<endl;
	write_polygon(outT_fr_P,x,y3,0);outT_fr_P<<"&"<<endl<<write_ave_line(x,y3)<<"&"<<endl;
	
	outP0<<"@type xydy"<<endl;
	outT<<"@type xydy"<<endl;
	outT_fr_P<<"@type xydy"<<endl;
	for(int ith=1;ith<nth_to_use;ith++)
	  {
	    outP0<<X2_ref[ith].med()<<" "<<FP_ref[ith]<<endl;
	    outT<<X2_ref[ith].med()<<" "<<FT_ref[ith]<<endl;
	    outT_fr_P<<X2_ref[ith].med()<<" "<<FT_fr_P_ref[ith]<<endl;
	  }
	outP0<<"&"<<endl;
	outT<<"&"<<endl;
	outT_fr_P<<"&"<<endl;
	for(int ith=0;ith<nth_to_use;ith++) outP0<<X2_ref[ith].med()<<" "<<F0_ref[ith]<<endl;
	outP0<<"&"<<endl;
      }
  
  //plot of fp and f0 with pole put back in place
  {	
    ofstream outP0("plots/fP0_chir_cont.xmg");
    ofstream outT("plots/fT_chir_cont.xmg");
    ofstream outTP("plots/fT_fr_P_chir_cont.xmg");
    double x[100];
    bvec y0(100,nboot,njack),y1(100,nboot,njack),y2(100,nboot,njack),y3(100,nboot,njack);
    for(int i=0;i<100;i++)
      {
	double q2=Q2_ph_max[stranged]/99*i;
	double p2=p2_fun(q2,stranged?MK_ph:MP_ph,MD_ph);
	
	double a2=0;
	x[i]=use_p2?p2:q2;

	fp0(y0[i],y1[i],pars.data,q2,p2,ml_phys_med,a2);
	ft(y2[i],pars.data,q2,p2,ml_phys_med,a2);
	if(reciproc_data)
	  {
	    y0[i]=1/y0[i];
	    y1[i]=1/y1[i];
	    y2[i]=1/y2[i];
	  }
	if(remove_from_data) put_pole(y0[i],y1[i],y2[i],q2);
	
	y3[i]=y2[i]/y0[i];
      }
	
    write_polygon(outP0,x,y0,0);outP0<<"&"<<endl<<write_ave_line(x,y0)<<"&"<<endl;
    write_polygon(outP0,x,y1,2);outP0<<"&"<<endl<<write_ave_line(x,y1)<<"&"<<endl;
    write_polygon(outT,x,y2,0);outT<<"&"<<endl<<write_ave_line(x,y2)<<"&"<<endl;
    write_polygon(outTP,x,y3,0);outTP<<"&"<<endl<<write_ave_line(x,y3)<<"&"<<endl;
  }
  
  //plot f0(q2mx)
  {	
    ofstream out0("plots/f0_q2max.xmg");
    out0<<"@type xydy"<<endl;
    for(int ib=0;ib<4;ib++)
      {
	for(int iens=0;iens<nens;iens++)
	  if(ibeta[iens]==ib)
	    {
	      int ith=0;
	      int ic=icombo(iens,ith);
	      out0<<sqr(ml[iens]).med()<<" "<<F0[ic]<<endl;
	    }
	out0<<"&"<<endl;
      }
  }
  
  //plot M2Pi and MD
  {	
    ofstream out_Pi("plots/M2Pi.xmg");
    ofstream out_D("plots/MD.xmg");
    out_Pi<<"@type xydy"<<endl;
    out_D<<"@type xydy"<<endl;
    for(int ib=0;ib<4;ib++)
      {
	for(int iens=0;iens<nens;iens++)
	  if(ibeta[iens]==ib) 
	    {
	      out_Pi<<1000*ml[iens].med()<<" "<<sqr(MP[iens])<<endl;
	      out_D<<ml[iens].med()<<" "<<MD[iens]<<endl;
	    }
	out_Pi<<"&"<<endl;
	out_D<<"&"<<endl;
      }
  }
  
  return 0;  
}
