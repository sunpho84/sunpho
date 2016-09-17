#include "common.cpp"

const int RE=0,IM=1;

char data_list_file[1024];

int tmin_P5P5_wl_loc,tmax_P5P5_wl_loc;
int tmin_P5P5_wl_sme,tmax_P5P5_wl_sme;
int tmin_P5P5_wh_loc,tmax_P5P5_wh_loc;
int tmin_P5P5_wh_sme,tmax_P5P5_wh_sme;

int tmin_MEL,tmax_MEL;

int nlight,njack;
int T,TH,TSEP,L;
int iml_un,ist_th;
int nmass,ntheta,ibeta;
char base_path[1024];
double *mass,*theta;
int issm;
int gen_im_S0,gen_im_S1;

int iopp[]={6,5,4,3,2,1,0};

jack Eth(jack M,int ith)
{
  double Qi=M_PI*theta[ith]/L;
  double Q2=3*sqr(Qi);
  //jack E0=sqrt(M*M+Q2);
  
  jack MH=M/2;
  double QH=sqrt(Q2)/2;
  
  jack SMH=sinh(MH);
  double SQH=sin(QH);
  
  jack S2MH=SMH*SMH;
  double S2QH=SQH*SQH;
  
  jack ARG=sqrt(S2MH+S2QH);
  
  jack E=2*asinh(ARG);
  
  return E;
}

void read_pars(const char *input)
{
  FILE *fin=open_file(input,"r");
  read_formatted_from_file_expecting(data_list_file,fin,"%s","data_list_file");
  read_formatted_from_file_expecting((char*)&issm,fin,"%d","is_smeared");
  read_formatted_from_file_expecting((char*)&gen_im_S0,fin,"%d","im_S0");
  read_formatted_from_file_expecting((char*)&gen_im_S1,fin,"%d","im_S1");
  
  if(issm)
    {
      read_formatted_from_file_expecting((char*)&tmin_P5P5_wl_loc,fin,"%d","tmin_P5P5_with_light_loc");
      read_formatted_from_file_expecting((char*)&tmax_P5P5_wl_loc,fin,"%d","tmax_P5P5_with_light_loc");
      read_formatted_from_file_expecting((char*)&tmin_P5P5_wl_sme,fin,"%d","tmin_P5P5_with_light_sme");
      read_formatted_from_file_expecting((char*)&tmax_P5P5_wl_sme,fin,"%d","tmax_P5P5_with_light_sme");
      read_formatted_from_file_expecting((char*)&tmin_P5P5_wh_loc,fin,"%d","tmin_P5P5_with_heavy_loc");
      read_formatted_from_file_expecting((char*)&tmax_P5P5_wh_loc,fin,"%d","tmax_P5P5_with_heavy_loc");
      read_formatted_from_file_expecting((char*)&tmin_P5P5_wh_sme,fin,"%d","tmin_P5P5_with_heavy_sme");
      read_formatted_from_file_expecting((char*)&tmax_P5P5_wh_sme,fin,"%d","tmax_P5P5_with_heavy_sme");
    }
  else
    {
      read_formatted_from_file_expecting((char*)&tmin_P5P5_wl_loc,fin,"%d","tmin_with_light_P5P5");
      read_formatted_from_file_expecting((char*)&tmax_P5P5_wl_loc,fin,"%d","tmax_with_light_P5P5");
      read_formatted_from_file_expecting((char*)&tmin_P5P5_wh_loc,fin,"%d","tmin_with_heavy_P5P5");
      read_formatted_from_file_expecting((char*)&tmax_P5P5_wh_loc,fin,"%d","tmax_with_heavy_P5P5");
    }
  
  read_formatted_from_file_expecting((char*)&tmin_MEL,fin,"%d","tmin_MEL");
  read_formatted_from_file_expecting((char*)&tmax_MEL,fin,"%d","tmax_MEL");
  
  fclose(fin);
}

int icombo_2pts_file(int r1,int im1,int r2,int im2,int ith2,int ri)
{
  return ri+2*(r1+2*(im1+nmass*(r2+2*(im2+nmass*ith2))));
}

int icombo_3pts_file(int im_S0,int ith_S0,int im_S1,int ith_S1,int ri)
{
  return ri+2*(im_S0+nmass*(ith_S0+ntheta*(im_S1+nmass*ith_S1)));
}

int icombo_2pts(int im1,int im2,int ith2)
{
  return im1+nmass*(im2+nmass*ith2);
}

//average ++ and --
jvec load_2pts(int sm,const char *name,int im1,int im2,int ith2,int ri)
{
  string path;
  
  if(issm) path=combine("%s/2pts_30_%02d_%s",base_path,sm,name);
  else     path=combine("%s/2pts_%s",base_path,name);
  
  int ia=icombo_2pts_file(0,im1,0,im2,ith2,ri);
  int ib=icombo_2pts_file(1,im1,1,im2,ith2,ri);
  
  jvec a=jvec_load(path.c_str(),T,njack,ia);
  jvec b=jvec_load(path.c_str(),T,njack,ib);
  
  return (a+b)/2;
}

//average only if in source P5
jvec load_3pts(const char *name,int im_S0,int ith_S0,int im_S1,int ith_S1,int ri)
{
  string path;
  
  if(issm) path=combine("%s/3pts_sp0_30_30_%s",base_path,name);
  else     path=combine("%s/3pts_%s",base_path,name);
  
  int ia=icombo_3pts_file(im_S0,ith_S0,im_S1,ith_S1,ri);
  jvec a=jvec_load(path.c_str(),T,njack,ia);
  
  return a;
}

//flip first half and second half
jvec flip_halves(jvec a)
{
  jvec b(T,njack);
  for(int t=0;t<T;t++)
    {
      int t1=(t<TSEP)?TSEP-t:(T-(t-TSEP))%T;
      b[t]=a[t1];
    }
  
  return b;
}

//average th and -th 
jvec load_P5P5(int sm,int im1,int im2,int ith2,int ri)
{
  return (load_2pts(sm,"P5P5",im1,im2,ith2,ri)+load_2pts(sm,"P5P5",im1,im2,iopp[ith2],ri))/2;
}

//load the four symmetric
void load_P5XP5(jvec &a,jvec &b,jvec &c,jvec &d,const char *name,int im_S0,int ith_S0,int im_S1,int ith_S1,int REIM)
{
  a=load_3pts(name,im_S0,ith_S0,im_S1,ith_S1,REIM);
  b=load_3pts(name,im_S0,iopp[ith_S0],im_S1,iopp[ith_S1],REIM);
  c=flip_halves(load_3pts(name,im_S1,ith_S1,im_S0,ith_S0,REIM));
  d=flip_halves(load_3pts(name,im_S1,iopp[ith_S1],im_S0,iopp[ith_S0],REIM));
}

//average th and -th
jvec load_P5V0P5(int im_S0,int ith_S0,int im_S1,int ith_S1)
{
  jvec a,b,c,d;
  
  load_P5XP5(a,b,c,d,"V0P5",im_S0,ith_S0,im_S1,ith_S1,RE);
  
  jvec e=(a+b+c+d)/4;
  
  if(TSEP==TH) e=(e-e.simmetric())/2;
  
  return e;

}

//average th and -th (minus sign!)
jvec load_P5VKP5(int im_S0,int ith_S0,int im_S1,int ith_S1)
{
  jvec a,b,c,d;
  load_P5XP5(a,b,c,d,"VKP5",im_S0,ith_S0,im_S1,ith_S1,IM);
  
  jvec e=(a+c-b-d)/4;
  
  if(TSEP==TH) e=(e+e.simmetric())/2;
  
  return e;
}

//average th and -th (minus sign!)
jvec load_P5TKP5(int im_S0,int ith_S0,int im_S1,int ith_S1)
{
  jvec a,b,c,d;
  load_P5XP5(a,b,c,d,"TKP5",im_S0,ith_S0,im_S1,ith_S1,IM);
  
  jvec e=(a+d-b-c)/4;
  
  if(TSEP==TH) e=(e-e.simmetric())/2;
  
  return e;
}

void fit_ZM_2pts(jack &M,jack &ZS,int im1,int im2,int ith2)
{
  jvec LOC=load_P5P5(0,im1,im2,ith2,0).simmetrized(1);
  
  int tmin_P5P5_loc,tmax_P5P5_loc;
  int tmin_P5P5_sme,tmax_P5P5_sme;
  
  if(im1<nlight||im2<nlight)
    {
      tmin_P5P5_loc=tmin_P5P5_wl_loc;
      tmax_P5P5_loc=tmax_P5P5_wl_loc;
      tmin_P5P5_sme=tmin_P5P5_wl_sme;
      tmax_P5P5_sme=tmax_P5P5_wl_sme;
    }
  else
    {
      tmin_P5P5_loc=tmin_P5P5_wh_loc;
      tmax_P5P5_loc=tmax_P5P5_wh_loc;
      tmin_P5P5_sme=tmin_P5P5_wh_sme;
      tmax_P5P5_sme=tmax_P5P5_wh_sme;
    }
  
  if(issm)
    {
      jvec SME=load_P5P5(30,im1,im2,ith2,0).simmetrized(1);
      jack ZL;
      two_pts_SL_fit(M,ZL,ZS,LOC,SME,tmin_P5P5_loc,tmax_P5P5_loc,tmin_P5P5_sme,tmax_P5P5_sme,combine("P5P5/%02d_%02d_%d.xmg",im1,im2,ith2).c_str());
    }
  else 
    {
      jack ZL2;
      two_pts_fit(M,ZL2,LOC,tmin_P5P5_loc,tmax_P5P5_loc,combine("P5P5/%02d_%02d_%d.xmg",im1,im2,ith2).c_str());
      ZS=sqrt(ZL2);
    }

}

int main()
{
  read_pars("analysis_pars");
  
  read_ensemble_pars(base_path,T,TSEP,ibeta,nmass,mass,iml_un,nlight,ntheta,theta,ist_th,njack,data_list_file);
  TH=L=T/2;
  
  //we stick to the light-heavy
  int im_spec=0;
  
  //load the 2 points and fits M and Z
  jack M[nmass];
  jack Z[nmass];
  for(int im_S0=0;im_S0<nmass;im_S0++)
    {
      int ith_S0=ist_th;
      fit_ZM_2pts(M[im_S0],Z[im_S0],im_spec,im_S0,ith_S0);
    }
  
  ofstream FT_file("FT.dat");
  ofstream FP_file("FP.dat");
  ofstream FM_file("FM.dat");
  ofstream F0_file("F0.dat");
  FILE *Q2_file=open_file("Q2_combo","w");
  
  //loop over all the theta combination
  for(int ith_S1=0;ith_S1<ntheta;ith_S1++)
    for(int ith_S0=ith_S1;ith_S0<ntheta;ith_S0++)
      {
	//compute energy
	jack E_source=Eth(M[gen_im_S0],ith_S0);
	jack E_sink=Eth(M[gen_im_S1],ith_S1);
	
	//construct time dependance
	jvec time_dep(T,njack);
	for(int t=0;t<T;t++)
	  {
	    jack coef=Z[gen_im_S0]*Z[gen_im_S1]/sqrt(2*E_source*2*E_sink);
	    jack forw,back;
	    if(t<TSEP)
	      {
		forw=exp(-E_source*t-E_sink*(TSEP-t));
		back=exp(-E_source*(T-t)-E_sink*(T-(TSEP-t)));
	      }
	    else
	      {
		forw=exp(-E_source*(T-t)-E_sink*(t-TSEP));
		back=exp(-E_source*t-E_sink*(T-(t-TSEP)));
	      }
	    time_dep[t]=coef*(forw+back);
	  }
	
	//loop over matrix element
	int iV0=0,iVK=1,iTK=2;
	char MEL_name[3][4]={"V0","VK","TK"};
	jack MEL[3];
	jvec MEL_corr[3];
	for(int iMEL=0;iMEL<3;iMEL++)
	  {
	    //load 3pts
	    jvec P5MELP5;
	    
	    switch(iMEL)
	      {
	      case 0:P5MELP5=load_P5V0P5(gen_im_S0,ith_S0,gen_im_S1,ith_S1);break;
	      case 1:P5MELP5=load_P5VKP5(gen_im_S0,ith_S0,gen_im_S1,ith_S1);break;
	      case 2:P5MELP5=load_P5TKP5(gen_im_S0,ith_S0,gen_im_S1,ith_S1);break;
	      }
	    
	    P5MELP5.subset(1,TSEP).print_to_file(combine("P5%sP5/%02d_%d_%02d_%d.xmg",MEL_name[iMEL],gen_im_S0,ith_S0,gen_im_S1,ith_S1).c_str());
	    
	    time_dep.subset(1,TSEP).print_to_file(combine("time_dep/%02d_%d_%02d_%d.xmg",gen_im_S0,ith_S0,gen_im_S1,ith_S1).c_str());
	    
	    //compute remotion of time dependance
	    MEL_corr[iMEL]=P5MELP5/time_dep;
	    MEL[iMEL]=constant_fit(MEL_corr[iMEL].subset(1,TSEP),tmin_MEL,tmax_MEL,combine("%s/%02d_%d_%02d_%d.xmg",MEL_name[iMEL],gen_im_S0,ith_S0,gen_im_S1,ith_S1).c_str());
	  }
	
	//ratios
	jvec VK_fr_TK_corr=MEL_corr[iVK]/MEL_corr[iTK];
	jack VK_fr_TK=constant_fit(VK_fr_TK_corr.subset(1,TSEP),tmin_MEL,tmax_MEL,combine("VK_fr_TK/%02d_%d_%02d_%d.xmg",gen_im_S0,ith_S0,gen_im_S1,ith_S1).c_str());
	jvec V0_fr_TK_corr=MEL_corr[iV0]/MEL_corr[iTK];
	jack V0_fr_TK=constant_fit(V0_fr_TK_corr.subset(1,TSEP),tmin_MEL,tmax_MEL,combine("V0_fr_TK/%02d_%d_%02d_%d.xmg",gen_im_S0,ith_S0,gen_im_S1,ith_S1).c_str());
	
	//////////////////////////// compute_form_factors ///////////////////////
	
	//2*Zv
	MEL[iVK]/=-2*1.69;
	MEL[iV0]/=-2*1.69;
	
	//impulses
	jack M_source=M[gen_im_S0],M_sink=M[gen_im_S1];
	double p_source=M_PI*theta[ith_S0]/L,p_sink=M_PI*theta[ith_S1]/L;
	jack P0=E_source+E_sink,Q0=E_source-E_sink;
	double Pi=p_source+p_sink,Qi=p_source-p_sink;
	jack Q2=Q0*Q0-3*Qi*Qi;
	jack P2=P0*P0-3*Pi*Pi;
	jack DM2=M_source*M_source-M_sink*M_sink;
	
	//TK=(E_source*p_sink-p_source*E_sink)*2FT/(M_source+M_sink)
	jack FT=MEL[iTK]*0.5*(M_source+M_sink)/(E_source*p_sink-p_source*E_sink);
	
	//FP,F0
	jack FM,FP,F0;
        if(ith_S0!=iopp[ith_S1])
          {
            jack DET=P0*Qi-Q0*Pi;
            //calculate (f+) and (f-)*Q2/(M_K^2-M_Pi^2)
            FP=(MEL[iV0]*Qi-Q0*MEL[iVK])/DET;
            FM=(P0*MEL[iVK]-MEL[iV0]*Pi)/DET;
            F0=FP+Q2/DM2*FM;
          }
        else
          {
            //calculate (f+) and (f-)*Q2/(M_K^2-M_Pi^2)
            FM=MEL[iVK]/Qi;
            FP=(MEL[iV0]-Q0*FM)/P0;
            F0=FP+Q2/DM2*FM;
          }
	
	FT_file<<Q2<<" "<<FT<<" "<<Pi<<endl;
	FP_file<<Q2<<" "<<FP<<" "<<Pi<<endl;
	FM_file<<Q2<<" "<<FM<<" "<<Pi<<endl;
	F0_file<<Q2<<" "<<F0<<" "<<Pi<<endl;
	
	fprintf(Q2_file,"th0=%f\t" "th1=%f\t" "Q2=%f\n",theta[ith_S0],theta[ith_S1],Q2.med());
    }

  fclose(Q2_file);
  
  return 0;
}
