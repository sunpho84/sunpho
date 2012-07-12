#include "common.cpp"

const int RE=0,IM=1;

char data_list_file[1024];

int nlight,njack;
int T,TH,TSEP,L;
int iml_un,ist_th;
int nmass,ntheta,ibeta;
char base_path[1024];
double *mass,*theta;

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
  string path=combine("%s/2pts_30_%02d_%s",base_path,sm,name);
  
  int ia=icombo_2pts_file(0,im1,0,im2,ith2,ri);
  int ib=icombo_2pts_file(1,im1,1,im2,ith2,ri);
  
  jvec a=jvec_load(path.c_str(),T,njack,ia);
  jvec b=jvec_load(path.c_str(),T,njack,ib);
  
  return (a+b)/2;
}

//average only if in source P5
jvec load_3pts(const char *name,int im_S0,int ith_S0,int im_S1,int ith_S1,int ri)
{
  string path=combine("%s/3pts_sp0_30_30_%s",base_path,name);
  
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
  
  return (a+b+c+d)/4;
}

//average th and -th (minus sign!)
jvec load_P5VKP5(int im_S0,int ith_S0,int im_S1,int ith_S1)
{
  jvec a,b,c,d;
  load_P5XP5(a,b,c,d,"VKP5",im_S0,ith_S0,im_S1,ith_S1,IM);

  return (a+c-b-d)/4;
}

//average th and -th (minus sign!)
jvec load_P5TKP5(int im_S0,int ith_S0,int im_S1,int ith_S1)
{
  jvec a,b,c,d;
  load_P5XP5(a,b,c,d,"TKP5",im_S0,ith_S0,im_S1,ith_S1,IM);
  
  return (a+d-b-c)/4;
}

void fit_ZM_2pts(jack &M,jack &ZS,int im1,int im2,int ith2)
{
  jvec LOC=load_P5P5(0,im1,im2,ith2,0).simmetrized(1);
  jvec SME=load_P5P5(30,im1,im2,ith2,0).simmetrized(1);
  jack ZL;
  two_pts_SL_fit(M,ZL,ZS,LOC,SME,12,24,12,19,combine("P5P5/%02d_%02d_%d.xmg",im1,im2,ith2).c_str());
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
  
  //we consider the second to be a pion (if im_spec==0)
  int im_S1=0;
  
  //loop over matrix element
  char MEL_name[3][4]={"V0","VK","TK"};
  for(int iMEL=0;iMEL<3;iMEL++)
    {
      //loop over all the theta combination
      for(int ith_S1=0;ith_S1<ntheta;ith_S1++)
	for(int ith_S0=0;ith_S0<ntheta;ith_S0++)
	  for(int im_S0=0;im_S0<nmass;im_S0++)
	    //for(int im_S1=0;im_S1<4;im_S1++)
	    {
	      //load 3pts
	      jvec P5MELP5;
	      
	      switch(iMEL)
		{
		case 0:P5MELP5=load_P5V0P5(im_S0,ith_S0,im_S1,ith_S1);break;
		case 1:P5MELP5=load_P5VKP5(im_S0,ith_S0,im_S1,ith_S1);break;
		case 2:P5MELP5=load_P5TKP5(im_S0,ith_S0,im_S1,ith_S1);break;
		}
	      
	      P5MELP5.print_to_file(combine("P5%sP5/%02d_%02d_%d_%d.xmg",MEL_name[iMEL],im_S0,im_S1,ith_S0,ith_S1).c_str());
	      
	      //compute energy
	      jack E_source=Eth(M[im_S0],ith_S0);
	      jack E_sink=Eth(M[im_S1],ith_S1);
	      
	      //construct time dependance
	      jvec time_dep(T,njack);
	      for(int t=0;t<T;t++)
		{
		  jack coef=Z[im_S0]*Z[im_S1]/sqrt(2*E_source*2*E_sink);
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
	      time_dep.print_to_file(combine("time_dep/%02d_%02d_%d_%d.xmg",im_S0,im_S1,ith_S0,ith_S1).c_str());
	      
	      //compute remotion of time dependance
	      jvec MEL=P5MELP5/time_dep;
	      jack M=constant_fit(MEL.subset(0,TSEP),9,11,combine("%s/%02d_%02d_%d_%d.xmg",MEL_name[iMEL],im_S0,im_S1,ith_S0,ith_S1).c_str());
	    }
    }
  
  return 0;
}
