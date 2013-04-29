#include <include.h>


//fit the mass and the matrix element in SS and SL combo
double *c_two_pts_aSL_fit[2],*e_two_pts_aSL_fit[2];
int TH_two_pts_aSL_fit;
int tmin_two_pts_aSL_fit[2];
int tmax_two_pts_aSL_fit[2];

double fun_two_pts_aSL_fit(double Z1,double Z2,double M,double t)
{return Z1*Z2*exp(-M*t)/(2*M);}

void ch2_two_pts_aSL_fit(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  double M=p[0];
  double ZL=p[1];
  double ZS=p[2];
  
  for(int t=tmin_two_pts_aSL_fit[0];t<=min(tmax_two_pts_aSL_fit[0],TH_two_pts_aSL_fit);t++)
    {
      double num=c_two_pts_aSL_fit[0][t];
      double teo=fun_two_pts_aSL_fit(ZL,ZS,M,t);
      double diff=num-teo;
      double err=e_two_pts_aSL_fit[0][t];
      double cont=sqr(diff/err);
      ch+=cont;
      if(flag==3) cout<<"SL, t="<<t<<", diff=("<<num<<"-"<<teo<<")="<<diff<<" err="<<err<<" cont="<<cont<<endl;
    }
  
  for(int t=tmin_two_pts_aSL_fit[1];t<=min(tmax_two_pts_aSL_fit[1],TH_two_pts_aSL_fit);t++)
    {
      double diff=c_two_pts_aSL_fit[1][t]-fun_two_pts_aSL_fit(ZS,ZS,M,t);
      double err=e_two_pts_aSL_fit[1][t];
      double cont=sqr(diff/err);
      ch+=cont;
      if(flag==3) cout<<"SS, t="<<t<<", diff="<<diff<<" err="<<err<<" cont="<<cont<<endl;
    }
}

void two_pts_aSL_fit(jack &M,jack &ZL,jack &ZS,jvec corrSL,jvec corrSS,int tminL,int tmaxL,int tminS,int tmaxS,const char *path1=NULL,const char *path2=NULL)
{
  jvec ecorrSL=aperiodic_effective_mass(corrSL);
  jvec ecorrSS=aperiodic_effective_mass(corrSS);
  
  jack ML=constant_fit(ecorrSL,tminL,tmaxL,NULL);
  //jack MS=constant_fit(ecorrSS,tminS,tmaxS,NULL);
  M=ML;//jack_weighted_average(ML,MS);
  jvec tempSL(corrSS.nel,corrSS.njack),tempSS(corrSS.nel,corrSS.njack);
  int TH=tempSS.nel-1;
  for(int t=0;t<=TH;t++)
    {
      tempSL[t]=corrSL[t]/exp(-M*t)/2;
      tempSS[t]=corrSS[t]/exp(-M*t)/2;
    }
  
  ZS=sqrt(constant_fit(tempSS,tminS,tmaxS,NULL));
  ZL=constant_fit(tempSL,tminL,tmaxL,NULL)/ZS;
  
  TMinuit minu;
  minu.SetPrintLevel(-1);
  minu.SetFCN(ch2_two_pts_aSL_fit);
  minu.DefineParameter(1,"ZL",ZL.med(),0.001,0,0);
  minu.DefineParameter(2,"ZS",ZS.med(),0.001,0,0);
  
  int njack=ML.njack;
  c_two_pts_aSL_fit[0]=new double[TH+1];
  c_two_pts_aSL_fit[1]=new double[TH+1];
  e_two_pts_aSL_fit[0]=new double[TH+1];
  e_two_pts_aSL_fit[1]=new double[TH+1];
  
  TH_two_pts_aSL_fit=TH;
  tmin_two_pts_aSL_fit[0]=tminL;
  tmin_two_pts_aSL_fit[1]=tminS;
  tmax_two_pts_aSL_fit[0]=tmaxL;
  tmax_two_pts_aSL_fit[1]=tmaxS;
  
  for(int iel=0;iel<=TH;iel++)
    {
      e_two_pts_aSL_fit[0][iel]=corrSL[iel].err();
      e_two_pts_aSL_fit[1][iel]=corrSS[iel].err();
    }
  
  for(int ijack_fit=0;ijack_fit<=njack;ijack_fit++)
    {
      minu.DefineParameter(0,"M",M[ijack_fit],0.001,0,0);
      minu.FixParameter(0);
      for(int iel=0;iel<=TH;iel++)
	{
	  c_two_pts_aSL_fit[0][iel]=corrSL[iel][ijack_fit];
	  c_two_pts_aSL_fit[1][iel]=corrSS[iel][ijack_fit];
	}
      minu.Migrad();
      double dum;
      minu.GetParameter(0,M.data[ijack_fit],dum);
      minu.GetParameter(1,ZL.data[ijack_fit],dum);
      minu.GetParameter(2,ZS.data[ijack_fit],dum);
    }
  
  double ch2,grad[3],par[3]={M[njack],ZL[njack],ZS[njack]};
  minu.Eval(3,grad,ch2,par,3);
  cout<<"ML: "<<ML<<", ch2: "<<ch2<<endl;
  
  if(path1!=NULL) write_constant_fit_plot(path1,ecorrSL,M,tminL,tmaxL);
  if(path2!=NULL) write_constant_fit_plot(path2,ecorrSS,M,tminS,tmaxS);
}


#define REAL 0
#define IMAG 1

int T,TH,L,tsep;
int njack;
int nmass,nspec;
int imass,ispec;

char base_path[1024]=".";
char corr_name[1024];

int tmin,tmax;

int icombo_2pts(int r1,int ispec,int r2,int imass,int reim=0)
{return reim+2*(r2+2*(ispec+nspec*(r1+2*imass)));}

jvec load_2pts_LR(const char *name,int ispec,int imass,const char *sl,int LR,int reim=0)
{
  if(ispec>=nspec||ispec<0) crash("check ispec %d in range [0,%d)",ispec,nspec);
  if(imass>=nmass||imass<0) crash("check imass %d in range [0,%d)",imass,nmass);
  
  int r=0;
  
  char LR_tag[2][3]={"L","R"};
  
  return (
	  jvec_load(combine("%s/2pts_%s_%s_%s",base_path,LR_tag[LR],name,sl).c_str(),T,njack,icombo_2pts(r,ispec,r,imass,reim))
	  +jvec_load(combine("%s/2pts_%s_%s_%s",base_path,LR_tag[LR],name,sl).c_str(),T,njack,icombo_2pts(!r,ispec,!r,imass,reim))
	  )/2;
}

jvec load_2pts(const char *name,int ispec,int imass,const char *sl,int reim=0)
{return (load_2pts_LR(name,ispec,imass,sl,0,reim)+load_2pts_LR(name,ispec,imass,sl,1,reim).shifted(-tsep))/2;}

void read_input()
{
  FILE *input_file=open_file("input","r");
  
  read_formatted_from_file_expecting(corr_name,input_file,"%s","corr_name");
  read_formatted_from_file_expecting((char*)(&njack),input_file,"%d","njack");
  read_formatted_from_file_expecting((char*)(&T),input_file,"%d","T");
  L=TH=T/2;
  read_formatted_from_file_expecting((char*)(&tsep),input_file,"%d","TSep");
  
  read_formatted_from_file_expecting((char*)(&nmass),input_file,"%d","nmass");
  read_formatted_from_file_expecting((char*)(&imass),input_file,"%d","imass");
  read_formatted_from_file_expecting((char*)(&nspec),input_file,"%d","nspec");
  read_formatted_from_file_expecting((char*)(&ispec),input_file,"%d","ispec");
  
  read_formatted_from_file_expecting((char*)(&tmin),input_file,"%d","tmin");
  read_formatted_from_file_expecting((char*)(&tmax),input_file,"%d","tmax");
  
  fclose(input_file);
}

int main()
{
  read_input();
  
  //load
  jvec C_25_00=load_2pts_LR(corr_name, ispec,imass, "25_00",0).simmetrized(1);
  jvec C_25_25=load_2pts_LR(corr_name, ispec,imass, "25_25",0).simmetrized(1);

  {
    ofstream corr("corr.xmg");
    corr<<"@type xydy"<<endl;
    corr<<C_25_00<<endl;
    corr<<"&"<<endl;
    corr<<C_25_25<<endl;
    corr<<"&"<<endl;
  }
  
  {
    ofstream eff_mass("eff_mass.xmg");
    eff_mass<<"@type xydy"<<endl;
    eff_mass<<aperiodic_effective_mass(C_25_00).subset(1,20)<<endl;
    eff_mass<<"&"<<endl;
    eff_mass<<aperiodic_effective_mass(C_25_25).subset(1,20)<<endl;
    eff_mass<<"&"<<endl;
  }
  
  int tmin=9,tmax=12;
  jack M,ZL,ZS;
  two_pts_aSL_fit(M,ZL,ZS,C_25_00,C_25_25,tmin,tmax,tmin,tmax,"temp_25_00.xmg","temp_25_25.xmg");
  
  cout<<"M: "<<smart_print(M)<<endl;
  cout<<"ZS: "<<smart_print(ZS)<<endl;
  cout<<"ZL: "<<smart_print(ZL)<<endl;
  
  return 0;
}
