#include <include.h>

#define REAL 0
#define IMAG 1

int T,TH,L,tsep;
int njack;
int nmass,nspec;
int imass,ispec;

char base_path[1024]="DATA";
char corr_name[1024];

int tmin,tmax;

int icombo_otto(int imass,int ir_combo,int reim)
{return reim+2*(ir_combo+8*imass);}

jvec load_otto_ir_combo(const char *name,int imass,const char *sl,int ir_combo,int reim=0)
{
  if(ir_combo>=8||ir_combo<0) crash("check ir_combo %d in range [0,8)",ir_combo);
  if(imass>=nmass||imass<0) crash("check imass %d in range [0,%d)",imass,nmass);
  
  return jvec_load(combine("%s/otto_wL_R_%s_%s",base_path,name,sl).c_str(),T,njack,icombo_otto(imass,ir_combo,reim));
}

jvec load_otto_CONN_DISC(const char *name,int imass,const char *sl,int CD,int reim=0)
{
  jvec otto[8];
  char CD_tag[2][20]={"CONNECTED","DISCONNECTED"};
  
  for(int ir_combo=0;ir_combo<8;ir_combo++)
    otto[ir_combo]=load_otto_ir_combo(combine("%s-%s",CD_tag[CD],name).c_str(),imass,sl,ir_combo,reim);
  
  jvec outA=otto[0]+otto[1]+otto[6]+otto[7];
  jvec outB=otto[2]+otto[3]+otto[4]+otto[5];
  
  for(int t=0;t<tsep;t++) outA[t]+=outB[tsep-t];
  for(int t=tsep;t<T;t++) outA[t]+=outB[(T-t+tsep)%T];
  
  return outA/8;
}

jvec load_otto(const char *name,int imass,const char *sl,int reim=0)
{return load_otto_CONN_DISC(name,imass,sl,0,reim)-load_otto_CONN_DISC(name,imass,sl,1,reim);}


int icombo_2pts(int r1,int ispec,int r2,int imass,int reim=0)
{return reim+2*(r2+2*(ispec+nspec*(r1+2*imass)));}

jvec load_2pts_LR(const char *name,int ispec,int imass,const char *sl,int LR,int rdiff,int reim=0)
{
  if(ispec>=nspec||ispec<0) crash("check ispec %d in range [0,%d)",ispec,nspec);
  if(imass>=nmass||imass<0) crash("check imass %d in range [0,%d)",imass,nmass);
  
  int r1=0,r2=r1+rdiff;
  
  char LR_tag[2][3]={"L","R"};
  
  return (
	  jvec_load(combine("%s/two_points_w%s_%s_%s",base_path,LR_tag[LR],name,sl).c_str(),T,njack,icombo_2pts(r1,ispec,r2,imass,reim))
	  +jvec_load(combine("%s/two_points_w%s_%s_%s",base_path,LR_tag[LR],name,sl).c_str(),T,njack,icombo_2pts(!r1,ispec,!r2,imass,reim))
	  )/2;
}

jvec load_2pts(const char *name,int ispec,int imass,const char *sl,int rdiff,int reim=0)
{return (load_2pts_LR(name,ispec,imass,sl,0,rdiff,reim)+load_2pts_LR(name,ispec,imass,sl,1,rdiff,reim).shifted(-tsep))/2;}

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
  jvec C_25_00_CH=load_2pts(corr_name, ispec,imass, "25_00",0,0).simmetrized(1);
  jvec C_25_25_CH=load_2pts(corr_name, ispec,imass, "25_25",0,0).simmetrized(1);
  jvec C_25_00_NE=load_2pts(corr_name, ispec,imass, "25_00",1,0).simmetrized(1);
  jvec C_25_25_NE=load_2pts(corr_name, ispec,imass, "25_25",1,0).simmetrized(1);

  {
    ofstream eff_mass("eff_mass_CH.xmg");
    eff_mass<<"@type xydy"<<endl;
    eff_mass<<aperiodic_effective_mass(C_25_00_CH).subset(1,20)<<endl;
    eff_mass<<"&"<<endl;
    eff_mass<<aperiodic_effective_mass(C_25_25_CH).subset(1,20)<<endl;
    eff_mass<<"&"<<endl;
  }
  {
    ofstream eff_mass("eff_mass_NE.xmg");
    eff_mass<<"@type xydy"<<endl;
    eff_mass<<aperiodic_effective_mass(C_25_00_NE).subset(1,20)<<endl;
    eff_mass<<"&"<<endl;
    eff_mass<<aperiodic_effective_mass(C_25_25_NE).subset(1,20)<<endl;
    eff_mass<<"&"<<endl;
  }
  
  int tmin=9,tmax=12;
  jack M_CH,ZL_CH,ZS_CH;
  two_pts_SL_fit(M_CH,ZL_CH,ZS_CH,C_25_00_CH,C_25_25_CH,tmin,tmax,tmin,tmax,"temp_25_00_CH.xmg","temp_25_25_CH.xmg");
  jack M_NE,ZL_NE,ZS_NE;
  two_pts_SL_fit(M_NE,ZL_NE,ZS_NE,C_25_00_NE,C_25_25_NE,tmin,tmax,tmin,tmax,"temp_25_00_NE.xmg","temp_25_25_NE.xmg");
  
  cout<<"M: "<<smart_print(M_CH)<<" "<<smart_print(M_NE)<<endl;
  cout<<"ZS: "<<smart_print(ZS_CH)<<" "<<smart_print(ZS_NE)<<endl;
  cout<<"ZL: "<<smart_print(ZL_CH)<<" "<<smart_print(ZL_NE)<<endl;
  
  
  jvec otto_corr(T,njack);
  otto_corr=0;
  char diag[8][3]={"A0","A1","A2","A3","V0","V1","V2","V3"};
  for(int idiag=0;idiag<8;idiag++)
    otto_corr+=2*load_otto(combine("P5%s",diag[idiag]).c_str(),1,"25_25",REAL);
  ofstream otto_file("otto.xmg");
  otto_file<<"@type xydy"<<endl;
  otto_file<<otto_corr;
  
  jvec dT_an(T,njack);
  jvec dT_nu(T,njack);
  for(int t=0;t<=tsep;t++)
    {
      dT_an[t]=ZL_CH*ZL_NE*ZS_CH*ZS_NE*exp(-M_CH*t)*exp(-M_NE*(tsep-t));
      dT_nu[t]=C_25_00_CH[t]*C_25_00_CH[tsep-t];
    }
      
  for(int t=tsep+1;t<T;t++)
    {
      dT_an[t]=ZL_CH*ZL_NE*ZS_CH*ZS_NE*exp(-M_NE*(t-tsep))*exp(-M_CH*(T-t));

      int t1=(t-tsep)%T;
      if(t1>=TH) t1=T-t1;

      int t2=(T-t)%T;
      if(t2>=TH) t2=T-t2;
      
      dT_nu[t]=C_25_00_CH[t1]*C_25_00_CH[t2];
    }  
  
  {
    ofstream dT_file("dT.xmg");
    dT_file<<"@type xydy"<<endl;
    dT_file<<dT_nu.subset(0,tsep+1);
    dT_file<<"&"<<endl;
    dT_file<<dT_an.subset(0,tsep+1);
  }
  
  {
    ofstream otto_file("ratio.xmg");
    otto_file<<"@type xydy"<<endl;
    otto_file<<(otto_corr/dT_nu).subset(0,tsep+1);
    otto_file<<"&"<<endl;
    otto_file<<(otto_corr/dT_an).subset(0,tsep+1);
  }
  
  return 0;
}
