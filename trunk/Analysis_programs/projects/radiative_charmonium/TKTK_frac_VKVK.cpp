#include <include.h>

#define REAL 0
#define IMAG 1

int T,TH,L,tsep;
int njack,nthetaS0,nthetaS1;

char base_path[1024];
int ibeta;
double lmass,cmass,theta;

int mode;

int tmin,tmax;

int icombo_2pts(int r1,int r2,int ith1,int reim)
{
  return reim+2*(r1+2*(r2+2*ith1));
}

jvec load_2pts(const char *name,int ith1,int r1,int r2,int reim,const char *sl)
{
  return (
	  jvec_load(combine("%s/2pts_%s_%s",base_path,name,sl).c_str(),T,njack,icombo_2pts(r1,r2,ith1,reim))
	  +jvec_load(combine("%s/2pts_%s_%s",base_path,name,sl).c_str(),T,njack,icombo_2pts(!r1,!r2,ith1,reim))
	  )/2;
}

void read_data_list(const char *path)
{
  FILE *input_file=open_file(path,"r");
  
  read_formatted_from_file_expecting(base_path,input_file,"%s","base_path");
  read_formatted_from_file_expecting((char*)(&T),input_file,"%d","T");
  L=TH=T/2;
  read_formatted_from_file_expecting((char*)(&njack),input_file,"%d","njack");
  read_formatted_from_file_expecting((char*)(&ibeta),input_file,"%d","ibeta");
  
  if(mode==2)
    {
      read_formatted_from_file_expecting((char*)(&nthetaS0),input_file,"%d","nthetaS0");
      read_formatted_from_file_expecting((char*)(&nthetaS1),input_file,"%d","nthetaS1");
    }
  else read_formatted_from_file_expecting((char*)(&nthetaS0),input_file,"%d","ntheta");
  
  read_formatted_from_file_expecting((char*)(&theta),input_file,"%lg","theta");
  read_formatted_from_file_expecting((char*)(&tsep),input_file,"%d","tsep");
  read_formatted_from_file_expecting((char*)(&lmass),input_file,"%lg","lmass");
  read_formatted_from_file_expecting((char*)(&cmass),input_file,"%lg","cmass");
}

void read_input()
{
  char data_list_path[1024];
  FILE *input_file=open_file("input","r");
  
  read_formatted_from_file_expecting(data_list_path,input_file,"%s","data_list_file");
  
  read_formatted_from_file_expecting((char*)(&mode),input_file,"%d","mode");
  read_formatted_from_file_expecting((char*)(&tmin),input_file,"%d","tmin");
  read_formatted_from_file_expecting((char*)(&tmax),input_file,"%d","tmax");
  
  read_data_list(data_list_path);
    
  fclose(input_file);
}

int main()
{
  read_input();
  
  ///////////////////////////// Load two points for standing cc //////////////////////////
  
  //load T
  jvec CT_sl=-load_2pts("TKTK",0, 0,0, 0, "30_00").simmetrized(1);
  jvec CT_ss=-load_2pts("TKTK",0, 0,0, 0, "30_30").simmetrized(1);
  
  //////////////////////////////////// Fit masses and Z for standing cc ////////////////////////////////////////
  
  //compute mass
  jack M,ZL,ZS;
  two_pts_SL_fit(M,ZL,ZS,CT_sl,CT_ss,tmin,tmax,tmin,tmax,"T_SL_fit.xmg","T_SS_fit.xmg");
  
  double Zt[4]={0.73,0.74,0.78,0.82};
  double Zv[4]={0.5816,0.6103,0.6451,0.686};

  jack fT=ZL/M*Zt[ibeta];
  cout<<fT<<endl;
  fT.write_to_binfile("fT");
  
  //load V
  jvec CV_sl=load_2pts("VKVK",0, 0,0, 0, "30_00").simmetrized(1);
  jvec CV_ss=load_2pts("VKVK",0, 0,0, 0, "30_30").simmetrized(1);
  
  
  jvec TfrV=CT_sl/CV_sl/sqrt(CT_ss/CV_ss)*Zt[ibeta]/Zv[ibeta];
  
  jack fT_fr_fV=constant_fit(TfrV,tmin,tmax,"fT_fr_fV.xmg");
  
  cout<<"fT/fV="<<fT_fr_fV<<", fT="<<fT_fr_fV*0.420<<endl;
  
  fT_fr_fV.write_to_binfile("fT_fr_fV");
  
  return 0;
}
