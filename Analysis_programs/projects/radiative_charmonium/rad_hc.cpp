#include <include.h>

#define REAL 0
#define IMAG 1

int T,TH,L,tsep;
double theta,lmass,cmass;
int njack,nthetaS0,nthetaS1;
const double Zv[4]={0.5816,0.6103,0.6451,0.686};
const double Zt[4]={0.73,0.74,0.78,0.82};

char base_path[1024];
int ibeta;

int tminL_B,tmaxL_B;
int tminS_B,tmaxS_B;
int tminL_P,tmaxL_P;
int tminS_P,tmaxS_P;
int tmin_f,tmax_f;

int icombo_2pts(int r1,int r2,int ith1,int reim)
{
  return reim+2*(r1+2*(r2+2*ith1));
}

int icombo_3pts(int ith1,int ith2,int reim)
{
  return reim+2*(ith1+nthetaS1*ith2);
}

jvec load_3pts(const char *name,int ith1,int ith2,int reim)
{return jvec_load(combine("%s/3pts_sp0_%s_30_30",base_path,name).c_str(),T,njack,icombo_3pts(ith1,ith2,reim));}

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
  read_formatted_from_file_expecting((char*)(&nthetaS0),input_file,"%d","nthetaS0");
  read_formatted_from_file_expecting((char*)(&nthetaS1),input_file,"%d","nthetaS1");
  read_formatted_from_file_expecting((char*)(&theta),input_file,"%lg","theta");
  read_formatted_from_file_expecting((char*)(&tsep),input_file,"%d","tsep");
  read_formatted_from_file_expecting((char*)(&lmass),input_file,"%lg","lmass");
  read_formatted_from_file_expecting((char*)(&cmass),input_file,"%lg","cmass");
  
  if(nthetaS0!=nthetaS1)
    {
      fprintf(stderr,"nthetaS0!=nthetaS1\n");
      exit(1);
    }
}

void read_input()
{
  char data_list_path[1024];
  FILE *input_file=open_file("input","r");
  
  read_formatted_from_file_expecting(data_list_path,input_file,"%s","data_list_file");
  
  read_formatted_from_file_expecting((char*)(&tminL_B),input_file,"%d","tminL_B");
  read_formatted_from_file_expecting((char*)(&tmaxL_B),input_file,"%d","tmaxL_B");
  read_formatted_from_file_expecting((char*)(&tminS_B),input_file,"%d","tminS_B");
  read_formatted_from_file_expecting((char*)(&tmaxS_B),input_file,"%d","tmaxS_B");
  
  read_formatted_from_file_expecting((char*)(&tminL_P),input_file,"%d","tminL_P");
  read_formatted_from_file_expecting((char*)(&tmaxL_P),input_file,"%d","tmaxL_P");
  read_formatted_from_file_expecting((char*)(&tminS_P),input_file,"%d","tminS_P");
  read_formatted_from_file_expecting((char*)(&tmaxS_P),input_file,"%d","tmaxS_P");
  
  read_formatted_from_file_expecting((char*)(&tmin_f),input_file,"%d","tmin_f");
  read_formatted_from_file_expecting((char*)(&tmax_f),input_file,"%d","tmax_f");
  
  read_data_list(data_list_path);
    
  fclose(input_file);
}

int main()
{
  read_input();
  
  ///////////////////////////// Load two points for standing D and D* //////////////////////////
  
  //load ss P5P5 for D
  jvec P5P5_ss=load_2pts("P5P5",0, 0,0, 0, "30_30");
  //load sl P5P5 for D
  jvec P5P5_sl=load_2pts("P5P5",0, 0,0, 0, "30_00");
  
  //load ss BKBK for D
  jvec BKBK_ss=load_2pts("BKBK",0, 0,0, 0, "30_30");
  //load sl BKBK for D
  jvec BKBK_sl=load_2pts("BKBK",0, 0,0, 0, "30_00");
  
  //////////////////////////////////// Fit masses and Z for standing D and D* ////////////////////////////////////////
  
  //compute D mass and Z
  jack M_P5,ZL_P5,ZS_P5;
  two_pts_SL_fit(M_P5,ZL_P5,ZS_P5,P5P5_sl.simmetrized(1),P5P5_ss.simmetrized(1),tminL_P,tmaxL_P,tminS_P,tmaxS_P,"MSL_P5.xmg","MSS_P5.xmg");
  cout<<"D mass: "<<M_P5<<", Z: "<<ZL_P5<<endl;
  
  //compute D* mass and Z
  jack M_BK,ZL_BK,ZS_BK;
  two_pts_SL_fit(M_BK,ZL_BK,ZS_BK,BKBK_sl.simmetrized(1),BKBK_ss.simmetrized(1),tminL_B,tmaxL_B,tminS_B,tmaxS_B,"MSL_BK.xmg","MSS_BK.xmg");
  
  //reconstuct moving D mass
  double qi=M_PI*theta/L;
  double q2=3*sqr(qi);
  jack Eth_P5=sqrt(sqr(M_P5)+q2);
  
  //reconstruct numerically and semi-analitically the time dependance of three points
  jvec Dth_DV_td(T,njack);
  for(int t=0;t<T;t++)
    {
      int dtsep=abs(tsep-t);
      
      Dth_DV_td[t]=(ZS_P5*ZS_BK)*
	(exp((-M_BK*t)+(-Eth_P5*dtsep))
	 +exp((-M_BK*(T-t)+(-Eth_P5*dtsep))))/
	(2*Eth_P5*2*M_BK);
    }
  
  //compare time dependance and its simmetric
  if(tsep==T/2)
    {
      ofstream out("Dth_DV_td_simm.xmg");
      out<<"@type xydy"<<endl;
      out<<Dth_DV_td<<endl;
    }
  
  jack Q2_fit=sqr(M_BK-Eth_P5)-q2;
  cout<<"theta optimal: "<<(M_BK*M_BK-M_P5*M_P5)/(2*M_BK)*L/sqrt(3)/M_PI<<endl;
  cout<<"Q2: "<<Q2_fit<<endl;
  
  //////////////////////////////// Calculate three points ////////////////////////////////

  int itheta=(nthetaS1==1 ? 0 : 1);
  
  //load corrs
  jvec P5thVIBJ_pt1=load_3pts("VIBJ_pt1",0,itheta,REAL);
  jvec P5thVIBJ_pt2=load_3pts("VIBJ_pt2",0,itheta,REAL);
  jvec P5thVIBJ_pt3=load_3pts("VIBJ_pt3",0,itheta,REAL);
  
  //compare the different parts
  {
    ofstream out("P5thVIBJ.xmg");
    out<<"@type xydy"<<endl;
    out<<P5thVIBJ_pt1<<endl;
    out<<"&/n@type xydy"<<endl;
    out<<P5thVIBJ_pt2<<endl;
    out<<"&/n@type xydy"<<endl;
    out<<P5thVIBJ_pt3<<endl;
  }
  
  
  ///////////////////////////// Determine the matrix element between D(th=+-) and D* ////////////////////////////
  
  //part 2 is 12 and ciclic
  jvec Dth_V_DV=(P5thVIBJ_pt1-P5thVIBJ_pt2+P5thVIBJ_pt3)/Dth_DV_td;
  
  //compare matrix element and its simmetric
  if(tsep==T/2)
    {
      ofstream out("Dth_V_DV_simm.xmg");
      out<<"@type xydy"<<endl;
      out<<Dth_V_DV<<endl;
      out<<"&\n@type xydy"<<endl;
      out<<-Dth_V_DV.simmetric()<<endl;
      out<<"&\n@type xydy"<<endl;
      out<<Dth_V_DV.simmetrized(-1)<<endl;
    }
  
  //fit matrix element
  jack R=constant_fit(Dth_V_DV.simmetrized(-1),tmin_f,tmax_f,"R.xmg");
  
  //compute the form factor at 0 recoil
  jack V=R/M_BK*Zv[ibeta];
  
  cout<<"F semi-analytical: "<<V<<endl;
  
  M_BK.write_to_binfile("M_BK");
  M_P5.write_to_binfile("M_P5");
  V.write_to_binfile("F");
  
  (ZL_BK/M_BK*Zt[ibeta]).write_to_binfile("ZHC");
  
  return 0;
}
