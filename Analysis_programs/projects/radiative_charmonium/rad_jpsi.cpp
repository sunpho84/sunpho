#include <include.h>
#include "radiative_common.cpp"

#define REAL 0
#define IMAG 1

int T,TH,L,tsep;
double theta,lmass,cmass;
int nr,ntheta,include_local;
const double Zvm[4]={0.5816,0.6103,0.6451,0.686};

char base_path[1024];
int ibeta;

int tminL_V,tmaxL_V;
int tminS_V,tmaxS_V;
int tminL_P,tmaxL_P;
int tminS_P,tmaxS_P;
int tmin_f,tmax_f;

double Zam[4]={0.746,0.746,0.772,0.780};
//double Za_err[4]={0.011,0.006,0.006,0.006};


int icombo_2pts(int r1,int r2,int ith1,int reim)
{
  return reim+nr*(r1+nr*(r2+2*ith1));
}

int icombo_3pts(int ith1,int ith2,int reim)
{
  return reim+2*(ith1+ntheta*ith2);
}

jvec load_3pts(const char *name,int ith1,int ith2,int reim)
{return jvec_load(combine("%s/3pts_sp0_%s_30_30",base_path,name).c_str(),T,njack,icombo_3pts(ith1,ith2,reim));}

jvec load_2pts(const char *name,int ith1,int r1,int r2,int reim,const char *sl)
{
  if(nr==2)
    return (
	    jvec_load(combine("%s/2pts_%s_%s",base_path,name,sl).c_str(),T,njack,icombo_2pts(r1,r2,ith1,reim))
	    +jvec_load(combine("%s/2pts_%s_%s",base_path,name,sl).c_str(),T,njack,icombo_2pts(!r1,!r2,ith1,reim))
	    )/2;
  else
    return jvec_load(combine("%s/2pts_%s_%s",base_path,name,sl).c_str(),T,njack,icombo_2pts(0,0,ith1,reim));
}

void read_data_list(const char *path)
{
  FILE *input_file=open_file(path,"r");
  
  read_formatted_from_file_expecting(base_path,input_file,"%s","base_path");
  read_formatted_from_file_expecting((char*)(&T),input_file,"%d","T");
  L=TH=T/2;
  read_formatted_from_file_expecting((char*)(&njack),input_file,"%d","njack");
  read_formatted_from_file_expecting((char*)(&ibeta),input_file,"%d","ibeta");
  read_formatted_from_file_expecting((char*)(&nr),input_file,"%d","nr");
  read_formatted_from_file_expecting((char*)(&include_local),input_file,"%d","include_local");
  read_formatted_from_file_expecting((char*)(&ntheta),input_file,"%d","ntheta");
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
  
  read_formatted_from_file_expecting((char*)(&tminL_V),input_file,"%d","tminL_V");
  read_formatted_from_file_expecting((char*)(&tmaxL_V),input_file,"%d","tmaxL_V");
  read_formatted_from_file_expecting((char*)(&tminS_V),input_file,"%d","tminS_V");
  read_formatted_from_file_expecting((char*)(&tmaxS_V),input_file,"%d","tmaxS_V");
  
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
  jvec P5P5_sl;
  if(include_local) P5P5_sl=load_2pts("P5P5",0, 0,0, 0, "30_00");
  
  //load ss P5P5 for D
  jvec P5P5_mov_ss;
  if(ntheta==2) P5P5_mov_ss=load_2pts("P5P5",1, 0,0, 0, "30_30");
  //load sl P5P5 for D
  jvec P5P5_mov_sl;
  if(ntheta==2 && include_local) P5P5_mov_sl=load_2pts("P5P5",1, 0,0, 0, "30_00");
  
  //load ss VKVK for D
  jvec VKVK_ss=load_2pts("VKVK",0, 0,0, 0, "30_30");
  //load sl VKVK for D
  jvec VKVK_sl;
  if(include_local) VKVK_sl=load_2pts("VKVK",0, 0,0, 0, "30_00");

  //////////////////////////////////// Fit masses and Z for standing D and D* ////////////////////////////////////////
  
  //compute D mass and Z
  jack M_P5,ZL_P5,ZS_P5;
  if(include_local)
    {
      two_pts_SL_fit(M_P5,ZL_P5,ZS_P5,P5P5_sl.simmetrized(1),P5P5_ss.simmetrized(1),tminL_P,tmaxL_P,tminS_P,tmaxS_P,
		     "MSL_P5.xmg","MSS_P5.xmg");
      cout<<"PSEUDO D mass: "<<M_P5<<", Z: "<<ZL_P5<<endl;
    }
  else
    {
      two_pts_fit(M_P5,ZS_P5,P5P5_ss.simmetrized(1),tminS_P,tmaxS_P,"MSS_P5.xmg");
      cout<<"PSEUDO D mass: "<<M_P5<<endl;
    }

  //compute D* mass and Z
  jack M_VK,ZL_VK,ZS_VK;
  if(include_local)
    {
      two_pts_SL_fit(M_VK,ZL_VK,ZS_VK,VKVK_sl.simmetrized(1),VKVK_ss.simmetrized(1),tminL_V,tmaxL_V,tminS_V,tmaxS_V,
		     "MSL_VK.xmg","MSS_VK.xmg");
      cout<<"VECT D mass: "<<M_VK<<", Z: "<<ZL_VK<<endl;
    }
  else
    {
      two_pts_fit(M_VK,ZS_VK,VKVK_ss.simmetrized(1),tminS_V,tmaxS_V,"MSS_VK.xmg");
      cout<<"VECT D mass: "<<M_VK<<endl;
    }
  
  //reconstruct moving D mass
  double qi=M_PI*theta/L;
  double q2=3*sqr(qi);
  jack Eth_P5=sqrt(sqr(M_P5)+q2);
  
  //reconstruct numerically and semi-analitically the time dependance of three points
  jvec Dth_DV_td_nu(T,njack);
  jvec Dth_DV_td_sa(T,njack);
  for(int t=0;t<T;t++)
    {
      int dtsep=abs(tsep-t);
      
      if(ntheta==2)
	{	
	  if(t<TH) Dth_DV_td_nu[t]=VKVK_ss[t]*P5P5_mov_ss[dtsep]/(ZS_P5*ZS_VK);
	  else     Dth_DV_td_nu[t]=VKVK_ss[T-t]*P5P5_mov_ss[dtsep]/(ZS_P5*ZS_VK);
	}
      
      Dth_DV_td_sa[t]=(ZS_P5*ZS_VK)*
	(exp((-M_VK*t)+(-Eth_P5*dtsep))
	 +exp((-M_VK*(T-t)+(-Eth_P5*dtsep))))/
	(2*Eth_P5*2*M_VK);
    }
  
  //compare time dependance and its simmetric
  if(tsep==T/2)
    {
      ofstream out("Dth_DV_td_sa_simm.xmg");
      out<<"@type xydy"<<endl;
      out<<Dth_DV_td_sa<<endl;
      out<<"&\n@type xydy"<<endl;
      out<<Dth_DV_td_nu<<endl;
    }
  
  jack Q2_fit=sqr(M_VK-Eth_P5)-q2;
  cout<<"theta optimal: "<<(M_VK*M_VK-M_P5*M_P5)/(2*M_VK)*L/sqrt(3)/M_PI<<endl;
  cout<<"Q2: "<<Q2_fit<<endl;
  
  //////////////////////////////// Calculate three points ////////////////////////////////

  //load corrs
  int ith;
  if(ntheta==2) ith=1;
  else ith=0;
  jvec P5thVIVJ_pt1= load_3pts("VIVJ_pt1",0,ith,IMAG);
  jvec P5thVIVJ_pt2=-load_3pts("VIVJ_pt2",0,ith,IMAG);
  jvec P5thVIVJ=(P5thVIVJ_pt1+P5thVIVJ_pt2)/2;

  //compare the different parts
  {
    ofstream out("P5thVIVJ_parts.xmg");
    out<<"@type xydy"<<endl;
    out<<P5thVIVJ<<endl;
  }

  ///////////////////////////// Determine the matrix element between D(th=+-) and D* ////////////////////////////
  
  jvec Dth_V_DV_sa=P5thVIVJ/Dth_DV_td_sa;
  jvec Dth_V_DV_nu=P5thVIVJ/Dth_DV_td_nu;
  
  //compare matrix element and its simmetric
  if(tsep==T/2)
    {
      ofstream out("Dth_V_DV_sa_simm.xmg");
      out<<"@type xydy"<<endl;
      out<<Dth_V_DV_sa<<endl;
      out<<"&\n@type xydy"<<endl;
      out<<-Dth_V_DV_sa.simmetric()<<endl;
      out<<"&\n@type xydy"<<endl;
      out<<Dth_V_DV_sa.simmetrized(-1)<<endl;
    }
  
  //load M1
  jack M1_P5(njack),M1_VK(njack);
  bool fit_M1=file_exists("M1_P5");
  
  if(fit_M1)
    {
      M1_P5.load("M1_P5",0);
      M1_VK.load("M1_VK",0);
      cout<<"D  M,M1 mass: "<<M_P5<<" "<<M1_P5<<endl;
      cout<<"D* M,M1 mass: "<<M_VK<<" "<<M1_VK<<endl;
    }
  
  //fit matrix element
  jack R_sa=constant_fit(Dth_V_DV_sa.simmetrized(-1),tmin_f,tmax_f,"R_sa.xmg");
  jack R_nu=constant_fit(Dth_V_DV_nu.simmetrized(-1),tmin_f,tmax_f,"R_nu.xmg");
  
  //compute the form factor at 0 recoil
  jack V_sa=R_sa/qi*(M_VK+M_P5)/(2*M_VK)*Zvm[ibeta];
  jack V1_sa;if(fit_M1) V1_sa=R_sa/qi*(M1_VK+M1_P5)/(2*M1_VK)*Zvm[ibeta];
  jack V_nu=R_nu/qi*(M_VK+M_P5)/(2*M_VK)*Zvm[ibeta];
  jack V1_nu;if(fit_M1) V1_nu=R_nu/qi*(M1_VK+M1_P5)/(2*M1_VK)*Zvm[ibeta];
  
  cout<<"F semi-analytical: "<<V_sa<<endl;
  if(fit_M1) cout<<"F1 semi-analytical: "<<V1_sa<<endl;
  cout<<"F numerical: "<<V_nu<<endl;
  if(fit_M1) cout<<"F1 numerical: "<<V1_nu<<endl;
  
  M_VK.write_to_binfile("M_VK");
  M_P5.write_to_binfile("M_P5");
  V_sa.write_to_binfile("V_sa");
  if(fit_M1) V1_sa.write_to_binfile("V1_sa");
  if(fit_M1) V1_nu.write_to_binfile("V1_nu");
  if(ntheta==2) V_nu.write_to_binfile("V_nu");
  
  if(include_local)
    {
      (ZL_VK/M_VK*Za_med[ibeta]).write_to_binfile("ZJPSI");
      (2*cmass*ZL_P5/M_P5/sinh(M_P5)).write_to_binfile("ZETAC");
    }      

  return 0;
}
