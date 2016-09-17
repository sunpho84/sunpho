#include "common.cpp"

int njacks=16,nth=5;
int tmin_K,tmax_K;
int tmin_H,tmax_H;
int im_spec,im_S0,im_S1;

const char path[2][60]={"../CORRECTIVE/2pts_P5P5_30_00","../CORRECTIVE/2pts_P5P5_30_30"};

int icombo_2pts(int r,int im1,int im2,int ith2)
{return 0+2*(im1+3*(r+2*(im2+3*ith2)));}

jvec load_2pts(int ism,int im1,int im2,int ith2)
{return (jvec_load(path[ism],T,njack,icombo_2pts(0,im1,im2,ith2))+
	 jvec_load(path[ism],T,njack,icombo_2pts(1,im1,im2,ith2)))/2;}

void read_input(const char *path)
{
  FILE *input_file=open_file(path,"r");

  read_formatted_from_file_expecting((char*)&im_spec,input_file,"%d","im_spec");
  read_formatted_from_file_expecting((char*)&im_S0,input_file,"%d","im_S0");
  read_formatted_from_file_expecting((char*)&im_S1,input_file,"%d","im_S1");

  read_formatted_from_file_expecting((char*)&tmin_K,input_file,"%d","tint_K");
  read_formatted_from_file((char*)&tmax_K,input_file,"%d","tint_K");
  
  read_formatted_from_file_expecting((char*)&tmin_H,input_file,"%d","tint_H");
  read_formatted_from_file((char*)&tmax_H,input_file,"%d","tint_H");
}

int main()
{
  read_set_pars("../data_pars");
  read_input("analysis_pars");
  jvec correctL[4],correctS[4];
  int il_case[4]={0,0,0,1},ic_case[4]={0,1,2,2};
  for(int icase=0;icase<4;icase++)
    {
      int il=il_case[icase];
      int ic=ic_case[icase];
      
      int tmin=(std::max(il,ic)==2)?tmin_H:tmin_K;
      int tmax=(std::max(il,ic)==2)?tmax_H:tmax_K;
      
      jvec M(nth,njacks),ZL(nth,njacks),ZS(nth,njacks);
      for(int ithc=0;ithc<nth;ithc++)
	{
	  jvec sl=load_2pts(0,il,ic,ithc).simmetrized(1);
	  jvec ss=load_2pts(1,il,ic,ithc).simmetrized(1);
	  
	  two_pts_SL_fit(M.data[ithc],ZL.data[ithc],ZS.data[ithc],sl,ss,tmin,tmax,tmin,tmax,
			 combine("../CORRECTIVE/th%d_%d%d_sl.xmg",ithc,il,ic).c_str(),
			 combine("../CORRECTIVE/th%d_%d%d_ss.xmg",ithc,il,ic).c_str());
	}
      
      correctL[icase]=ZL/ZL[0];
      correctS[icase]=ZS/ZS[0];
      
      ofstream out(combine("../CORRECTIVE/o_%d.xmg",icase).c_str());
      out<<"@type xydy"<<endl; 
      out<<correctL[icase]<<endl;
      out<<"&"<<endl;
      out<<correctS[icase]<<endl;
    }
  
  char path_out[]="../CORRECTIVE/corrective_factor";
  correctS[0].write_to_binfile(path_out);
  correctS[1].append_to_binfile(path_out);
  correctS[2].append_to_binfile(path_out);
  correctS[3].append_to_binfile(path_out);

  ofstream out("../CORRECTIVE/prodamir.txt");
  for(int ith=0;ith<nth;ith++)
    {
      out<<endl;
      for(int icase=0;icase<4;icase++)
	{
	  for(int ijack=0;ijack<=njacks;ijack++) out<<correctS[icase][ith][ijack]<<" ";
	  out<<endl;
	}
    }
  
  return 0;
}
