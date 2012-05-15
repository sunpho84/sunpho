#include "../HH_common.cpp"

int nbeta=4;

double aM0[4]={1.264,1.131,0.942,0.780};
double aM1[4]={1.71,1.42,1.08,0.84};

double eM0=0.001;
double eM1[4]={0.01,0.02,0.01,0.01};

int main(int narg,char **arg)
{
  init_latpars();
  
  bvec M0(4,nboot,njack);
  bvec M1(4,nboot,njack);
  
  ifstream in("input");
  if(!in.good()) crash("opening input");
  
  for(int ibeta=0;ibeta<nbeta;ibeta++)
    {
      string path_m0,path_m1;
      
      in>>path_m0>>path_m1;
      if(!in.good()) crash("Loading ibeta %02d",ibeta);
      
      jack tM0(njack);tM0.load(path_m0.c_str());
      jack tM1(njack);tM1.load(path_m1.c_str());
      
      M0[ibeta].fill_gauss(tM0.med(),tM0.err(),67832854+ibeta);
      M1[ibeta].fill_gauss(tM1.med(),tM1.err(),832854+ibeta);
      
      M0[ibeta]/=lat[ibeta];
      M1[ibeta]/=lat[ibeta];
    }
  
  
  ofstream out("check_fermilab.xmg");
  out<<"@type xydy"<<endl;
  
  for(int ibeta=0;ibeta<nbeta;ibeta++)
    out<<(lat[ibeta]*lat[ibeta]).med()<<" "<<M0[ibeta]<<endl;
  
  
  out<<"&"<<endl;
  
  for(int ibeta=0;ibeta<nbeta;ibeta++)
    out<<(lat[ibeta]*lat[ibeta]).med()<<" "<<M1[ibeta]<<endl;
  
  out.close();

  //Mh_chir_cont.write_to_binfile("results");
  
  return 0;
}
