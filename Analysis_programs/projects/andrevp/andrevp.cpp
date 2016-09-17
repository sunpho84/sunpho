#include "../../src/include.h"

const int nlevls=2;
const int njacks=88;
const int T=96,TH=T/2,L=48;
const int t0=1;

jvec load(const char *path)
{
  jvec out(T,njacks);

  FILE *fin=open_file(path,"r");
  for(int iconf=0;iconf<njacks;iconf++)
    for(int t=0;t<T;t++)
      read_formatted_from_file((char*)(&(out.data[t][iconf])),fin,"%lg",combine("%s iconf %d t=%d",path,iconf,t).c_str());
  fclose(fin);
  out.clusterize();
  
  return out;
}

int main(int narg,char **arg)
{
  gevp_pars_t g(nlevls,njacks,TH,t0);

  debug_load=false;
  
  jvec data[4];
  const char name[4][3]={"PP","AP","PA","AA"};
  for(int i=0;i<4;i++)
    {
      data[i]=load(combine("%s.txt",name[i]).c_str());
      data[i]/=L*L*L;
      data[i].print_to_file("%s.xmg",name[i]);
    }
  const char path[]="data.dat";
  data[0].write_to_binfile(path);
  for(int i=1;i<4;i++)
    data[i].append_to_binfile(path);
    
  int map[2]={0,1},nlevls_sto=2;
  g.load_raw_data("raw_data.xmg",path,map,nlevls_sto,0);

  ////////////////////////////// finished reading input ///////////////////////////
  
  //resolve gevp
  g.gevp();
  g.check_orthogonality();
  
  //reorder
  g.reorder_eig();
  //g.convert_to_full_eig_ve();
  
  for(int t=0;t<TH;t++)
    {
      cout.precision(16);
      cout<<t<<"  ";
      for(int ilev=0;ilev<nlevls;ilev++) cout<</*smart_print*/(g.eig_va[ilev][t])<<" ";
      cout<<endl;
    }
  cout<<endl;
  
  {
    ofstream out("gevp.xmg");
    out<<"@type xydy"<<endl;
    for(int ilev=0;ilev<nlevls;ilev++)
      {
        jvec m=aperiodic_effective_mass(g.eig_va[ilev]);
        out<<m<<endl;
        out<<"&"<<endl;
      }
  }

  
  return 0;
}
