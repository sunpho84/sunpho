#include <include.h>

int nmass=9;
int njack=15;
int T=96;

int ranges[3][2]={{0,2},{2,5},{5,9}};
double residues[9]={1.e-10,1.e-20,1.e-10,1.e-20,1.e-30,1.e-10,1.e-20,1.e-30,1.e-40};

int icombo(int im1,int im2,int r1,int r2,int ri)
{return ri+2*(r1+2*(im1+nmass*(r2+2*im2)));}

jvec load_combo(int im1,int im2)
{
  int ic1=icombo(im1,im2,0,0,0);
  int ic2=icombo(im1,im2,1,1,0);
  
  return (jvec_load("P5P5_30_00",T,njack,ic1)+
          jvec_load("P5P5_30_00",T,njack,ic2)).simmetrized(1)/2;
}

void write_combos(const char *outpath,int *range1,int *range2)
{
  int iset=0;
  ofstream out(outpath);
  out<<"@type xydy"<<endl;
  for(int i1=range1[0];i1<range1[1];i1++)
    {
      int i2;
      if(range1==range2) i2=i1;
      else i2=range2[0];
      do
	{
	  out<<effective_mass(load_combo(i1,i2));
	  out<<"@s"<<iset<<" legend \""<<residues[i1]<<" "<<residues[i2]<<"\""<<endl;
	  out<<"&"<<endl;
	  i2++;
	  iset++;
	}
      while(i2<range2[1]);
    }
}

int main()
{
  write_combos("pion.xmg",ranges[0],ranges[0]);
  write_combos("eta.xmg",ranges[1],ranges[1]);
  write_combos("bottomonium.xmg",ranges[2],ranges[2]);
  write_combos("B.xmg",ranges[0],ranges[2]);
  write_combos("Bs.xmg",ranges[1],ranges[2]);
  
  return 0;
}
