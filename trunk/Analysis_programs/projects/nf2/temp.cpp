#include "common.cpp"

int main()
{
  bvec MD_LOW(13,nboot,njack);
  bvec MD_HIGH(13,nboot,njack);
  
  MD_LOW.load("/home/francesco/QCD/LAVORI/NF2/RUN_LOW/ANALYSIS/P5P5/D/c_interpolation/interpolated_M_D",0);
  MD_HIGH.load("/home/francesco/QCD/LAVORI/NF2/RUN_HIGH/ANALYSIS/P5P5/D/c_interpolation/interpolated_M_D",0);
  
  for(int iens=0;iens<13;iens++)
    {
      boot s=MD_LOW[iens]+MD_HIGH[iens];
      boot d=MD_LOW[iens]-MD_HIGH[iens];
      cout<<iens<<" "<<MD_LOW[iens]<<" "<<MD_HIGH[iens]<<" "<<d/s<<endl;
    }
  
  return 0;
}
