#include <include.h>

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
  jvec C_05_00=load_2pts_LR(corr_name, ispec,imass, "05_00",0).simmetrized(1);
  jvec C_05_05=load_2pts_LR(corr_name, ispec,imass, "05_05",0).simmetrized(1);

  {
    ofstream corr("corr.xmg");
    corr<<"@type xydy"<<endl;
    corr<<C_05_00<<endl;
    corr<<"&"<<endl;
    corr<<C_05_05<<endl;
    corr<<"&"<<endl;
  }
  
  {
    ofstream eff_mass("eff_mass.xmg");
    eff_mass<<"@type xydy"<<endl;
    eff_mass<<aperiodic_effective_mass(C_05_00).subset(1,20)<<endl;
    eff_mass<<"&"<<endl;
    eff_mass<<aperiodic_effective_mass(C_05_05).subset(1,20)<<endl;
    eff_mass<<"&"<<endl;
  }
  
  cout<<constant_fit(aperiodic_effective_mass(C_05_00).subset(1,20),8,12,"temp_05_00.xmg")<<endl;
  cout<<constant_fit(aperiodic_effective_mass(C_05_05).subset(1,20),8,12,"temp_05_05.xmg")<<endl;
  
  return 0;
}
