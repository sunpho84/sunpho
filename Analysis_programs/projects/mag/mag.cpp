#include "include.h"
#include "../nf2/common_pars.cpp"

int T,TH;
int njacks=16;

const double Zt[4]={0.73,0.74,0.78,0.82};

boot boot_from_jack(jack in,int *iboot_jack)
{
  boot out(nboot,njacks);
  
  int njack=in.njack;
  for(int iboot=0;iboot<nboot;iboot++) out.data[iboot]=in.data[iboot_jack[iboot]];
  out.data[nboot]=in.data[njack];
  
  return out;
}

jvec load(const char *name,int r,int ri,int par)
{
  jvec a(T,njacks);
  
  a.load(combine("corrs/2pts_%s_00_00",name).c_str(),ri+2*r);
  
  return a.simmetrized(par);
}

jvec load(const char *name,int par=-1)
{
  return (load(name,0,0,par)+load(name,1,0,par))/2;
}

void load_iboot(int *iboot_jack,const char *path)
{
   
  FILE *fiboot=fopen(path,"r");
  if(fiboot==NULL)
    {
      perror(combine("Error opening file iboot %s",path).c_str());
      exit(1);
    }
  int nr=fread(iboot_jack,sizeof(int),100,fiboot);
  if(nr!=100)
    {
      perror(combine("Error loading iboot data %s",path).c_str());
      exit(1);
    }
  fclose(fiboot);
}

jvec fun(jvec corr)
{
  jvec f(corr);
  
  for(int t=0;t<=TH;t++)
    f[t]*=2*t;
  
  return f;
}

jvec integrate(jvec corr)
{
  jvec summ(TH+1,njacks);
  
  summ[0]=0;
  for(int t=1;t<=TH;t++)
    summ[t]=summ[t-1]+(corr[t]+corr[t-1])/2;
  
  return summ;
}

jack integrate2(jack Z2,jack E,int t)
{
  return Z2*exp(-E*t)*(1+E*t)/(pow(E,3));
  //return Z2*exp(-E*t)*(1+E*t)/E/E/E/2;
}

boot conv(jack in,int ibeta,const char *path,int VT=0)
{
  int iboot_jack[100];
  load_iboot(iboot_jack,path);
  
  boot out=boot_from_jack(in,iboot_jack);
  
  out*=Zt[ibeta]/lat[ibeta];
  if(VT==0) out*=Zv[ibeta];
  else      out*=Zt[ibeta];
  
  return out;
}

boot conv2(jack in,int ibeta,const char *path)
{
  int iboot_jack[100];
  load_iboot(iboot_jack,path);
  
  boot out=boot_from_jack(in,iboot_jack);
  out/=lat[ibeta];
  
  return out;
}

void antisimm_two_pts_fit(jack &E,jack &Z2,jvec corr,int tmin,int tmax,const char *path1=NULL,const char *path2=NULL)
{
  E=constant_fit(effective_mass(corr,TH,-1),tmin,tmax,path1);
  jvec temp(corr.nel,corr.njack);
  int TH=temp.nel-1;
  for(int t=0;t<=TH;t++) temp[t]=corr[t]/exp(-E*TH)/sinh(E*(TH-t))*E;
  Z2=constant_fit(temp,tmin,tmax,path2);

  cout<<"E: "<<E<<endl;
}


int main(int narg,char **arg)
{
  FILE *fin=open_file("analysis_pars","r");
  read_formatted_from_file_expecting((char*)&T,fin,"%d","T");
  TH=T/2;
  
  int ibeta;
  read_formatted_from_file_expecting((char*)&ibeta,fin,"%d","ibeta");
  
  char iboot_path[1024];
  read_formatted_from_file_expecting(iboot_path,fin,"%s","iboot_path");

  int tmin_TKVK,tmax_TKVK,inte_limit;
  read_formatted_from_file_expecting((char*)&tmin_TKVK,fin,"%d","tmin");
  read_formatted_from_file_expecting((char*)&tmax_TKVK,fin,"%d","tmax");
  read_formatted_from_file_expecting((char*)&inte_limit,fin,"%d","inte_limit");

  jvec corr_VKTK=(load("V1T1")+load("V2T2")+load("V3T3"))/3;
  jvec corr_TKVK=-(load("T1V1")+load("T2V2")+load("T3V3"))/3;

  jvec corr_TKTK=(load("T1T1",1)+load("T2T2",1)+load("T3T3",1))/3;
  jvec corr_BKBK=(load("B1B1",1)+load("B2B2",1)+load("B3B3",1))/3;
  corr_TKTK[0]=corr_BKBK[0]=0;
  {
    ofstream out("plots/TBKTBK.xmg");
    out<<"@type xydy"<<endl;
    out<<fun(corr_TKTK)<<"&"<<endl;
    out<<fun(corr_BKBK)<<"&"<<endl;
  }

  jvec corr_P5P5=load("P5P5",1);
  {
    ofstream out("plots/corr_P5P5.xmg");
    for(int ijack=0;ijack<=njacks;ijack++)
      {
	out<<"@type xy"<<endl;
	for(int t=0;t<=TH;t++)
	  out<<t+1.5*ijack/njacks<<" "<<corr_P5P5[t][ijack]<<endl;;
	out<<"&"<<endl;
      }
  }
    
  jack P5P5;
  int tmin[4]={10,12,16,19};
  {
    P5P5=constant_fit(effective_mass(corr_P5P5),tmin[ibeta],TH+1,"plots/P5P5.xmg");
    cout<<"P5P5: "<<smart_print(P5P5)<<endl;
  }
  
  jack E_VKTK(njacks),Z2_VKTK(njacks);
  {    
    ofstream out("plots/VKTK.xmg");
    out<<"@type xydy"<<endl;
    out<<aperiodic_effective_mass(corr_VKTK);
    //out<<corr_VKTK<<endl<<"&"<<endl;
    two_pts_fit(E_VKTK,Z2_VKTK,corr_VKTK,tmin_TKVK,tmax_TKVK,"plots/VKTK_M.xmg","plots/VKTK_Zxmg");
    //jvec test(TH+1,njacks);
    //for(int t=0;t<=TH;t++) test[t]=Z2_VKTK*exp(-t*E_VKTK)/(2*E_VKTK);
    //out<<test<<endl<<"&"<<endl;    
  }

  //plot bare corr
  {
    ofstream out("plots/corr_TKVK.xmg");
    out<<"@type xydy"<<endl;
    out<<fun(corr_TKVK);
  }
  
  jack E_TKVK(njacks),Z2_TKVK(njacks);
  {
    ofstream out("plots/TKVK.xmg");
    out<<"@type xydy"<<endl;
    out<<aperiodic_effective_mass(corr_TKVK);
    two_pts_fit(E_TKVK,Z2_TKVK,corr_TKVK,tmin_TKVK,tmax_TKVK,"plots/TKVK_M.xmg","plots/TKVK_Zxmg");
  }

  {
    ofstream out("plots/KK.xmg");
    out<<"@type xydy"<<endl;
    out<<(corr_TKVK+corr_VKTK)/2;
  }

  jvec sum_TKTK=integrate(corr_TKTK);
  jvec sum_BKBK=integrate(corr_BKBK);
  jvec sum_TKVK=integrate(fun(corr_TKVK));
  jvec sum_VKTK=integrate(fun(corr_VKTK));
  jvec sum_VKTK_2(sum_VKTK);
  jvec sum_TKVK_2(sum_VKTK);
  for(int t=inte_limit+1;t<=TH;t++)
    {
      sum_VKTK_2[t]=sum_VKTK_2[inte_limit]+integrate2(Z2_VKTK,E_VKTK,inte_limit)-integrate2(Z2_VKTK,E_VKTK,t);
      sum_TKVK_2[t]=sum_TKVK_2[inte_limit]+integrate2(Z2_TKVK,E_TKVK,inte_limit)-integrate2(Z2_TKVK,E_TKVK,t);
    }
  
  {
    ofstream inte_VKTK("plots/inte_VKTK.xmg");
    ofstream inte_TKVK("plots/inte_TKVK.xmg");
    ofstream inte_TKTK("plots/inte_TKTK.xmg");
    ofstream inte_BKBK("plots/inte_BKBK.xmg");
    inte_VKTK<<"@type xydy"<<endl;
    inte_TKVK<<"@type xydy"<<endl;
    inte_TKTK<<"@type xydy"<<endl;
    inte_BKBK<<"@type xydy"<<endl;
    inte_TKVK<<sum_TKVK<<"&";
    inte_VKTK<<sum_VKTK<<"&";
    inte_TKTK<<sum_TKTK<<"&";
    inte_BKBK<<sum_BKBK<<"&";
    inte_TKVK<<sum_VKTK_2<<"&";
    inte_VKTK<<sum_TKVK_2<<"&";
  }
  
  init_latpars();
  
  cout<<"TKVK: "<<smart_print(sum_TKVK[TH])<<endl;
  cout<<"VKTK: "<<smart_print(sum_VKTK[TH])<<endl;
  cout<<"TKVK2: "<<smart_print(sum_TKVK_2[TH])<<endl;
  cout<<"VKTK2: "<<smart_print(sum_VKTK_2[TH])<<endl;
  boot MPI=conv2(P5P5,ibeta,iboot_path);
  cout<<"MPi: "<<smart_print(MPI)<<endl;
  
  boot TKVK=conv(sum_TKVK[TH],ibeta,iboot_path);
  boot VKTK=conv(sum_VKTK[TH],ibeta,iboot_path);
  boot TVKVTK=(VKTK+TKVK)*0.5;
  
  boot BKBK=-conv(sum_BKBK[TH],ibeta,iboot_path,1);
  boot TKTK=-conv(sum_TKTK[TH],ibeta,iboot_path,1);
  boot BTKBTK=-conv(sum_TKTK[TH]+sum_BKBK[TH],ibeta,iboot_path,1)/2;
  cout<<TVKVTK<<endl;
  cout<<BTKBTK<<endl;
  
  TKVK.write_to_binfile("TKVK");
  VKTK.write_to_binfile("VKTK");
  sqr(MPI).write_to_binfile("M2PI");

  TVKVTK.write_to_binfile("TVKVTK");
  BTKBTK.write_to_binfile("BTKBTK");
  BKBK.write_to_binfile("BKBK");
  TKTK.write_to_binfile("TKTK");
  
  return 0;
}
