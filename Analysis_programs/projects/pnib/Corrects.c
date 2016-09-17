#include <stdarg.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const int njacks=15;
const int TH=24;

//return the element (t,ijack)
int ind(int t,int ijack)
{return ijack+(njacks+1)*t;}

//return average jacknife
double ave(double *j)
{return j[njacks];}

//return jacknife error
double err(double *j)
{
  double sx=0,s2x=0;
  
  for(int ijack=0;ijack<njacks;ijack++)
    {
      double x=j[ijack];
      sx+=x;
      s2x+=x*x;
    }
  sx/=njacks;
  s2x/=njacks;
  s2x-=sx*sx;
  
  return sqrt(fabs(s2x)*(njacks-1));
}

//print grace style
void print_to_file(const char *path,double *corr,int nel)
{
  FILE *fout=fopen(path,"w");
  if(fout==NULL)
    {
      fprintf(stderr,"opening %s for writing\n",path);
      exit(0);
    }
  
  fprintf(fout,"@type xydy\n");
  for(int i=0;i<nel;i++)
    {
      double a=ave(corr+i*(njacks+1));
      double e=err(corr+i*(njacks+1));
      if(!isnan(a)&&!isnan(e)) fprintf(fout,"%d %.16lg %.16lg\n",i,a,e);
    }
}

//read a contraction from file
void read_contr(double *data,const char *path)
{
  FILE *fin=fopen(path,"r");
  if(fin==NULL)
    {
      fprintf(stderr,"error opening %s\n",path);
      exit(0);
    }
  
  for(int t=0;t<=TH;t++)
    for(int ijack=0;ijack<=njacks;ijack++)
      {
	int read_t,read_ijack;
	int rc=fscanf(fin,"%d %d %lg",&read_t,&read_ijack,&(data[ind(t,ijack)]));
	//check read
	if(rc!=3)
	  {
	    fprintf(stderr,"reading t=%d, ijack=%d from file %s\n",t,ijack,path);
	    exit(0);
	  }
	//check on t
	if(read_t!=t)
	  {
	    fprintf(stderr,"expecting t=%d, read %d\n",t,read_t);
	    exit(0);
	  }
	//check on ijack
	if(read_ijack!=ijack)
	  {
	    fprintf(stderr,"expecting ijack=%d, read %d\n",ijack,read_ijack);
	    exit(0);
	  }
      }
}

//read the mass insertion
void read_mass_insertion(double *corr)
{
  //allocate direct and exchange
  double *direct=(double*)malloc(sizeof(double)*(TH+1)*(njacks+1));
  double *exchan1=(double*)malloc(sizeof(double)*(TH+1)*(njacks+1));
  double *exchan2=(double*)malloc(sizeof(double)*(TH+1)*(njacks+1));
  double *exchan3=(double*)malloc(sizeof(double)*(TH+1)*(njacks+1));
  
  //read them
  read_contr(direct,"direct_Lls.txt");
  read_contr(exchan1,"exchange_Lls.txt");
  read_contr(exchan2,"exchange_Lsl.txt");
  read_contr(exchan3,"exchange_Sll.txt");
  
  //combine the two
  for(int t=0;t<=TH;t++)
    for(int ijack=0;ijack<=njacks;ijack++)
      {
	int i=ind(t,ijack);
	corr[i]=direct[i]-exchan1[i]+exchan2[i]-exchan3[i];
      }
  
  free(direct);
  free(exchan1);
  free(exchan2);
  free(exchan3);
}

//read the nucleon, combining direct-exchange
void read_nucleon(double *corr)
{
  //allocate direct and exchange
  double *direct=(double*)malloc(sizeof(double)*(TH+1)*(njacks+1));
  double *exchan=(double*)malloc(sizeof(double)*(TH+1)*(njacks+1));
  
  //read them
  read_contr(direct,"direct_Lll.txt");
  read_contr(exchan,"exchange_Lll.txt");
  
  //combine the two
  for(int t=0;t<=TH;t++)
    for(int ijack=0;ijack<=njacks;ijack++)
      {
	int i=ind(t,ijack);
	corr[i]=direct[i]-exchan[i];
      }
  
  free(direct);
  free(exchan);
}

int main()
{
  //load the correlation function
  double *corr=(double*)malloc(sizeof(double)*(njacks+1)*(TH+1));
  read_nucleon(corr);
  
  //compute the effective mass
  double *eff=(double*)malloc(sizeof(double)*(njacks+1)*TH);
  for(int t=0;t<TH;t++)
    for(int ijack=0;ijack<=njacks;ijack++)
      eff[ind(t,ijack)]=log(corr[ind(t,ijack)]/corr[ind(t+1,ijack)]);
  
  print_to_file("nucleon.xmg",corr,TH+1);
  print_to_file("nucleon_eff_mass.xmg",eff,TH);
  
  //load the mass correction
  double *mass_corr=(double*)malloc(sizeof(double)*(njacks+1)*(TH+1));
  read_mass_insertion(mass_corr);
  
  double *mass_corr_ratio=(double*)malloc(sizeof(double)*(njacks+1)*(TH+1));
  for(int t=0;t<=TH;t++)
    for(int ijack=0;ijack<=njacks;ijack++)
      {
	int i=ind(t,ijack);
	mass_corr_ratio[i]=mass_corr[i]/corr[i];
      }
  
  print_to_file("mass_corr.xmg",mass_corr,TH+1);
  print_to_file("mass_corr_ratio.xmg",mass_corr_ratio,TH+1);
  
  double *mass_corr_ratio_der=(double*)malloc(sizeof(double)*(njacks+1)*TH);
  for(int t=0;t<TH;t++)
    for(int ijack=0;ijack<=njacks;ijack++)
      mass_corr_ratio_der[ind(t,ijack)]=mass_corr_ratio[ind(t+1,ijack)]-mass_corr_ratio[ind(t,ijack)];
  
  print_to_file("mass_corr_ratio_der.xmg",mass_corr_ratio_der,TH);
  
  free(corr);
  free(mass_corr);
  free(mass_corr_ratio);
  free(mass_corr_ratio_der);
  free(eff);
  
  return 0;
}
