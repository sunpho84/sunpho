#include <include.h>
#include <iostream>
#include <sstream>

using namespace std;

int nmass,nlights,iml_un;
int njack=16;
double *mass;
jvec c[2];
int T,TH;
char base_path[1024];
int tSEP;
int ibeta;
double afm[4]={0.0977,0.0847,0.0671,0.0536};

double fun_fit(double Z1,double Z2,double M,int t)
{
  return Z1*Z2*exp(-M*TH)*cosh(M*(TH-t))/M;
}

int icombo(int im1,int im2,int r1,int r2,int iwall,int ri)
{
  if(im1>=nlights)
    {
      cerr<<"Error, requiring uncomputed combo!"<<endl;
      exit(1);
    }
  
  int ic=ri+2*(iwall+3*(r1+2*(im1+nlights*(r2+2*im2))));

  return ic;
}

jvec load_corr(char *base_path,int sm_lev,int im1,int im2)
{
  char path[1024];
  int tag[2]={0,30};
  sprintf(path,"%s/P5P5_30_%02d",base_path,tag[sm_lev]);

  jvec l=(jvec_load(path,T,njack,icombo(im1,im2,0,0,0,0))+
	  jvec_load(path,T,njack,icombo(im1,im2,1,1,0,0)))/2;
  
  jvec r=(jvec_load(path,T,njack,icombo(im1,im2,0,0,1,0))+
	  jvec_load(path,T,njack,icombo(im1,im2,1,1,1,0)))/2;
  
  jvec r2=(jvec_load(path,T,njack,icombo(im1,im2,0,0,2,0))+
	   jvec_load(path,T,njack,icombo(im1,im2,1,1,2,0)))/2;
  
  
  for(int tL=0;tL<T;tL++)
    {
      int tR=tL+16;
      if(tR>=T) tR-=T;
      
      int tR2=tL+24;
      if(tR2>=T) tR2-=T;
      
      l.data[tL]+=r.data[tR]+r2.data[tR2];
    }
  
  
  return l/3;
}

void read_ensemble_pars(char *base_path,int &T,int &ibeta,int &nmass,double *&mass,int &nlights,const char *data_list_file)
{
  FILE *input=open_file(data_list_file,"r");
  
  read_formatted_from_file_expecting(base_path,input,"%s","base_path");
  read_formatted_from_file_expecting((char*)&T,input,"%d","T");
  read_formatted_from_file_expecting((char*)&ibeta,input,"%d","Beta");
  
  read_formatted_from_file_expecting((char*)&nmass,input,"%d","nmass");

  expect_string_from_file(input,"mass_list");
  mass=(double*)malloc(sizeof(double)*nmass);
  for(int imass=0;imass<nmass;imass++) read_formatted_from_file((char*)&(mass[imass]),input,"%lg","mass");
  read_formatted_from_file_expecting((char*)&iml_un,input,"%d","iml_un");
  read_formatted_from_file_expecting((char*)&nlights,input,"%d","nlights");
  read_formatted_from_file_expecting((char*)&tSEP,input,"%d","tSEP");
  
  fclose(input);
}

void read_pars(const char *input)
{
  FILE *fin=open_file(input,"r");
  char data_list_file[1024];
  read_formatted_from_file_expecting(data_list_file,fin,"%s","data_list_file");
  fclose(fin);
  
  read_ensemble_pars(base_path,T,ibeta,nmass,mass,nlights,data_list_file);
  
  TH=T/2;
}

string rectangle(int a,int b,jack c)
{
  ostringstream out;
  
  out<<a<<" "<<c.med()-c.err()<<endl<<b<<" "<<c.med()-c.err()<<endl;
  out<<b<<" "<<c.med()+c.err()<<endl<<a<<" "<<c.med()+c.err()<<endl;
  out<<a<<" "<<c.med()-c.err()<<endl;	    
  
  return out.str();
}

int main()
{
  read_pars("input");
  
  jack mass_estim[2];
  jack Z_estim[2];
  int tmin_390[2]={9,7};                                                                                                     
  //loop over heavier mass
  int ncombo=nlights*nmass;
  jvec Z(ncombo,njack),aM(ncombo,njack);
  //WARNING!!!! here order reverted!!!!
  for(int im1=0;im1<nlights;im1++)
    for(int im2=0;im2<nmass;im2++)
      {
	int ic=im1*nmass+im2;
	
	jvec eff[2];
	for(int sm_lev=0;sm_lev<2;sm_lev++)
	  {
	    int tmin=(int)(tmin_390[sm_lev]*afm[ibeta]/afm[1]+0.5);
	    //load the current combo
	    c[sm_lev]=load_corr(base_path,sm_lev,im1,im2);
	    c[sm_lev].print_to_file(combine("out_%02d_%02d_%d.xmg",im1,im2,sm_lev).c_str());
	    
	    //compute effective mass
	    eff[sm_lev]=effective_mass(c[sm_lev].simmetrized(1));
	    
	    //compute eff mass subtracting t+1 with t. Find the point
	    //where the fit with 0 is less than 1
	    {
	      jvec temp(TH-1,njack);
	      for(int t=0;t<TH-1;t++)
		{
		  temp[t]=sqr(eff[sm_lev][t+1]-eff[sm_lev][t]);
		  double ch=temp[t].med()/temp[t].err();
		  //if(ch<1 && t>4 && tmin]==-1) tmin[sm_lev]=t;
		}
	      ofstream out(combine("eff_diff_%02d_%02d_%d.xmg",im1,im2,sm_lev).c_str());
	      out<<temp;
	    }
	    
	    //cout<<im1<<" "<<im2<<" "<<sm_lev<<" "<<tmin[sm_lev]+2<<endl;
	    
	    //fit mass and print it
	    {
	      mass_estim[sm_lev]=constant_fit(eff[sm_lev],tmin,TH);
	      ofstream out(combine("eff_%02d_%02d_%d.xmg",im1,im2,sm_lev).c_str());
	      out<<"@type xydy"<<endl<<"@s0 line type 0"<<endl;
	      out<<eff[sm_lev];
	      out<<"&\n@type xy\n";
	      out<<rectangle(tmin,TH,mass_estim[sm_lev])<<endl;
	    }
	    
	    //compute the corr func removed of the temporal dependency,
	    //fit Z and print it
	    {
	      ofstream out(combine("rem_tdep_%02d_%02d_%d.xmg",im1,im2,sm_lev).c_str());
	      out<<"@type xydy"<<endl<<"@s0 line type 0"<<endl;
	      jvec temp=c[sm_lev].simmetrized(1);
	      for(int iel=0;iel<=TH;iel++)
		for(int ijack=0;ijack<njack+1;ijack++)
		  temp[iel].data[ijack]=temp[iel].data[ijack]/fun_fit(1,1,mass_estim[sm_lev][ijack],iel);
	      Z_estim[sm_lev]=constant_fit(temp,tmin,TH);
	      out<<temp<<"&\n@type xy\n";
	      out<<rectangle(tmin,TH,Z_estim[sm_lev])<<endl;
	    }
	  }
	
	//compute Zl and M
	aM[ic]=(mass_estim[0]+mass_estim[1])/2;
	Z[ic]=Z_estim[0]/sqrt(Z_estim[1]);
	
	//print M and f
	double a=afm[ibeta]/0.197;
	jack M=aM[ic]/a;
	jack af=Z[ic]*(mass[im1]+mass[im2])/aM[ic]/sinh(aM[ic]);
	jack f=af/a;
	
	cout<<im1<<" "<<im2<<" "<<M<<" "<<f<<endl;
      }
  
  aM.write_to_binfile("results");
  (Z*Z).append_to_binfile("results");
  
  return 0;
}
