#include "include.h"

const int deb=0;

const int njacks=16;
const int T=48;
const int L=24;
double glbl_coeff;
int parity;
int nlevls;
int nml=2,iml=0,imc=1;
int nlevls_sto;
jvec *data;
char corr_file_path[100];

template<class to> to two(to Z1,to Z2,to M,int t)
{return Z1*Z2*exp(-M*L)*cosh(M*(L-t))/M;}

jvec load(const char *path,int itheta,int ism_so_lv,int ism_si_lv)
{
  int ri=0;
  
  jvec a(T,njacks);
  a.load(path,ri+2*(iml+nml*(imc+nml*(ism_si_lv+nlevls_sto*ism_so_lv))));
  
  return glbl_coeff*a;
}

void load_raw_data(int itheta,int *map,const char *raw_data_path=NULL)
{
  //open output if path passed
  ofstream raw_data,eff_mass_raw_data;
  if(raw_data_path!=NULL)
    {
      raw_data.open(raw_data_path);
      raw_data<<"@type xydy"<<endl;
      
      eff_mass_raw_data.open(combine("eff_mass_%s",raw_data_path).c_str());
      eff_mass_raw_data<<"@type xydy"<<endl;
    }
  
  //load data and write the file with all the corrs used
  for(int ism_so=0;ism_so<nlevls;ism_so++)
    for(int ism_si=0;ism_si<=ism_so;ism_si++)
      {
	jvec temp=load(corr_file_path,itheta,map[ism_so],map[ism_si]);
	data[ism_so*nlevls+ism_si]=temp.simmetrized(parity);
	
	//write raw data
	if(raw_data_path!=NULL)
	  {
	    raw_data<<temp<<"&"<<endl;
	    if(parity==1) eff_mass_raw_data<<effective_mass(data[ism_so*nlevls+ism_si])<<"&"<<endl;
	    else          eff_mass_raw_data<<aperiodic_effective_mass(data[ism_so*nlevls+ism_si])<<"&"<<endl;
	  }
      }
  
  //change sink smearead with more precise smeared-source corrs
  for(int ism_so=0;ism_so<nlevls;ism_so++)
    for(int ism_si=ism_so+1;ism_si<nlevls;ism_si++)
      data[ism_so*nlevls+ism_si]=data[ism_si*nlevls+ism_so];
  
  //close output
  if(raw_data_path!=NULL)
    {
      raw_data.close();
      eff_mass_raw_data.close();
    }
}

int main(int narg,char **arg)
{
  if(narg<2) crash("use %s input" ,arg[0]);
  FILE *fin=open_file(arg[1],"r");
  
  //file path
  read_formatted_from_file_expecting(corr_file_path,fin,"%s","corr_file_path");
  
  //read global coeff and parity
  read_formatted_from_file_expecting((char*)&glbl_coeff,fin,"%lg","glbl_coeff");
  read_formatted_from_file_expecting((char*)&parity,fin,"%d","parity");
  
  //allocate nlevls-depending stuff
  read_formatted_from_file_expecting((char*)&nlevls_sto,fin,"%d","nlevls_sto");
  read_formatted_from_file_expecting((char*)&nlevls,fin,"%d","nlevls");
  int *map=(int*)malloc(nlevls*sizeof(int));
  data=(jvec*)malloc(sizeof(jvec)*nlevls*nlevls);
  
  //read levels to be used
  for(int ilev=0;ilev<nlevls;ilev++)
    read_formatted_from_file((char*)&(map[ilev]),fin,"%d","ilevl");
    
  //timeslice for normalization
  int tinv;
  read_formatted_from_file_expecting((char*)&tinv,fin,"%d","tinv");
  
  //timeslice for fit
  int tfit;
  read_formatted_from_file_expecting((char*)&tfit,fin,"%d","tfit");
  
  //interval to fit first state
  int tfit_first_min,tfit_first_max;
  read_formatted_from_file_expecting((char*)&tfit_first_min,fin,"%d","tfit_first");
  read_formatted_from_file((char*)&tfit_first_max,fin,"%d","tfit_first");
  
  fclose(fin);
  
  ////////////////////////////// finish reading input ///////////////////////////
  
  //load data
  load_raw_data(0,map,"raw_data.xmg");
  
  //find eigenvectors and eigenstates
  double diag_m[nlevls*nlevls];
  find_diagonalizing_matrix(diag_m,tinv,tfit,data,nlevls);

  //diagonalize using diag matrix
  jvec corr_d[nlevls];
  for(int ilev=0;ilev<nlevls;ilev++) corr_d[ilev]=jvec(L+1,njacks);
  for(int t=0;t<=L;t++)
    {
      //left mult with transposed
      jvec temp(nlevls*nlevls,njacks);
      for(int ilev=0;ilev<nlevls;ilev++)
	for(int jlev=0;jlev<nlevls;jlev++)
	  {
	    temp[ilev*nlevls+jlev]=0;
	    for(int klev=0;klev<nlevls;klev++)
	      temp[ilev*nlevls+jlev]+=diag_m[klev*nlevls+ilev]*data[klev*nlevls+jlev][t];
	  }
      
      //right mult
      for(int ilev=0;ilev<nlevls;ilev++)
	{
	  corr_d[ilev][t]=0;
	  for(int klev=0;klev<nlevls;klev++)
	    corr_d[ilev][t]+=temp[ilev*nlevls+klev]*diag_m[klev*nlevls+ilev];
	}
    }
  
  //write the coefficients to be used
  cout<<endl;
  cout<<"Taking coefficients at time: "<<tfit<<endl;
  for(int ilev=0;ilev<nlevls;ilev++)
    {
      for(int ieig=0;ieig<nlevls;ieig++)
	cout<<diag_m[ilev*nlevls+ieig]<<"\t";
      cout<<endl;
    }
  
  //write ground state
  ofstream ground_GEVP("ground_GEVP.xmg");
  ground_GEVP<<"@type xydy"<<endl;
  //ground_GEVP<<effective_mass(eig[0])<<"&"<<endl;
  ground_GEVP<<effective_mass(corr_d[0])<<"&"<<endl;
  
  //fit first excited mass
  cout<<endl<<"Ground state mass: ";
  for(int ilev=0;ilev<nlevls;ilev++) cout<<map[ilev]<<" ";
  cout<<"\t"<<constant_fit(effective_mass(corr_d[0]),11,20)<<endl;
  
  //fit first excited mass
  cout<<endl<<"First excited mass: ";
  for(int ilev=0;ilev<nlevls;ilev++) cout<<map[ilev]<<" ";
  cout<<"\t"<<constant_fit(effective_mass(corr_d[1]),tfit_first_min,tfit_first_max,"first_GEVP.xmg")<<endl;
  
  //second excited mass
  if(nlevls>2)
    {
      cout<<endl<<"Second excited mass: ";
      for(int ilev=0;ilev<nlevls;ilev++) cout<<map[ilev]<<" ";
      cout<<"\t"<<constant_fit(effective_mass(corr_d[2]),4,6,"second_GEVP.xmg")<<endl;
    }

  return 0;
}
