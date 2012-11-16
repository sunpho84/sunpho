#include <include.h>

int T;
int njack=16;

int *icombo;
char **data_path;
int tmin,tmax;

int ndata;

jvec load_2pts(int idata)
{
  return jvec_load(data_path[idata],T,njack,icombo[idata]);
}

void read_input()
{
  FILE *input_file=open_file("input","r");

  read_formatted_from_file_expecting((char*)(&ndata),input_file,"%d","ndata");
  data_path=(char**)malloc(ndata*sizeof(char*));
  icombo=(int*)malloc(ndata*sizeof(int));
  for(int idata=0;idata<ndata;idata++)
    {
      data_path[idata]=(char*)malloc(1024);
      read_formatted_from_file_expecting(data_path[idata],input_file,"%s","data_path");
      read_formatted_from_file((char*)&(icombo[idata]),input_file,"%d","icombo");
    }
  
  read_formatted_from_file_expecting((char*)(&T),input_file,"%d","T");
  
  read_formatted_from_file_expecting((char*)(&tmin),input_file,"%d","tmin");
  read_formatted_from_file_expecting((char*)(&tmax),input_file,"%d","tmax");
  
  fclose(input_file);
}

int main()
{
  read_input();
  
  ///////////////////////////// Load two points for standing D and D* //////////////////////////
  
  jvec CSL(T,njack);
  jvec CSS(T,njack);
  
  jack M,ZL,ZS;
  int idata_SL=0;
  for(int idata=0;idata<ndata;idata++)
    {
      CSL=(CSL*idata+load_2pts(idata))/(idata+1);
      data_path[idata][strlen(data_path[idata])-2]='3';
      if(file_exists(data_path[idata]))
	{
	  CSS=(CSS*idata_SL+load_2pts(idata))/(idata_SL+1);
	  idata_SL++;
	}
      
      //compute mass
      two_pts_SL_fit(M,ZL,ZS,CSL.simmetrized(1),CSS.simmetrized(1),tmin,tmax,tmin,tmax,combine("M_%02d.xmg",idata).c_str());
      cout<<"nsources: "<<idata<<" mass: "<<M<<endl;
      
      if(idata==ndata-1)
	{
	  M.write_to_binfile("M");
	  ZL.write_to_binfile("ZL");
	}
    }
  
  CSL.write_to_binfile("CSL");
  
  return 0;
}
