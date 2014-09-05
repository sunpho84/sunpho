#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "/Users/francesco/Prace/sunpho/Analysis_programs/src/auto.cpp"

autocorr_data_t data(1);

//exit with error message
void crash(const char *temp,...)
{
  char buffer[1024];
  va_list args;

  va_start(args,temp);
  vsprintf(buffer,temp,args);
  va_end(args);

  cerr<<"ERROR: "<<buffer<<endl;
  exit(1);
}

//read allocating
void read(const char *path)
{
  //scan
  ifstream fin(path);
  double t;
  vector<double> buf;
  while(fin>>t) buf.push_back(t);
  data.get_from(buf);
  fin.close();
}

int main(int narg,char **arg)
{
  if(narg<2) crash("Use: %s filein [clust_size]",arg[0]);
  
  //load
  read(arg[1]);
  if(data.size==0) crash("empty file %s",arg[1]);
  else cout<<"data size: "<<data.size<<endl;
  
  if(narg<3) data.clust_size=1;
  else data.clust_size=atoi(arg[2]);
  if(data.clust_size<0) crash("suggested negative clust_size %u",data.clust_size);
  if(data.clust_size>=data.size) crash("suggested too large clust_size %u>buf_size %u",
					 data.clust_size,data.size);
  
  //compute average and error
  double ave,err;
  data.ave_err(ave,err);
  
  //write errors
  double med_tint,err_tint;
  data.compute_tint(med_tint,err_tint,"/tmp/autocorr.xmg");
  cout<<"med_tint: "<<med_tint<<" "<<err_tint<<endl;
  double tau=(2*med_tint-1)/2;
  cout.precision(8);
  cout<<"value: "<<ave<<" +- "<<err*sqrt(2*tau+1)<<endl;
  
  return 0;
}
