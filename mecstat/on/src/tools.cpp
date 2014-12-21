#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <cstdarg>
#include <cstdio>
#include <fstream>
#include <iostream>

#include "geometry.hpp"
#include "parameters.hpp"
#include "random.hpp"

using namespace std;

//crash promptin error message
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

//check if a file exists
bool file_exists(const char *path)
{
  ifstream file(path);
  bool ex=file.good();
  file.close();
  return ex;  
}

//write the configuration to disk
void write_conf(const char *path,int isweep)
{
  ofstream conf_file(path);
  if(!conf_file.good()) crash("opening %s to write conf",path);
  
  if(!(conf_file.write((char*)&isweep,sizeof(int)))) crash("writing isweep");
  //write the rnd status
  for(int site=0;site<=V;site++) conf_file<<gen[site]<<endl;
  if(!conf_file.good()) crash("writing random number generator");
  
  conf_file.close();
}

//read the configuration to disk
void read_conf(int &isweep,const char *path)
{
  ifstream conf_file(path);
  if(!conf_file.good()) crash("opening %s to read conf",path);
  
  if(!(conf_file.read((char*)&isweep,sizeof(int)))) crash("reading isweep");
  //read the rnd status
  for(int site=0;site<=V;site++) conf_file>>gen[site];
  if(!conf_file.good()) crash("reading random number generator");
  
  conf_file.close();
}
