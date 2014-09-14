#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <cstdarg>
#include <cstdio>
#include <fstream>
#include <iostream>

#include "data.hpp"
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

//write the configuration to disk
void write_conf(const char *path)
{
  ofstream conf_file(path);
  if(!conf_file.good()) crash("opening %s to write conf",path);
  
  //write the two pieces
  if(!(conf_file.write((char*)zeta,sizeof(dcomplex)*N*V))) crash("writing zeta");
  if(!(conf_file.write((char*)lambda,sizeof(dcomplex)*V*NDIMS))) crash("writing lambda");
  //write the rnd status
#ifndef GOOD_GENERATOR
  if(!(conf_file.write((char*)(&gen),sizeof(rnd_gen)))) crash("writing random generator");
#endif  
  
  conf_file.close();
}

//read the configuration to disk
void read_conf(const char *path)
{
  ifstream conf_file(path);
  if(!conf_file.good()) crash("opening %s to read conf",path);
  
  //read the two pieces
  if(!(conf_file.read((char*)zeta,sizeof(dcomplex)*N*V))) crash("reading zeta");
  if(!(conf_file.read((char*)lambda,sizeof(dcomplex)*V*NDIMS))) crash("reading lambda");
#ifndef GOOD_GENERATOR
  if(!(conf_file.read((char*)(&gen),sizeof(rnd_gen)))) crash("reading random generator");
#endif  
  
  conf_file.close();
}
