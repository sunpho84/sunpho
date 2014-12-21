#ifndef _TOOLS_HPP
#define _TOOLS_HPP

#include <fstream>
#include <string>

using namespace std;

void crash(const char *temp,...);
bool file_exists(const char *path);
void read_conf(int &itraj,const char *path);
void write_conf(const char *path,int itraj);

//read an element from input file
template <class T> void read(T &out,ifstream &in,string is)
{
  string s;
  if(!(in>>s)) crash("impossible to read expect string \"%s\"",is.c_str());
  if(s!=is) crash("obtained %s while reading %s",s.c_str(),is.c_str());
  if(!(in>>out)) crash("reading data");
}


#endif
