#include <fstream>
#include <iostream>
#include <string.h>
#include <vector>

#include "debug.hpp"

using namespace std;

//read a column in a vector
vector<double> read_column_in_file(const char *path,int ic)
{
  vector<double> out;
     
  //open file
  ifstream fin(path);
  if(!fin.good()) CRASH("troubles opening file %s",path);
     
  //scan line by line
  char line[256];
  int nl=0;
  while(fin.getline(line,256))
    {
      //check that the line is not empty
      int acc=0;
      for(int len=strlen(line),i=0;i<len;i++)
	{
	  acc|=line[i]!=' ';
	  acc|=line[i]!='\t';
	}
         
      //if the line is acceptable
      if(acc)
	{
	  int nc=0;
	  int found=0;
	  for(char *brk,sep[]=" \t",*str=strtok_r(line,sep,&brk);str!=NULL;str=strtok_r(NULL,sep,&brk))
	    {
	      if(nc==ic)
		{
		  //if found push it back and crash if failed
		  double t;
		  found=sscanf(str,"%lg",&t);
		  if(found) out.push_back(t);
		  else CRASH("entry %d (\"%s\") on line %d of file %s non readable",nc,str,nl,path);
		}
                 
	      nc++;
	    }
             
	  //check to have found something
	  if(!found) CRASH("not found entry %d on line %d of file %s, \"%s\"",ic,nl,path,line);
	}
           
      nl++;
    };
  fin.close();
     
  return out;
}
