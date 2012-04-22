#include <include.h>
#include <iostream>

#include "prace_common.cpp"

using namespace std;

int main()
{
  read_data_list();
  int ith2=ith_spec;
  int im2=0;
  int ru=1;
  
  jvec c(T,njack);
  c=read_two_points(ith2,im2,ru,im_spec,ru,0);
  cout<<c<<endl;
  
  return 0;
}
