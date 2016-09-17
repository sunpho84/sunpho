#include "common.cpp"

int main(int narg,char **arg)
{
  if(narg<2) crash("Use %s file1",arg[0]);
  
  int njacks=16;
  int nm=9;
  
  jvec A(nm,njacks);
  
  A.load(arg[1],0);
  
  for(int i=0;i<8;i++)
    {
      jack M=A[i+1]/A[i];
      double triv_err=sqrt(sqr(A[i+1].err()/A[i].med())+sqr(A[i+1].med()/sqr(A[i].med())*A[i].err()));
      cout<<" "<<M<<" "<<triv_err/M.err()<<endl;
    }
  
  return 0;
}
