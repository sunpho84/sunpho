#include "on.hpp"

int isweep;
double *phi,*pi;

//allocate all fields
void allocate()
{
  phi=new double[V*N];
  pi=new double[V*N];
}

//deallocate fields
void delocate()
{
  delete[] phi;
  delete[] pi;
}

//init according to start condition
void setup_conf()
{
  if(!file_exists("conf"))
    {
      isweep=0;
      switch(start_cond)
	{
	case HOT:
	  //init_system_to_hot();
	  break;
	case COLD:
	  //init_system_to_cold();
	  break;
	case LOAD: crash("conf not exists!");break;
	}
    }
  else read_conf(isweep,"conf"); //remember rnd gen reinit
}

int main(int narg,char **arg)
{
  //read input and initialize 
  read_input("input");
  init();
  
  return 0;
}
