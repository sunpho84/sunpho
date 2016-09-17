#include <iostream>

#include "cpn.hpp"

int main()
{
  read_pars_t read_pars;
  read_input(read_pars,"input");
  
  int jsweep=0;
  if(use_topo_pot==2)
    {
      ifstream topology_file("topology");
      chrono_topo.init();
      
      do
	{
	  for(int jlev=0;jlev<=nstout_lev;jlev++)
	    {
	      double topo_num,topo_sim;
	      int isweep,ilev;
	      if(topology_file>>isweep>>ilev>>topo_sim>>topo_num)
		{
		  if(ilev!=jlev) crash("reading ilev");
		  if(ilev==nstout_lev) chrono_topo.update(isweep,topo_num);
		  if(isweep!=jsweep) crash("obtained wrong sweep, %d when expecting %d",isweep,jsweep);
		}
	    }
	  
	  jsweep++;
	}
      while(topology_file.good() && !topology_file.eof());
      
      chrono_topo.save();
    }
  
  return 0;
}
