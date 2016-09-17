#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

using namespace std;

struct N_t
{
  int N;
  int N0;
};

int V=1728;
vector<N_t> data;
const double meta_sigma_N=10;
const double meta_sigma_N0=10;
const double meta_coeff=1000;

//increase the meta-potential relative to a certain occupation
double compute_meta_potential(int occ_N,int occ_N0,int itraj)
{
  double P=0;
  double pref=meta_coeff/(meta_sigma_N*meta_sigma_N0*sqrt(2*M_PI));
  for(int jtraj=0;jtraj<itraj;jtraj++)
    {
      double dN=(occ_N-data[jtraj].N)/meta_sigma_N;
      double dN0=(occ_N0-data[jtraj].N0)/meta_sigma_N0;
      
      double d2=dN*dN+dN0*dN0;
      if(d2<9) P+=exp(-d2/2);
    }
  P*=pref;
  
  return P;
}

int main()
{
  ifstream in("/tmp/data.xmg");
  N_t temp;
  int min_N=V,max_N=0;
  int min_N0=V,max_N0=-V;
  while(in>>temp.N0>>temp.N)
    {
      data.push_back(temp);
      min_N=std::min(min_N,temp.N);
      min_N0=std::min(min_N0,temp.N0);
      max_N=std::max(max_N,temp.N);
      max_N0=std::max(max_N0,temp.N0);
    }
  int ntraj=data.size();
  
  double Mpot=0;
  int delta=1;
  for(int N=min_N;N<=max_N;N+=delta)
    {
      for(int N0=std::max(-V+N,min_N0);N0<=std::min(V-N,max_N0);N0+=delta)
	{
	  double pot=compute_meta_potential(N,N0,ntraj);
	  Mpot=std::max(pot,Mpot);
	}
      cout<<N<<" "<<Mpot<<endl;
    }

  cout<<Mpot<<endl;
  int niso=10;
  double del_Mpot=Mpot/niso;

  std::map<int,std::vector<pair<int,int> > > potiso;
  for(int N=min_N;N<=max_N;N+=delta)
    {
      int ipot=0;
      for(int N0=std::max(-V+N,min_N0);N0<=std::min(V-N,max_N0);N0+=delta)
	{
	  double pot=compute_meta_potential(N,N0,ntraj);
	  int jpot=pot/del_Mpot;
	  //if(ipot!=jpot) 
	  //potiso[std::min(ipot,jpot)].push_back(std::make_pair(N,N0));
	  ipot=jpot;
	  potiso[jpot].push_back(std::make_pair(N,N0));
	}
    }
  
  for(std::map<int,std::vector<pair<int,int> > >::iterator it=potiso.begin();it!=potiso.end();it++)
    {
      for(std::vector<pair<int,int> >::iterator jt=it->second.begin();jt!=it->second.end();jt++)
	cerr<<jt->first<<" "<<jt->second<<endl;
      cerr<<"&"<<endl;
    }
  
  return 0;
}
