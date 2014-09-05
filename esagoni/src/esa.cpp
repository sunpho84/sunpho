#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "plot.hpp"
#include "pos.hpp"
#include "tile.hpp"
#include "types.hpp"

#include <cmath>
#include <iostream>
#include <unistd.h>

const int debug=0;
const bool debug_new=false;
const int draw_all=0;

using namespace std;

const int write_each=1000;

#define CHECK_ORDER
#define CHECK_TAB_ORIE

//find a tile jumping already chosen ones
int find_tile_from_other_tile(int *is_chosen,int other_tile)
{
  int tile=0;

  while(is_chosen[tile]||(other_tile!=0))
    {
      if(!is_chosen[tile]) other_tile--;
      tile++;
    }
  is_chosen[tile]=1;
  
  return tile;
}

//find the top (tile occupying position)
void unmap_lop_from_conf(top_t lop,conf_t conf)
{
  //find the entry among the leftovers
  for(int ipos=18;ipos>=0;ipos--)
    {
      //mark it in the position and diminish configuration
      conf_t temp=conf/(19-ipos);
      lop[ipos]=conf-temp*(19-ipos);
      conf=temp;
    }
}

//write down a top
void write_top(top_t top)
{
  for(int ipos=0;ipos<19;ipos++) cout<<top[ipos]<<" ";
  cout<<endl;
}

//adjust for alredy chosen
void unmap_top_from_lop(top_t top,top_t lop)
{
  int is_chosen[19]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  for(int ipos=0;ipos<19;ipos++)
    {
      int tile=find_tile_from_other_tile(is_chosen,lop[ipos]);
      top[ipos]=tile;
      is_chosen[tile]=1;
    }
  
  if(debug) write_top(top);
}

//check if a top conf is valid
int check_valid_top(top_t top,int *orie)
{
  int has_been_3=0,has_been_5=0,has_been_9=0;
  
  int vertex[54]={0};
  
  //find if the current configuration works
  bool work=true;
  //loop over all positions
  int ipos=0;
  do
    {
      int itile=top[ipos];
      
#ifdef CHECK_ORDER
      //mark 3,5 and 9
      switch(itile)
	{
	case 3: has_been_3=1;break;
	case 5: has_been_5=1;break;
	case 9: has_been_9=1;break;
	case 10: if(!has_been_5) {work=0;if(debug_new) cout<<"Failing to satisfy 5-10"<<endl;}break;
	case 11: if(!has_been_3) {work=0;if(debug_new) cout<<"Failing to satisfy 3-11"<<endl;}break;
	case 16: if(!has_been_9) {work=0;if(debug_new) cout<<"Failing to satisfy 9-16"<<endl;}break;
	}      
#endif
      
#ifdef CHECK_TAB_ORIE
      if(work && (itile==pivot_tile && !(ipos==sext_00||ipos==sext_01||ipos==sext_11||(ipos==sext_22 && orie[ipos]==0))))
	{
	  work=0;
	  if(debug_new) cout<<"Failing to put in a correct place tile 1"<<endl;
	}
#endif
      
      //loop on corners
      int icor=0;
      if(work)
	do
	  {
	    //find the vertex associated to the corner of the the position
	    int ver=grid[ipos].v[icor];
	    
	    //find old and proposed occupation
	    int old=vertex[ver];
	    int pro=tiles[itile].fc[(6-orie[ipos]+icor)%6]+1; //add 1: 0=empty, 1=white, 2=red
	    
	    //accept it if corner was empty
	    if(old!=0) work=(old==pro);
	    
	    //if still working mark it down
	    if(!work && debug_new) cout<<"Failed to satisfy vertex "<<ver<<" of pos "<<ipos<<endl;
	    if(work) vertex[ver]=pro;	  
	    if(debug) cout<<" ipos: "<<ipos<<", vertex: "<<ver<<" old: "<<old<<" pro: "<<pro<<", work: "<<work<<endl;
	    
	    //if still working increment the position
	    if(work) icor++;
	  }
	while(work && icor<6);
      
      //if still working increment position
      if(work) ipos++;
    }
  while(work && ipos<19);
  
  if(!work) if(debug) cout<<ipos<<endl;
  
  return ipos;
}

//increment the leftovers
void increment_lop(top_t lop,int ipos)
{
  lop[ipos]++;
  for(int jpos=ipos+1;jpos<19;jpos++) lop[jpos]=0;
}

//convert back from lop to conf
conf_t map_lop_to_conf(top_t lop)
{
  conf_t conf=0;
  for(int ipos=0;ipos<19;ipos++) conf=conf*(19-ipos)+lop[ipos];

  return conf;
}

//increment the configuration
bool increment_conf(conf_t &conf,top_t lop,top_t top,int *orie_pos,int nworking)
{
  bool finished_all=false;
  for(int itile=nworking+1;itile<19;itile++) orie_pos[itile]=lop[itile]=0;
  
  //loop until working conf
  bool working;

  do
    {
      working=true;
      //if can, rotate it
      if(orie_pos[nworking]!=(tiles[top[nworking]].norie_poss-1)) orie_pos[nworking]++;
      else
	{
	  //rotate to zero
	  orie_pos[nworking]=0;
	  
	  //if not finished all the tries for the site, take another tile
	  if(lop[nworking]!=18-nworking) increment_lop(lop,nworking);
	  else
	    if(nworking!=0)
	      {
		lop[nworking]=orie_pos[nworking]=0;
		nworking--;
		working=false;
	      }
	    else finished_all=true;
	}
    }
  while(!working);
  
  conf=map_lop_to_conf(lop);

  return finished_all;
}

int main(int narg,char **arg)
{
  //compute the number of plausable configurations
  conf_t nplaus=1;
  for(int i=1;i<=19;i++) nplaus*=i;
  cout<<"nplaus: "<<nplaus<<endl;
  
  //initialize the graphic
#ifdef PLOT
  init_graph();
#endif
  
  //store the configuration
  conf_t conf=0;
  int orie_pos[19]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  cout<<"Start: "<<conf<<" ";
  {
    top_t lop,top;
    unmap_lop_from_conf(lop,conf);
    unmap_top_from_lop(top,lop);
    for(int i=0;i<19;i++) cout<<orie_pos[i]<<","<<top[i]<<"|";
    cout<<" "<<int(10000*(double)conf/nplaus)/100.0<<"% "<<endl;
  }
  
  //number of working
  int nworking=0;
  
  //loop over all the configurations
  long long int nfound=0;
  bool finished_all=false;
  do
    {      
      //loop until find a valid conf
      top_t lop,top;
      do
	{
	  //unmap the configuration
	  unmap_lop_from_conf(lop,conf);
	  unmap_top_from_lop(top,lop);

#ifdef PLOT
	  if(draw_all)
	    {
	      draw_conf(top,orie_pos);
	      sleep(1);
	    }
#endif
	  
	  //increment 
	  nworking=check_valid_top(top,orie_pos);
	  finished_all=increment_conf(conf,lop,top,orie_pos,nworking);
	  
	  //cout<<"Increased: "<<conf<<" |";
	  //for(int i=0;i<19;i++) cout<<orie_pos[i]<<","<<top[i]<<"|";
	  //cout<<" "<<int(10000*(double)conf/nplaus)/100.0<<"% "<<nfound<<endl;
	  
	  //draw it
#ifdef PLOT
	  if(nworking==19 && nfound%write_each==0) draw_conf(top,orie_pos);
#endif

	}
      while(nworking!=19 && !finished_all);
      
      //mark as found
      if(!finished_all)
	{
	  if(nworking==19 && nfound%write_each==0)
	    {
	      cout<<"Found: "<<conf<<" |";
	      for(int i=0;i<19;i++) cout<<orie_pos[i]<<","<<top[i]<<"|";
	      cout<<" "<<int(10000*(double)conf/nplaus)/100.0<<"% "<<nfound<<endl;
	    }
	  
	  //if still plausible increase
	  finished_all=increment_conf(conf,lop,top,orie_pos,18);
	  nfound++;
	}
    }
  while(!finished_all);
  
  cout<<"Finished,nfound: "<<nfound<<endl;
  
#ifdef PLOT
  myapp.Run();
#endif

  return 0;
}
