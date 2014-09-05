#ifndef _TYPES_HPP
#define _TYPES_HPP

#define WHI kBlue
#define RED kRed

//structur to hold a tile
struct tile_t
{
  bool color;
  int norie_poss;
  bool fc[6];
  tile_t(bool color,int norie_poss,bool c0,bool c1,bool c2,bool c3,bool c4,bool c5) :
    color(color),norie_poss(norie_poss)
  {
    fc[0]=c0;
    fc[1]=c1;
    fc[2]=c2;
    fc[3]=c3;
    fc[4]=c4;
    fc[5]=c5;
  };
};

//holds a position
struct pos_t
{
  int shift[2];
  int v[6];
  pos_t(int s0,int s1,int v0,int v1,int v2,int v3,int v4,int v5)
  {
    shift[0]=s0;
    shift[1]=s1;
    v[0]=v0;
    v[1]=v1;
    v[2]=v2;
    v[3]=v3;
    v[4]=v4;
    v[5]=v5;
  };    
};

typedef long long int conf_t;
typedef int top_t[19];

#endif
