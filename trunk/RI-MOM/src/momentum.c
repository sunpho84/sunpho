#pragma once

void define_canonical_momentum_set()
{
  int iP[2][4][2]={{{0,3},{0,2},{0,2},{0,2}},{{4,7},{2,3},{2,3},{2,3}}};
  int lP[4],size[4]={T,L,L,L},imom=0;
  
  for(int iset=0;iset<2;iset++)
    {
      for(lP[0]=iP[iset][0][0];lP[0]<=iP[iset][0][1];lP[0]++)
	for(lP[1]=iP[iset][1][0];lP[1]<=iP[iset][1][1];lP[1]++)
	  for(lP[2]=iP[iset][2][0];lP[2]<=iP[iset][2][1];lP[2]++)
	    for(lP[3]=iP[iset][3][0];lP[3]<=iP[iset][3][1];lP[3]++)
	      {
		P2[imom]=SinP2[imom]=SinP4[imom]=0;
		for(int idir=0;idir<4;idir++)
		  {
		    if(idir==0) size[idir]=T;
		    else size[idir]=L;
		    
		    P[idir][imom]=2*PI*(lP[idir]+bc[idir]*0.5)/size[idir];
		    P2[imom]+=pow(P[idir][imom],2);
		    SinP[idir][imom]=sin(P[idir][imom]);
		    SinP2[imom]+=pow(SinP[idir][imom],2);
		    SinP4[imom]+=pow(SinP[idir][imom],4);
		  }
		imom++;
	      }
    }
}
