#include "../../src/include.h"

const int nmoms=19;
const int tsep=24;
const int T=128;
const int nconfs=77;
const char templ[]="3pt_ms_0.02661_ms_0.02661_ml_0.000678_maxmom_2_Dt24_G%d";
const int iS=0,iVX=1,iVY=2,iVZ=4,iV0=8,iTX=9,iTY=10,iTZ=12;
const int RE=0,IM=1;
const int moms[nmoms][3]={{-1,-1,+0},{-1,+0,-1},{-1,+0,+0},{-1,+0,+1},{-1,+1,+0},
			  {+0,-1,-1},{+0,-1,+0},{+0,-1,+1},{+0,+0,-1},{+0,+0,+0},
			  {+0,+0,+1},{+0,+1,-1},{+0,+1,+0},{+0,+1,+1},{+1,-1,+0},
			  {+1,+0,-1},{+1,+0,+0},{+1,+0,+1},{+1,+1,+0}};
const int nind_moms=3;
const int ind_moms[nind_moms][3]={{0,0,0},{0,0,1},{0,1,1}};
int ind_of_mom[nmoms];
int deg_of_ind_mom[nind_moms];
void fill_ind_of_mom()
{
  for(int ind=0;ind<nind_moms;ind++) deg_of_ind_mom[ind]=0;
  for(int imom=0;imom<nmoms;imom++)
    {
      int n=0;
      for(int i=0;i<3;i++) n+=abs(moms[imom][i]);
      ind_of_mom[imom]=n;
      deg_of_ind_mom[n]++;
    }
}

void load(jvec *out,const char *in,int ig,int ri)
{
  //open allocate and read
  FILE *fin=open_file(combine(in,ig).c_str(),"r");
  double *data=new double[2*T*nconfs*nmoms];
  if(fread(data,sizeof(double),2*T*nconfs*nmoms,fin)!=2*T*nconfs*nmoms) crash("error reading");
  for(int imom=0;imom<nmoms;imom++) out[imom]=jvec(T,nconfs);
  
  //put in place
  for(int imom=0;imom<nmoms;imom++)
    for(int t=0;t<T;t++)
      for(int iconf=0;iconf<nconfs;iconf++)
	out[imom][t][iconf]=data[ri+2*(iconf+nconfs*(t+T*imom))];
  
  //clusterize
  for(int imom=0;imom<nmoms;imom++)
    out[imom].clusterize();
  
  //free and close
  delete[] data;
  fclose(fin);
}

void print(const char *path,jvec *data)
{
  ofstream out(path);
  out<<"@type xydy"<<endl;
  for(int ind=0;ind<nind_moms;ind++)
    {
      out<<data[ind]<<endl;
      out<<"&"<<endl;
    }
}

int main()
{
  fill_ind_of_mom();
  for(int ind=0;ind<nind_moms;ind++) cout<<ind<<" "<<deg_of_ind_mom[ind]<<endl;
  
  jvec tS[nmoms];
  jvec tVX[nmoms];
  jvec tVY[nmoms];
  jvec tVZ[nmoms];
  jvec tV0[nmoms];
  jvec tTX[nmoms];
  jvec tTY[nmoms];
  jvec tTZ[nmoms];
  load(tS,templ,iS,RE);
  load(tVX,templ,iVX,IM);
  load(tVY,templ,iVY,IM);
  load(tVZ,templ,iVZ,IM);
  load(tV0,templ,iV0,RE);
  load(tTX,templ,iTX,IM);
  load(tTY,templ,iTY,IM);
  load(tTZ,templ,iTZ,IM);
  
  //average
  jvec S[nind_moms];
  jvec V0[nind_moms];
  jvec VK[nind_moms];
  jvec TK[nind_moms];
  for(int ind=0;ind<nind_moms;ind++) S[ind]=V0[ind]=VK[ind]=TK[ind]=jvec(T,nconfs);
  for(int imom=0;imom<nmoms;imom++)
    {
      int ind=ind_of_mom[imom];
      S[ind]+=-tS[imom];
      V0[ind]+=tV0[imom];
      VK[ind]+=tVX[imom]*moms[imom][2];
      VK[ind]+=tVY[imom]*moms[imom][1];
      VK[ind]+=tVZ[imom]*moms[imom][0];
      TK[ind]+=tTX[imom]*moms[imom][2];
      TK[ind]+=tTY[imom]*moms[imom][1];
      TK[ind]+=tTZ[imom]*moms[imom][0];
    }
  
  //truncate on first part
  
  //normalize
  for(int ind=0;ind<nind_moms;ind++)
    {
      S[ind]/=deg_of_ind_mom[ind];
      V0[ind]/=deg_of_ind_mom[ind];
      if(ind)
	{
	  VK[ind]/=deg_of_ind_mom[ind]*ind;
	  TK[ind]/=deg_of_ind_mom[ind]*ind;
	}
    }
  
  //truncate
  for(int ind=0;ind<nind_moms;ind++)
    {
      S[ind]=S[ind].subset(0,tsep+1);
      V0[ind]=V0[ind].subset(0,tsep+1);
      VK[ind]=VK[ind].subset(0,tsep+1);
      TK[ind]=TK[ind].subset(0,tsep+1);
    }

  //print
  print("/tmp/S.xmg",S);
  print("/tmp/V0.xmg",V0);
  print("/tmp/VK.xmg",VK);
  print("/tmp/TK.xmg",TK);

  //write
  for(int ind=0;ind<nind_moms;ind++)
    {
      S[ind].append_to_binfile("/tmp/W");
      V0[ind].append_to_binfile("/tmp/W");
      VK[ind].append_to_binfile("/tmp/W");
      TK[ind].append_to_binfile("/tmp/W");
    }
  
  return 0;
}
