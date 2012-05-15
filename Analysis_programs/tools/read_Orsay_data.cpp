#include <stdio.h>
#include <stdlib.h>

typedef double complex[2]; 

typedef complex color[3];
typedef color spincolor[4];
typedef spincolor colorspincolor[3];
typedef colorspincolor spincolorspincolor[4];

//fucking tool to revert the endianess of doubles
void doubles_to_doubles_changing_endianess(double *dest,double *sour,int ndoubles)
{
    char *cdest,*csour;
    char temp;
  
    if(dest==sour)
	for(int idouble=0;idouble<ndoubles;idouble++)
	{
	    cdest=(char*)(dest+idouble);
	    csour=(char*)(sour+idouble);
        
	    temp=csour[7];
	    csour[7]=cdest[0];
	    cdest[0]=temp;

	    temp=csour[6];
	    csour[6]=cdest[1];
	    cdest[1]=temp;

	    temp=csour[5];
	    csour[5]=cdest[2];
	    cdest[2]=temp;

	    temp=csour[4];
	    csour[4]=cdest[3];
	    cdest[3]=temp;
	}
    else
	for(int idouble=0;idouble<ndoubles;idouble++)
	{
	    cdest=(char*)(dest+idouble);
	    csour=(char*)(sour+idouble);
        
	    cdest[0]=csour[7];
	    cdest[1]=csour[6];
	    cdest[2]=csour[5];
	    cdest[3]=csour[4];
	    cdest[4]=csour[3];
	    cdest[5]=csour[2];
	    cdest[6]=csour[1];
	    cdest[7]=csour[0];

	}
}

int imom_true(int t,int x,int y,int z,int L2P1)
{return x+L2P1*(y+L2P1*(z+L2P1*t));}

int imom(int t,int x,int y,int z,int L,int T,int L2P1)
{return imom_true(t+T,x+L,y+L,z+L,L2P1);}

int main()
{
    int LTOT=24;
    int TTOT=48;
    
    int L=LTOT/4;
    int T=TTOT/4;
    
    int L2P1=2*L+1;
    int T2P1=2*T+1;
    
    int V=L2P1*L2P1*L2P1*T2P1;
    
    spincolorspincolor *A=malloc(V*sizeof(spincolorspincolor));
    if(A==NULL)
    {
	fprintf(stderr,"Not allocated\n");
	exit(1);
    }
    
    ////////////////////////////////////////////////////////////////////////////////
    
    char path[1024]="propagator.mass7.2770.fft";
    
    FILE *fin=fopen(path,"r");
    if(fin==NULL)
    {
	fprintf(stderr,"Could not open file %s\n",path);
	exit(1);
    }
    
    int n=fread(A,sizeof(spincolorspincolor),V,fin);
    if(n!=V)
    {
	fprintf(stderr,"Error, read %d instead than %d momenta\n",n,V);
	exit(1);
    }
    
    doubles_to_doubles_changing_endianess((double*)A,(double*)A,V);
    
    int tmin[2]={0,4};
    int tmax[2]={3,9};
    int smin[2]={0,2};
    int smax[2]={2,5};

    FILE *fout=fopen("temp","w");
    if(fout==NULL)
    {
	fprintf(stderr,"Error opening output!\n");
	exit(1);
    }
    
    int imo=0;

    for(int inte=0;inte<2;inte++)
	for(int t=tmin[inte];t<=tmax[inte];t++)
	    for(int y=smin[inte];y<=smax[inte];y++)
		for(int z=smin[inte];z<=smax[inte];z++)
		    for(int x=smin[inte];x<=smax[inte];x++)
		    {		
			for(int ic_so=0;ic_so<3;ic_so++)
			    for(int id_so=0;id_so<4;id_so++)
				for(int ic_si=0;ic_si<3;ic_si++)
				    for(int id_si=0;id_si<4;id_si++)
					for(int ri=0;ri<2;ri++)
					{
					    int i=imom(t,x,y,z,L,T,L2P1);
					    
					    int nw=fwrite(&(A[i][id_so][ic_so][id_si][ic_si][ri]),sizeof(double),1,fout);
					    if(nw!=1)
					    {
						fprintf(stderr,"Error in writing\n");
						exit(1);
					    }
					}
			imo++;
		    }
    printf("imo=%d\n",imo);
    
    fclose(fout);
    
    fclose(fin);
    
    ////////////////////////////////////////////////////////////////////////////////
    
    free(A);
    
    return 0;
}
