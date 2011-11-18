#pragma once

void check_endianess()
{
  big_endian=1;
  big_endian=(int)(*(char*)(&big_endian));
}

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

//open file
FILE* open_file(const char *path,const char *mode)
{
  FILE *out=fopen(path,mode);
  if(out==NULL)
    {
      fprintf(stderr,"Error trying to open file: %s\n",path);
      exit(1);
    }
  else printf("Opened file: %s for %s\n",path,mode);
  
  return out;
}

//read a bunch of doubles
void read_doubles(double *out,FILE *fin,int n)
{
  double *buf=big_endian ? malloc(sizeof(double)*n) : out;
  int nr=fread(buf,sizeof(double),n,fin);
  fread(out,1,1,fin);
  if(!feof(fin))
    {
      fprintf(stderr,"Error, not reached end of file!\n");
      exit(1);
    }
  if(nr!=n)
    {
      fprintf(stderr,"Error reading from file!\n");
      exit(1);
    }
  if(big_endian)
    {
      doubles_to_doubles_changing_endianess(out,buf,nr);
      free(buf);
    }

  for(int i=0;i<n;i++) printf("%lg\n",out[i]);
}

//read a whole ape propagators
void read_ape_propagators(ccss_propagator **out,const char *path,int nmass)
{
  FILE *fin=open_file(path,"r");
  
  ape_propagator buf[2][nmass];
  read_doubles((double*)buf,fin,2*nmass*sizeof(ape_propagator)/sizeof(double));
  for(int imass=0;imass<nmass;imass++)
    for(int r=0;r<2;r++)
      for(int imom=0;imom<nmom;imom++)
	for(int ic_so=0;ic_so<3;ic_so++)
	  for(int id_so=0;id_so<4;id_so++)
	    for(int ic_si=0;ic_si<3;ic_si++)
	      for(int id_si=0;id_si<4;id_si++)
		{
		  memcpy(out[r][imass][imom][ic_si][ic_so][id_si][id_so],buf[r][imass][imom][ic_so][id_so][ic_si][id_si],sizeof(complex));
		  printf("%d %d %d %d %d %d %d %lg %lg %lg %lg\n",imass,r,imom,ic_so,id_so,ic_si,id_si,out[r][imass][imom][ic_si][ic_so][id_si][id_so][0],out[r][imass][imom][ic_si][ic_so][id_si][id_so][1],buf[r][imass][imom][ic_so][id_so][ic_si][id_si][0],buf[r][imass][imom][ic_so][id_so][ic_si][id_si][1]);
		}
}

//read a list of propagator file clusterizing them
void read_ccss_propagator_set(ccss_propagator ***out,const char *path_template,int init_conf_num,int nconf,int skip,int njack,int nmass)
{
  int clust_size=nconf/njack;
  
  for(int ijack=0;ijack<=njack;ijack++)
    for(int r=0;r<2;r++)
      for(int imass=0;imass<nmass;imass++)
	memset(out[ijack][r][imass],0,sizeof(ccss_propagator));
  
  for(int iconf=0;iconf<nconf;iconf++)
    {
      char path[1024];
      sprintf(path,path_template,init_conf_num+iconf*skip);
      read_ape_propagators(out[njack],path,nmass); //use njack el as buf
      int iclus=iconf/clust_size;
      printf("%d %d\n",iclus,clust_size);
      for(int r=0;r<2;r++)
	for(int imass=0;imass<nmass;imass++)
	  ccss_propagator_summassign(out[iclus][r][imass],out[njack][r][imass]);
    }
  
  for(int r=0;r<2;r++)
    for(int imass=0;imass<nmass;imass++)
      memset(out[njack][r][imass],0,sizeof(ccss_propagator));
  
  for(int r=0;r<2;r++)
    for(int imass=0;imass<nmass;imass++)
      {
	/*
	for(int iclus=0;iclus<njack;iclus++)
	  ccss_propagator_summassign(out[njack][r][imass],out[iclus][r][imass]);

	for(int iclus=0;iclus<njack;iclus++)
	  ccss_propagator_subt(out[iclus][r][imass],out[njack][r][imass],out[iclus][r][imass]);
	for(int ijack=0;ijack<njack;ijack++)
	  ccss_propagator_prodassign_real(out[ijack][r][imass],1.0/(nconf-clust_size));
	ccss_propagator_prodassign_real(out[njack][r][imass],1.0/nconf);
	*/
      }
}
