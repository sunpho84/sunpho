#include "common.cpp"

//configurations to skip
const int nskip=8;
int skip_conf[nskip]={920,1080,1400,1660,1700,1720,1740,1880};
vector<int> skip_iconf;
bool check_blacklisted(int conf)
{
  int is_blacklist=0;
  for(int iskip=0;iskip<nskip;iskip++)
    if(skip_conf[iskip]==conf)
      is_blacklist=1;
  
  return is_blacklist;
}

//parity and real/imag
const int RE=0,IM=1;
const int EVN=1,ODD=-1;

//parameters to load a correlator with insertion
struct corr_t
{
  int ig;
  int ri;
  int par;
  corr_t(int ig,int ri,int par) : ig(ig),ri(ri),par(par) {}
};
corr_t raw_S0(0,RE,EVN);
corr_t raw_VX(1,IM,ODD);
corr_t raw_VY(2,IM,ODD);
corr_t raw_VZ(4,IM,ODD);
corr_t raw_V0(8,RE,EVN);
corr_t raw_TX(9,IM,ODD);
corr_t raw_TY(10,IM,ODD);
corr_t raw_TZ(12,IM,ODD);

//list the momentum
int nmom=19;
vector<array<int,4>> mom(nmom);

//check or store conf id
vector<int> conf_list(njacks,-1);
void check_store_conf(int iconf,int conf_id)
{
  if(conf_list[iconf]==-1) conf_list[iconf]=conf_id;
  //else if(conf_list[iconf]!=conf_id) crash("for conf %d/%d expected %d obtained %d",iconf,njacks,conf_list[iconf],conf_id);
}

//load meson
jvec load_meson(semi_t semi,int imes)
{
  jvec out(T,njacks);
  string path=semi.get_name()+"_mes"+to_string(imes)+"_"+semi.get_mes_name(imes)+"_pp";
  cout<<"Opening "<<path<<endl;
  FILE *fin=open_file(path.c_str(),"r");
  
  int iline=0;
  int ijack=0,iconf_read=0;
  
  do
    {
      //check if was asked to skip
      int skip=0;
      for(auto &s : skip_iconf) if(iconf_read==s) skip=1;
      
      for(int t=0;t<T;t++)
	{
	  int t_read;
	  double reim[2];
	  int rc=fscanf(fin,"%d %lg %lg",&t_read,&reim[0],&reim[1]);
	  if(rc!=3) crash("error, expecting 3 obtained %d at line %d",rc,iline);
	  iline++;
	  
	  if(t_read!=t) crash("obtained t %d when expecting %d",t_read,t);
	  if(!skip) out[t][ijack]=reim[RE];
	}
      
      //increase iconf_read and if not skipped, also ijack
      if(!skip) ijack++;
      // else cout<<"Skipping conf_read "<<iconf_read<<endl;
      iconf_read++;
    }
  while(ijack<njacks);
  
  out.clusterize();
  
  return out.simmetrized(1);
}
//load a whole file
vector<jvec> load(string path,int tsep,int ri)
{
  vector<jvec> data(nmom,jvec(T,njacks));
  skip_iconf.clear();
  
  //open and skip line
  cerr<<"Opening file "<<path<<endl;
  FILE *fin=popen(("grep -v \\# "+path).c_str(),"r");
  
  int iline=0;
  int iconf=0,iconf_read=0,last_conf_read=-1;
  
  do
    {
      //int exp_conf=min_conf+conf_offset*iconf;
      for(int imom=0;imom<nmom;imom++)
	for(int t=0;t<T;t++)
	  {
	    int conf_read,zmom,ymom,xmom,t_read;
	    double reim[2];
	    int rc=fscanf(fin,"%d %d %d %d %d %lg %lg",&conf_read,&zmom,&ymom,&xmom,&t_read,&reim[0],&reim[1]);
	    if(feof(fin)) crash("reached EOF prematurely after %d confs",iconf);
	    if(rc!=7) crash("error, expecting 7 obtained %d at line %d",rc,iline);
	    iline++;
	    
	    //increas iconf if different from previous one and not in the blacklist
	    int is_blacklisted=check_blacklisted(conf_read);
	    if(last_conf_read!=conf_read)
	      {
		if(last_conf_read!=-1) iconf_read++;
		
		if(is_blacklisted) skip_iconf.push_back(iconf_read);
		else
		  if(last_conf_read!=-1)
		    iconf++;
	      }
	    last_conf_read=conf_read;
	    
	    // if(imom==0 && t==0)  cout<<"ijack "<<iconf<<", iconf_read "<<iconf_read<<", conf name"<<conf_read<<endl;
	    
	    //cout<<exp_conf<<" "<<imom<<" "<<t<<" "<<conf_read<<" "<<xmom<<" "<<ymom<<" "<<zmom<<" "<<t_read<<" "<<re<<" "<<im<<endl;
	    check_store_conf(iconf,conf_read);
	    if(t_read!=t) crash("obtained t %d when expecting %d",t_read,t);
	    
	    if(!is_blacklisted) data[imom][t][iconf]=reim[ri];
	    mom[imom][0]=xmom;
	    mom[imom][1]=ymom;
	    mom[imom][2]=zmom;
	    
	    int ind_mom=sqr(xmom)+sqr(ymom)+sqr(zmom);
	    mom[imom][3]=ind_mom;
	  }
    }
  while(iconf<njacks);
  
  //clusterize, print, invert and takes only first tsep
  for(auto &d : data)
    {
      d.clusterize();
      
      //d.print_to_file("plots/clusters/"+path+".xmg");
      d=d.simmetric().subset(0,tsep+1);
    }
  
  // cout<<"List of index of skipped conf for file "<<path<<endl;
  // for(auto s : skip_iconf) cout<<s<<endl;
  
  return data;
}

//load using info stored into corr_t
vector<jvec> load(semi_t semi,int tsep,corr_t c)
{
  string base_path=semi.get_name()+"_deltat_"+to_string(tsep)+"_MaxComp_2_gammaop_";
  return load(base_path+to_string(c.ig),tsep,c.ri);
}

//load a time correlator
vector<jvec> load_0(semi_t semi,int tsep,corr_t S)
{
  vector<jvec> temp=load(semi,tsep,S);
  
  //average space and norm
  vector<jvec> c0(nind_mom,jvec(tsep+1,njacks));
  vector<int> norm(nind_mom,0);
  
  //sum all non-ind mom
  for(int imom=0;imom<nmom;imom++)
    {
      int ind_mom=mom[imom][3];
      norm[ind_mom]++;
      c0[ind_mom]+=temp[imom];
    }
  
  //normalize
  for(int ind_mom=0;ind_mom<nind_mom;ind_mom++)
    c0[ind_mom]/=norm[ind_mom];
  
  return c0;
}
vector<jvec> load_S0(semi_t semi,int tsep) {return load_0(semi,tsep,raw_S0);}
vector<jvec> load_V0(semi_t semi,int tsep) {return load_0(semi,tsep,raw_V0);}

//load a spatial correlator
vector<jvec> load_K(semi_t semi,int tsep,corr_t X,corr_t Y,corr_t Z)
{
  vector<jvec> c[3];
  c[0]=load(semi,tsep,X);
  c[1]=load(semi,tsep,Y);
  c[2]=load(semi,tsep,Z);
  
  //average space and norm
  vector<jvec> cK(nind_mom,jvec(tsep+1,njacks));
  vector<int> norm(nind_mom,0);
  
  for(int imom=0;imom<nmom;imom++)
    for(int mu=0;mu<3;mu++)
      if(mom[imom][mu])
	  {
	    int ind_mom=mom[imom][3];
	    norm[ind_mom]++;
	    cK[ind_mom]+=c[mu][imom]*mom[imom][mu];
	  }
  
  //normalize
  for(int ind_mom=0;ind_mom<nind_mom;ind_mom++)
    if(norm[ind_mom])
      cK[ind_mom]/=norm[ind_mom];
  
  return cK;
}
vector<jvec> load_VK(semi_t semi,int tsep) {return load_K(semi,tsep,raw_TX,raw_TY,raw_TZ);}
vector<jvec> load_TK(semi_t semi,int tsep) {return load_K(semi,tsep,raw_VX,raw_VY,raw_VZ);}
vector<jvec> (*load_semi[4])(semi_t semi,int tsep)={load_S0,load_VK,load_V0,load_TK};

void convert_all_corrs()
{
  for(int isemi=0;isemi<nsemi;isemi++)
    {
      //load 3pts
      for(int itsep=0;itsep<ntseps;itsep++)
      	for(int icur=0;icur<ncurrs;icur++)
      	  {
      	    vector<jvec> s=load_semi[icur](semi[isemi],tsep_list[itsep]);
      	    for(int imom=0;imom<nind_mom;imom++) corr_3pts[icorr_3pts(isemi,itsep,icur,imom)]=s[imom];
      	  }
      //load 2pts
      for(int i12=1;i12<=2;i12++)
	{
	  semi_t &s=semi[isemi];
	  corr_2pts[icorr_2pts(isemi,i12)]=load_meson(s,i12);
	}
    }
   write_all_corrs();
}
