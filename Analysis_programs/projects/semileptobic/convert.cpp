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
  double norm;
  corr_t(int ig,int ri,int par,double norm) : ig(ig),ri(ri),par(par),norm(norm) {}
};
corr_t raw_S0(0,RE,EVN,1);
corr_t raw_VX(1,IM,ODD,1);
corr_t raw_VY(2,IM,ODD,1);
corr_t raw_VZ(4,IM,ODD,1);
corr_t raw_V0(8,RE,EVN,-1);
corr_t raw_TX(9,IM,ODD,1);
corr_t raw_TY(10,IM,ODD,1);
corr_t raw_TZ(12,IM,ODD,1);

//list the momentum
int nmom=19;
vector<array<int,4>> mom(nmom);

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
  int iconf_used=0,iconf_read=0,last_conf_label_read=-1;
  
  do
    {
      //int exp_conf=min_conf+conf_offset*iconf;
      for(int imom=0;imom<nmom;imom++)
	for(int t=0;t<T;t++)
	  {
	    int conf_label,zmom,ymom,xmom,t_read;
	    double reim[2];
	    int rc=fscanf(fin,"%d %d %d %d %d %lg %lg",&conf_label,&zmom,&ymom,&xmom,&t_read,&reim[0],&reim[1]);
	    if(feof(fin)) crash("reached EOF prematurely after %d confs read",iconf_read);
	    if(rc!=7) crash("error, expecting 7 obtained %d at line %d",rc,iline);
	    iline++;
	    
	    //increas iconf if different from previous one and not in the blacklist
	    int is_blacklisted=check_blacklisted(conf_label);
	    if(last_conf_label_read!=conf_label)
	      {
		if(last_conf_label_read!=-1) iconf_read++;
		
		if(is_blacklisted) skip_iconf.push_back(iconf_read);
		else
		  if(last_conf_label_read!=-1)
		    iconf_used++;
	      }
	    last_conf_label_read=conf_label;
	    
	    // if(imom==0 && t==0)  cout<<"ijack "<<iconf<<", iconf_read "<<iconf_read<<", conf name"<<conf_read<<endl;
	    
	    //cout<<exp_conf<<" "<<imom<<" "<<t<<" "<<conf_read<<" "<<xmom<<" "<<ymom<<" "<<zmom<<" "<<t_read<<" "<<re<<" "<<im<<endl;
	    if(t_read!=t) crash("obtained t %d when expecting %d",t_read,t);
	    
	    if(!is_blacklisted) data[imom][t][iconf_used]=reim[ri]*T;
	    mom[imom][0]=xmom;
	    mom[imom][1]=ymom;
	    mom[imom][2]=zmom;
	    
	    int ind_mom=sqr(xmom)+sqr(ymom)+sqr(zmom);
	    mom[imom][3]=ind_mom;
	  }
      iconf_read++;
    }
  while(iconf_used<njacks);
  
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
      c0[ind_mom]+=temp[imom]*S.norm;
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
  corr_t M[3]={X,Y,Z};
  vector<jvec> c[3];
  for(int i=0;i<3;i++) c[i]=load(semi,tsep,M[i]);
  
  //average space and norm
  vector<jvec> cK(nind_mom,jvec(tsep+1,njacks));
  vector<int> norm(nind_mom,0);
  
  for(int imom=0;imom<nmom;imom++)
    for(int mu=0;mu<3;mu++)
      if(mom[imom][mu])
	{
	  int ind_mom=mom[imom][3];
	  norm[ind_mom]++;
	  cK[ind_mom]+=c[mu][imom]*mom[imom][mu]*M[mu].norm;
	}
  
  //normalize
  for(int ind_mom=0;ind_mom<nind_mom;ind_mom++)
    if(norm[ind_mom])
      cK[ind_mom]/=norm[ind_mom];
  
  return cK;
}
vector<jvec> load_TK(semi_t semi,int tsep) {return load_K(semi,tsep,raw_TX,raw_TY,raw_TZ);}
vector<jvec> load_VK(semi_t semi,int tsep) {return load_K(semi,tsep,raw_VX,raw_VY,raw_VZ);}
vector<jvec> (*load_semi[ncurrs])(semi_t semi,int tsep);

void convert_all_corrs()
{
  load_semi[S0]=load_S0;
  load_semi[V0]=load_V0;
  load_semi[VK]=load_VK;
  load_semi[TK]=load_TK;
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

vector<vector<jvec>> read_unaveraged(string path,int ireim)
{
  int ntso=T/2;
  //tso  combo
  vector<vector<jvec> > data(ntso);
  for(auto &d : data) d.resize(nmom,jvec(T,njacks));
  
  string binarized_path="unaveraged/"+path+"_binarized";
  if(!file_exists(binarized_path))
    {
      FILE *fin=popen(("grep -v \\# unaveraged/"+path).c_str(),"r");
      for(int ijack=0;ijack<njacks;ijack++)
	for(int itso=0;itso<ntso;itso++)
	  for(int imom=0;imom<nmom;imom++)
	    for(int t=0;t<T;t++)
	      {
		int iconf,tso,mel,tre,bo;
		double reim[2];
		int rc=fscanf(fin,"%d %d %d %d %d %d %d %d %d %lg %lg",&iconf,&tso,&tre,&mel,&bo,&bo,&bo,&bo,&bo,&reim[0],&reim[1]);
		if(rc!=11) crash("%d",rc);
		if(feof(fin)) crash("premature end");
		
		data[itso][imom][t][ijack]=reim[ireim];
	      }
      
      FILE *fout=open_file(binarized_path,"w");
      for(auto &d : data)
	for(auto &da : d)
	  {
	    da.clusterize();
	    da.write_to_binfile(fout);
	  }
      fclose(fout);
    }
  else
    {
      FILE *fin=open_file(binarized_path,"r");
      for(auto &d : data) for(auto &da : d) da.load(fin);
      fclose(fin);
    }
  
  // ofstream out("/tmp/temp.xmg");
  // out<<"@type xydy"<<endl;
  // jvec ave(T,njacks);
  // for(int itso=0;itso<ntso;itso++)
  //   {
  //     ave+=data[itso][0];
  //     out<<data[itso][0]<<endl;
  //   }
  // out<<ave/ntso<<endl;
  
  return data;
}

vector<jvec> read_unaveraged_meson(string path,int ireim)
{
  int ntso=T/2;
  //tso  combo
  vector<jvec> data(ntso,jvec(T,njacks));
  
  string binarized_path="unaveraged/"+path+"_binarized";
  if(!file_exists(binarized_path))
    {
      FILE *fin=popen(("grep -v \\# unaveraged/"+path).c_str(),"r");
      for(int ijack=0;ijack<njacks;ijack++)
	for(int itso=0;itso<ntso;itso++)
	  for(int t=0;t<T;t++)
	    {
	      int tre,bo;
	      double reim[2];
	      int rc=fscanf(fin,"%d %d %d %d %d %d %lg %lg",&tre,&bo,&bo,&bo,&bo,&bo,&reim[0],&reim[1]);
	      if(rc!=8) crash("%d",rc);
	      if(feof(fin)) crash("premature end");
	      
	      data[itso][t][ijack]=reim[ireim];
	    }
      
      FILE *fout=open_file(binarized_path,"w");
      for(auto &d : data)
	{
	  d.clusterize();
	  d.write_to_binfile(fout);
	}
      fclose(fout);
    }
  else
    {
      FILE *fin=open_file(binarized_path,"r");
      for(auto &d : data) d.load(fin);
      fclose(fin);
    }
  
  // ofstream out("/tmp/temp.xmg");
  // out<<"@type xydy"<<endl;
  // jvec ave(T,njacks);
  // for(int itso=0;itso<ntso;itso++)
  //   {
  //     ave+=data[itso];
  //     out<<data[itso]<<endl;
  //   }
  // out<<ave/ntso<<endl;
  
  return data;
}
