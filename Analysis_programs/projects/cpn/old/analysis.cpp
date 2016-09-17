#include "include.h"

int njacks=50;
int term=10000;
int L,tmin;

const double chrono_topo_coeff=0.06;
const double chrono_topo_width=0.125;
const double chrono_topo_barr=3.5;
const double chrono_topo_force_out=0;
//const double chrono_topo_well_tempering=0;
const double chrono_topo_bend=0.0;
const bool symm=true;

double qnongeo_min,qnongeo_max;
vector<double> chrono_topo_past_values;

bool reco_top_flag;
double reco_topo_dq,reco_topo_min,reco_topo_max;
int reco_topo_nq;
double *reco_topo_q,*reco_topo_V;

//read the magnetization
void read_mag_ingredients(jack &num,jack &den)
{
  ifstream in("mag");
  vector<double> n_buf;
  vector<double> d_buf;

  int iconf;
  double tn,td;
  while(in>>iconf>>tn>>td)
    if(iconf>term)
      {
	n_buf.push_back(tn);
	d_buf.push_back(td);
      }
  
  int size=n_buf.size();
  if(size==0) crash("null mag buffer");
  if(size%njacks) crash("njacks %d does not divide mag size %d",njacks,size);
  //if(njacks==-1) njacks=sqrt(size);
  //else if(njacks!=(int)sqrt(size)) crash("njacks for mag does not match: %d vs. %d!",njacks,(int)sqrt(size));
  cout<<"Data size: "<<size<<", NJacks: "<<njacks<<endl;
  jack n(njacks),d(njacks);
  int clust_size=size/njacks;
  for(int ijack=0;ijack<njacks;ijack++)
    {
      n[ijack]=0;
      d[ijack]=0;
      for(int iel=ijack*clust_size;iel<(ijack+1)*clust_size;iel++)
	{
	  n[ijack]+=n_buf[iel];
	  d[ijack]+=d_buf[iel];
	}
      n[ijack]/=clust_size;
      d[ijack]/=clust_size;
    }
  
  n.clusterize();
  d.clusterize();

  num=n;
  den=d;
}  

jack compute_mag(jack &n,jack &d)
{
  return sqrt(n/d-1)/(2*sin(M_PI/L));
}

//read the energy
jack read_energy()
{
  ifstream in("energy");
  vector<double> energy_buf;

  int iconf;
  double tenergy;
  while(in>>iconf>>tenergy)
    if(iconf>term) energy_buf.push_back(tenergy);
  
  int size=energy_buf.size();
  if(size==0) crash("null energy buffer");
  if(size%njacks) crash("njacks %d does not divide energy size %d",njacks,size);
  //if(njacks==-1) njacks=sqrt(size);
  //else if(njacks!=(int)sqrt(size)) crash("njacks for energy does not match: %d vs. %d!",njacks,(int)sqrt(size));
  
  jack energy(njacks);
  int clust_size=size/njacks;
  for(int ijack=0;ijack<njacks;ijack++)
    {
      energy[ijack]=0;
      for(int iel=ijack*clust_size;iel<(ijack+1)*clust_size;iel++)
	energy[ijack]+=energy_buf[iel];
      energy[ijack]/=clust_size;
    }
  
  energy.clusterize();
  
  return energy;
}

//read the corr
jvec read_corr(const char *path)
{
  ifstream in(path);
  vector<double> corr_buf;

  int iconf;
  double d,c;
  while(in>>iconf>>d>>c)
    if(iconf>term)
      corr_buf.push_back(c);
  if(corr_buf.size()%(L/2+1)) crash("correlation file size %u not multiple of L/2+1=%d",corr_buf.size(),L/2+1);

  int size=corr_buf.size()/(L/2+1);
  
  if(size==0) crash("null corr buffer");
  if(size%njacks) crash("njacks %d does not divide corr size %d",njacks,size);
  
  jvec corr(L/2+1,njacks);
  corr=0;
  
  int clust_size=size/njacks;
  for(int d=0;d<L/2+1;d++)
    for(int ijack=0;ijack<njacks;ijack++)
      {
	corr[d][ijack]=0;
	for(int iel=ijack*clust_size;iel<(ijack+1)*clust_size;iel++)
	  corr[d][ijack]+=corr_buf[iel*(L/2+1)+d];
      }
  corr/=clust_size;
  corr.clusterize();
  
  return corr;
}

//read the topology
void read_topology(double &med_tint,double &err_tint,jack &topo,jack &topo_susc)
{
  ifstream in("topology");
  map<int,vector<double>> topo_buf;
  
  int iconf,istep,maxstep=0;
  double geo,nongeo;
  qnongeo_min=+1e300;
  qnongeo_max=-1e300;
  map<int,double> sum,sum2;
  while(in>>iconf>>istep>>geo>>nongeo)
    {
      if(istep==2) chrono_topo_past_values.push_back(nongeo);
      
      if(iconf>term)
	{
	  if(istep>maxstep) maxstep=istep;
	  if(istep==0) topo_buf[-1].push_back(geo);
	  topo_buf[istep].push_back(nongeo);
	  sum[istep]+=(nongeo-geo);
	  sum2[istep]+=(nongeo-geo)*(nongeo-geo);
	  
	  qnongeo_min=min(qnongeo_min,nongeo);
	  qnongeo_max=max(qnongeo_max,nongeo);
	}
    }
  
  for(istep=-1;istep<=maxstep;istep++)
    {
      string path=combine("plots/topo_%d.xmg",istep);
      ofstream topo_file(path.c_str());
      if(!topo_file.good()) crash("opening \"%s\"",path.c_str());
      for(auto it=topo_buf[istep].begin();it!=topo_buf[istep].end();it++)
	topo_file<<*it<<endl;
      
      sum[istep]/=topo_buf[istep].size();
      sum2[istep]/=topo_buf[istep].size();
      cout<<"StdVar istep "<<istep<<": "<<sqrt(sum2[istep]-sum[istep]*sum[istep])<<endl;
    }
  
  int size=topo_buf[-1].size();
  if(size==0) crash("null topo buffer");
  if(size%njacks) crash("njacks %d does not divide topo size %d",njacks,size);
  //if(njacks==-1) njacks=sqrt(size);
  //else if(njacks!=(int)sqrt(size)) crash("njacks for topo does not match: %d vs. %d!",njacks,(int)sqrt(size));
  
  topo=jack(njacks);
  topo_susc=jack(njacks);
  int clust_size=size/njacks;
  for(int ijack=0;ijack<njacks;ijack++)
    {
      topo[ijack]=0;
      topo_susc[ijack]=0;
      for(int iel=ijack*clust_size;iel<(ijack+1)*clust_size;iel++)
	{
	  topo[ijack]+=topo_buf[-1][iel];
	  topo_susc[ijack]+=topo_buf[-1][iel]*topo_buf[-1][iel];
	}
    }
  topo/=clust_size;
  topo_susc/=clust_size;
  topo.clusterize();
  topo_susc.clusterize();
  
  //compute_tint(med_tint,err_tint,topo_buf[-1],sqrt(topo_buf[-1].size()));
}

//compute the topodynamical potential using past history
double compute_theta_pot(double Q,bool ave=false)
{
  if(chrono_topo_past_values.size()==0) crash("not loaded topology");
  
  //compute 
  double harm=0;
  if(Q<-chrono_topo_barr) harm=chrono_topo_force_out*sqr(-Q-chrono_topo_barr)/2;
  if(Q>+chrono_topo_barr) harm=chrono_topo_force_out*sqr(+Q-chrono_topo_barr)/2;
  
  //put inside the barrier
  if(Q<-chrono_topo_barr) Q=-chrono_topo_barr;
  if(Q>+chrono_topo_barr) Q=+chrono_topo_barr;
  
  double topotential=0;
  int size=chrono_topo_past_values.size();
  int weight=(ave?(size-term-1):1);
  harm*=-weight;
  int i=0;
  for(auto it=chrono_topo_past_values.begin();it!=chrono_topo_past_values.end();it++)
    {
      double q=(*it);
      double diff=Q-q,f=diff/chrono_topo_width;
      double cont=exp(-f*f/2+chrono_topo_bend*Q*Q/2);
      topotential+=cont*weight;
      if(ave && i>term) weight--;
      i++;
    }
  
  //check
  if(ave && weight!=0) 
    {
      cerr<<"Error: "<<weight<<endl;
      exit(1);
    }

  topotential*=chrono_topo_coeff;
    
  return (topotential+harm)/(ave?(size-term-1):1);
}

double interp_potential(double xint)
{
  if(!reco_top_flag) crash("reconstruct first the potential");
  
  int iq=(xint-reco_topo_min+reco_topo_dq)/reco_topo_dq;
  if(iq<1||iq>(reco_topo_nq-1)) crash("iq=%d, nq=%d",iq,reco_topo_nq);
  
  double x[3]={reco_topo_q[iq-1],reco_topo_q[iq],reco_topo_q[iq+1]};
  double y[3]={reco_topo_V[iq-1],reco_topo_V[iq],reco_topo_V[iq+1]};
  
  return parabolic_spline(x,y,xint);
}

//read the correlation and separate in different topological sectors
void read_correlation_per_topology()
{
  ifstream intop("topology");
  
  //store the distances
  double dist[31];

  //read topology
  map<int,int> topo_map;
  map<int,double> topo_map_nongeo;
  int max_g=0;
  while(intop.good())
    {
      int tconf,istep;
      double geo,nongeo;
      intop>>tconf>>istep>>geo>>nongeo;
      
      if(tconf>term)
	{
	  int g=(int)(fabs(geo)+0.5);
	  topo_map[tconf]=g;
	  if(istep==2) topo_map_nongeo[tconf]=nongeo;
	  max_g=max(max_g,g);
	}
    }

  //allocate
  jvec *corr_per_topo=new jvec[max_g+1];
  int nconf_per_topo[max_g+1];
  int iconf_per_topo[max_g+1];
  int clust_size_per_topo[max_g+1];
  for(int g=0;g<=max_g;g++)
    {
      corr_per_topo[g]=jvec(L/2+1,njacks);
      corr_per_topo[g]=0;
      nconf_per_topo[g]=0;
      iconf_per_topo[g]=0;
    }
  
  //count topo per sector
  ifstream tmpcorrd("corrd");
  int ntot_conf=0;
  while(tmpcorrd.good())
    {
      //read corrd
      int cconf;
      for(int t=0;t<L/2+1;t++)
	{
	  double d,c;
	  tmpcorrd>>cconf>>d>>c;
	}
      
      if(tmpcorrd.good() && cconf>term)
	{
	  nconf_per_topo[topo_map[cconf]]++;
	  ntot_conf++;
	}
    }
  tmpcorrd.close();
  
  //compute cluster size
  for(int g=0;g<=max_g;g++)
    {
      clust_size_per_topo[g]=nconf_per_topo[g]/njacks;
      nconf_per_topo[g]=clust_size_per_topo[g]*njacks;
      cout<<"sector: "<<g<<" "<<nconf_per_topo[g]<<endl;
    }
  
  //load correlation functions
  ifstream incorrd("corrd");
  double pot0=compute_theta_pot(0,true);
  double reweighted_corr[L/2+1],reweighted_corr_err[L/2+1];
  double tot_weight=0;
  for(int t=0;t<L/2+1;t++) reweighted_corr[t]=reweighted_corr_err[t]=0;

  int jtot_conf=0;
  int nconsidered=0;
  map<int,double> histoQ;
  while(incorrd.good())
    {
      //read corrd
      int gup=-1; //g to up
      int cconf;
      double weight=0;
      for(int t=0;t<L/2+1;t++)
	{
	  double d,c;
	  incorrd>>cconf>>d>>c;
	  
	  if(cconf>term && fabs(topo_map_nongeo[cconf])<chrono_topo_barr*0.8)
	    {
	      //take topological sector and store dist
	      int g=topo_map[cconf];
	      dist[t]=d;
	      
	      //reweight
	      if(t==0)
		{
		  double pot=interp_potential(topo_map_nongeo[cconf]);
		  weight=exp(pot-pot0);
		  histoQ[topo_map[cconf]]+=weight;
		  tot_weight+=weight;
		  nconsidered++;
		}
	      reweighted_corr[t]+=c*weight;
	      reweighted_corr_err[t]+=c*c*weight;
	      
	      //put in the appropriate topological sector
	      if(iconf_per_topo[g]<nconf_per_topo[g] && clust_size_per_topo[g])
		{
		  int iclust=iconf_per_topo[g]/clust_size_per_topo[g];
		  corr_per_topo[g][t][iclust]+=c;
		  gup=g;
		}
	    }
	}
      
      //increase the conf
      if(gup!=-1) iconf_per_topo[gup]++;
      if(cconf>term) jtot_conf++;
    }
  
  cout<<"Reweighted topocharge histogram"<<endl;
  for(auto it=histoQ.begin();it!=histoQ.end();it++) cout<<it->first<<" "<<it->second/tot_weight<<endl;
  
  cout<<"Assuming geo=nongeo"<<endl;
  map<int,double> histo_assume;
  double norm=0;
  for(auto it=histoQ.begin();it!=histoQ.end();it++) norm+=(histo_assume[it->first]=exp(interp_potential(it->first)-interp_potential(0)));
  for(auto it=histo_assume.begin();it!=histo_assume.end();it++) cout<<it->first<<" "<<it->second/norm<<endl;
  
  //write the clust
  ofstream plot_per_topo("plots/corrd_per_topo.xmg");
  plot_per_topo<<"@type xydy"<<endl;
  jvec totcorr(L/2+1,njacks);
  int ntot=0;
  totcorr=0;
  for(int g=0;g<=max_g;g++)
    if(nconf_per_topo[g])
      {
	corr_per_topo[g].clusterize(clust_size_per_topo[g]);
	//cout<<corr_per_topo[g][0]<<endl;
	//plot_per_topo<<effective_mass(corr_per_topo[g])*sqrt(2)<<endl;
	for(int t=0;t<L/2+1;t++) plot_per_topo<<corr_per_topo[g][t]<<endl;
	plot_per_topo<<"&"<<endl;
	
	totcorr+=corr_per_topo[g]*nconf_per_topo[g];
	ntot+=nconf_per_topo[g];
      }
  //plot the average
  totcorr/=ntot;
  //plot_per_topo<<effective_mass(totcorr)*sqrt(2)<<endl;
  for(int t=0;t<L/2+1;t++) plot_per_topo<<totcorr[t]<<endl;
  plot_per_topo<<"&"<<endl;
  //plot the reweighted
  for(int t=0;t<L/2+1;t++)
    {
      reweighted_corr[t]/=tot_weight;
      reweighted_corr_err[t]/=tot_weight;
      reweighted_corr_err[t]-=reweighted_corr[t]*reweighted_corr[t];
      reweighted_corr_err[t]=sqrt(reweighted_corr_err[t]/tot_weight);
      
      cout<<t<<" "<<reweighted_corr[t]<<" "<<reweighted_corr_err[t]<<endl;
      plot_per_topo<<t<<" "<<reweighted_corr[t]<<" "<<reweighted_corr_err[t]<<endl;
    }
}

void draw_topotential(bool ave=false)
{
  ofstream out("plots/topotential.xmg");
  
  reco_top_flag=true;
  reco_topo_dq=chrono_topo_width/2;
  reco_topo_min=qnongeo_min*1.1;
  reco_topo_max=qnongeo_max*1.1;
  reco_topo_nq=(reco_topo_max-reco_topo_min)/reco_topo_dq;
  reco_topo_q=new double[reco_topo_nq+1];
  reco_topo_V=new double[reco_topo_nq+1];
  for(int iq=0;iq<=reco_topo_nq;iq++)
    {
      reco_topo_q[iq]=reco_topo_min+iq*reco_topo_dq;
      reco_topo_V[iq]=compute_theta_pot(reco_topo_q[iq],ave);
      if(symm) reco_topo_V[iq]=(reco_topo_V[iq]+compute_theta_pot(-reco_topo_q[iq],ave))/2;
      out<<reco_topo_q[iq]<<" "<<reco_topo_V[iq]<<endl;
    }
  
  out.close();
}

int main(int narg,char **arg)
{
  if(narg<3) crash("Use: %s L term tmin",arg[0]);
  L=atoi(arg[1]);
  term=atoi(arg[2]);
  tmin=atoi(arg[3]);
  
  //read correlation functions
  jvec corrw=read_corr("corr");
  jvec corrd=read_corr("corrd");

  corrw.print_to_file("plots/corrw_raw.xmg");
  corrd.print_to_file("plots/corrd_raw.xmg");
  
  corrd*=sqrt(2);
  
  //fit E and A
  debug_fit=0;
  jack E,Z2;
  two_pts_fit(E,Z2,corrw,tmin,L/2,"plots/corrw.xmg");
  jack A_w=Z2/E,xi_w=1/E;
  two_pts_fit(E,Z2,corrd,tmin,L/2);
  E*=sqrt(2);
  write_constant_fit_plot("plots/corrd.xmg",effective_mass(corrd)*sqrt(2),E,tmin,L/2);
  jack A_d=Z2/E,xi_d=1/E;
  cout<<"Xi_w: "<<smart_print(xi_w)<<", Xi_d: "<<smart_print(xi_d)<<", Xi_d/Xi_w: "<<smart_print(xi_d/xi_w)<<endl;
  cout<<"A_w: "<<smart_print(A_w)<<", A_d: "<<smart_print(A_d)<<", A_d/A_w: "<<smart_print(A_d/A_w)<<endl;
  cout<<"Side in units of Xi_w: "<<smart_print(L/xi_w)<<endl;

  //read magnetization ingredients and compute it
  jack num,den;
  read_mag_ingredients(num,den);
  jack chi_m=num;
  jack xi_g=compute_mag(num,den);
  cout<<"Chi_m: "<<smart_print(chi_m)<<endl;
  cout<<"Xi_g: "<<smart_print(xi_g)<<endl;
  jack A_g=chi_m/sqr(xi_g)*xi_w;
  cout<<"A_g: "<<smart_print(A_g)<<endl;
  
  //read the energy
  jack energy=read_energy();
  cout<<"Energy: "<<smart_print(energy)<<endl;
  
  //read the topological charge
  autocorr_debug=false;
  jack topo,topo_susc;
  double topo_med_tint,topo_err_tint;
  read_topology(topo_med_tint,topo_err_tint,topo,topo_susc);
  cout<<"Topological charge: "<<smart_print(topo)<<endl;
  cout<<"Topological susceptibility: "<<smart_print(topo_susc)<<endl;
  //cout<<"Topological autocorrelation time: "<<smart_print(topo_med_tint,topo_err_tint)<<endl;

  draw_topotential(true);
  
  read_correlation_per_topology();
  
  return 0;
}
