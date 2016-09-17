#include "lib.hpp"

vector<vector<double>> non_geo_Q;
vector<int> geo_Q;

int max_geo_Q;
int min_geo_Q;

//read topology
void read_topology()
{
  cout<<"Reading topology"<<endl;
  
  //allocate data
  chrono_topo.setup();
  non_geo_Q.resize(nstout_lev+1);
  for(int ilev=0;ilev<=nstout_lev;ilev++) non_geo_Q[ilev].resize(nsweep);
  geo_Q.resize(nsweep);
  
  FILE *fintop=open_file("topology","r");
  
  min_geo_Q=+100000;
  max_geo_Q=-100000;
  for(int isweep=0;isweep<nsweep;isweep++)
    for(int jlev=0;jlev<=nstout_lev;jlev++)
      {
	int iconf,ilev;
	double geo,non_geo;
	//if(fscanf(fintop,"%d %d %lg %lg",&iconf,&ilev,&geo,&non_geo)!=4) crash("reading topo at isweep %d ilev %d",isweep,ilev);
	char siconf[16],silev[3],sgeo[20],snon_geo[20];
	if(fscanf(fintop,"%s %s %s %s",siconf,silev,sgeo,snon_geo)==0) crash("reading topo at isweep %d ilev %",isweep,jlev);
	iconf=atoi(siconf);
	ilev=atoi(silev);
	geo=strtod(sgeo,NULL);
	non_geo=strtod(snon_geo,NULL);
	
	if(isweep!=iconf) crash("isweep=%d!=iconf=%d",isweep,iconf);
	if(ilev!=jlev) crash("jlev=%d!=ilev=%d",jlev,ilev);
	
	//approx to closest integer the geometrical
	geo_Q[isweep]=round(geo);
	
	//store non geo
	non_geo_Q[ilev][isweep]=non_geo;
	if(ilev==nstout_lev) chrono_topo.coll[isweep]=non_geo;
	
	//find min/max
	min_geo_Q=min(min_geo_Q,geo_Q[isweep]);
	max_geo_Q=max(max_geo_Q,geo_Q[isweep]);
      }
  fclose(fintop);
}

//read all the data files
void read_data()
{
  int init_read_time=time(0);
  
  read_energy();
  read_magnetization();
  read_correlation();
  read_topology();  
  
  cout<<">>Reading time: "<<time(0)-init_read_time<<" s"<<endl;
}

//================================= analysis ================================

//compute the autocorrelation time of the topological charge
template <class T> double compute_topo_tint(const char *path,T &Q)
{
  fftw_complex *in=(fftw_complex*)fftw_malloc((nsweep-nterm)*sizeof(fftw_complex));
  fftw_complex *out=(fftw_complex*)fftw_malloc((nsweep-nterm)*sizeof(fftw_complex));
  
  //fill and take fft
  for(int i=0;i<nsweep-nterm;i++)
    {
      in[i][0]=Q[i+nterm];
      in[i][1]=0;
    }
  fftw_plan pf=fftw_plan_dft_1d(nsweep-nterm,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_execute(pf);
  fftw_destroy_plan(pf);
  
  //take the square, remove zero mode and put normalization
  double norm=1.0/sqr(nsweep-nterm);
  for(int i=0;i<nsweep-nterm;i++)
    {
      in[i][0]=(out[i][0]*out[i][0]+out[i][1]*out[i][1])*norm;
      in[i][1]=0;
    }
  in[0][0]=in[0][1]=0;
  
  //compute back
  fftw_plan pb=fftw_plan_dft_1d(nsweep-nterm,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
  fftw_execute(pb);
  fftw_destroy_plan(pb);
  
  //print the corr
  ofstream out_corr(path);
  for(int i=0;i<(nsweep-nterm)/2;i++) out_corr<<out[i][0]/out[0][0]<<endl;
  
  free(in);
  free(out);
  
  return 0;
}

//compute the autocorrelation time of the mag
template <class T> double compute_mag_tint(const char *path,T &M)
{
  int n=nsweep/compute_corr_each;
  n-=n%3;
  fftw_complex *in=(fftw_complex*)fftw_malloc(2*n/3*sizeof(fftw_complex));
  fftw_complex *out=(fftw_complex*)fftw_malloc(2*n/3*sizeof(fftw_complex));
  
  //fill and take fft
  for(int i=n/3;i<n;i++)
    {
      in[i-n/3][0]=M[i];
      in[i-n/3][1]=0;
    }
  fftw_plan pf=fftw_plan_dft_1d(2*n/3,in,out,FFTW_FORWARD,FFTW_ESTIMATE);  
  fftw_execute(pf);
  fftw_destroy_plan(pf);
  
  //take the square, remove zero mode and put normalization
  double norm=1.0/sqr(2*n/3);
  for(int i=0;i<2*n/3;i++)
    {
      in[i][0]=(out[i][0]*out[i][0]+out[i][1]*out[i][1])*norm;
      in[i][1]=0;
    }
  in[0][0]=in[0][1]=0;
    
  //compute back
  fftw_plan pb=fftw_plan_dft_1d(2*n/3,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
  fftw_execute(pb);
  fftw_destroy_plan(pb);
  
  //print the corr
  ofstream out_corr(path);
  for(int i=0;i<2*n/3/2;i++) out_corr<<i*compute_corr_each<<" "<<out[i][0]/out[0][0]<<endl;
  
  free(in);
  free(out);
  
  return 0;
}

//compute bare histogram of Q and reweighted ones
void draw_Q_histograms()
{
  for(int rew=0;rew<2;rew++)
    {
      map<int,boot_t> histo_Q;
      for(int Q=min_geo_Q;Q<=max_geo_Q;Q++)
	histo_Q[Q].populate([Q,rew](int isweep)->double{return (geo_Q[isweep]==Q)*(rew?(chrono_topo.weight[isweep]):1);});
      
      //normalize factor
      boot_t norm=0;
      for(auto &it : histo_Q) norm+=it.second;
      for(auto &it : histo_Q) it.second/=norm;
      
      //write the histo
      ofstream histo_geo_Q(rew?"plots/histo_geo_Q_rew.xmg":"plots/histo_geo_Q.xmg");
      histo_geo_Q<<"@type xydy"<<endl;
      for(auto it : histo_Q) histo_geo_Q<<it.first<<" "<<it.second.ave_err()<<endl;
    }
}

//compute the non-integer part of the noise
void plot_noise(double *Q_reno)
{
  cout<<">>Writing noise"<<endl;
  
  ofstream out_sigma("plots/sigma_noise.xmg");
  
  //allocate noise and plan
  fftw_complex *noise=(fftw_complex*)fftw_malloc((nsweep-nterm)*sizeof(fftw_complex));
  fftw_plan pf=fftw_plan_dft_1d(nsweep-nterm,noise,noise,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_plan pb=fftw_plan_dft_1d(nsweep-nterm,noise,noise,FFTW_BACKWARD,FFTW_ESTIMATE);
  
  //extract the noise from each corr
  FILE *fout_noise=open_file("plots/noise.xmg","w");
  FILE *fout_power=open_file("plots/noise_power_spectrum.xmg","w");
  FILE *fout_corre=open_file("plots/noise_correlation.xmg","w");
  for(int ilev=0;ilev<=nstout_lev;ilev++)
    {
      //compute the noise and print it
      for(int isweep=nterm;isweep<nsweep;isweep++)
	{
	  double n=geo_Q[isweep]-non_geo_Q[ilev][isweep]/Q_reno[ilev];
	  while(n<-0.5) n+=1;
	  while(n>+0.5) n-=1;
	  noise[isweep-nterm][0]=n;
	  noise[isweep-nterm][1]=0;
	  fprintf(fout_noise,"%lg\n",n);
	}
      fprintf(fout_noise,"&\n");
      
      //take fft and print
      fftw_execute(pf);
      for(int isweep=0;isweep<nsweep-nterm;isweep++)
	{
	  double n2=sqr(noise[isweep][0])+sqr(noise[isweep][1]);
	  noise[isweep][0]=n2/(nsweep-nterm);
	  noise[isweep][1]=0;
	  
	  if(isweep<2000) fprintf(fout_power,"%lg\n",n2);
	}
      fprintf(fout_power,"&\n");
            
      //take back fft and print
      fftw_execute(pb);
      for(int isweep=0;isweep<2000;isweep++) fprintf(fout_corre,"%lg\n",noise[isweep][0]/noise[0][0]);
      fprintf(fout_corre,"&\n");
      
      out_sigma<<ilev<<" "<<sqrt(noise[0][0]/(nsweep-nterm))<<endl;
    }
  
  out_sigma.close();
  
  fclose(fout_noise);
  fclose(fout_power);
  fclose(fout_corre);
  fftw_destroy_plan(pf);
  fftw_destroy_plan(pb);
  free(noise);
}

int main()
{
  int init_time=time(0);
  read_input(false);
  
  //read the data and reconstruct reweigthing factor
  read_data();
  cout<<"NSweeps: "<<nsweep<<endl;
  
  //draw random indices
  prepare_boot_sample();
  
  //perform the reweighting
  chrono_topo.mega_reweight();
  
  //compute autocorrelation time
  compute_topo_tint("plots/topo_correlation.xmg",geo_Q);
  compute_topo_tint("plots/nongeo_topo_correlation.xmg",non_geo_Q[nstout_lev]);
  
  //===================== now analysis can start ==================
  
  draw_Q_histograms();
  
  //write the topology
  if(nsweep<10000000)
    {
      cout<<"Writing geo topology"<<endl;
      ofstream out_topo("plots/geo_topology.xmg");
      for(int i=nterm;i<nsweep;i++) out_topo<<geo_Q[i]<<endl;
      ofstream out_topo_non_geo("plots/non_geo_topology.xmg");
      for(int i=nterm;i<nsweep;i++) out_topo_non_geo<<non_geo_Q[nstout_lev][i]<<endl;
    }
  
  //compute average renormalization factors and Q2
  double Q_reno[nstout_lev+1];
  //boot_t geo_Q4([](int iconf)->double{return sqr(sqr(geo_Q[iconf]));});
  //boot_t geo_Q6([](int iconf)->double{return pow(geo_Q[iconf],6);});
  //cout<<"Geo topo chi4: "<<((geo_Q4-3*sqr(geo_Q2))/geo_Q2).ave_err()<<endl;
  //cout<<"Geo topo chi6: "<<((geo_Q6-15*geo_Q2*geo_Q4+30*pow(geo_Q2,3))/geo_Q2).ave_err()<<endl;
  
  //susceptibility
  boot_t geo_Q1([](int iconf)->double{return geo_Q[iconf];});
  boot_t geo_Q2([](int iconf)->double{return sqr(geo_Q[iconf]);});
  boot_t geo_susc=(geo_Q2)/(L*L)*1e4;
  cout<<"Geo topo charge: "<<smart_print(geo_Q1.ave_err())<<endl;
  cout<<"Geo topo susceptibility: "<<smart_print(geo_susc.ave_err())<<endl;
  
  //susceptibility reweighted
  boot_t geo_Q1_rew([](int iconf)->double{return geo_Q[iconf]*chrono_topo.weight[iconf];});
  boot_t geo_Q2_rew([](int iconf)->double{return sqr(geo_Q[iconf])*chrono_topo.weight[iconf];});
  geo_Q1_rew/=chrono_topo.conf_weight;
  geo_Q2_rew/=chrono_topo.conf_weight;
  boot_t geo_susc_rew=(geo_Q2_rew)/(L*L)*1e4;
  cout<<"Geo topo susceptibility rew using average potential not prop err:       "<<smart_print(geo_susc_rew.ave_err())<<endl;
  
  //susceptibility reweighted in slices
  boot_t geo_susc_rew_sliced=geo_susc_rew;
  boot_t geo_susc_rew_slice[npart];
  if(chrono_topo.coeff!=0)
    {
      for(int ipart=0;ipart<npart;ipart++)
	{
	  boot_t geo_Q1_rew_slice([ipart](int iconf)->double{return geo_Q[iconf]*chrono_topo.weight_slice[ipart][iconf];});
	  boot_t geo_Q2_rew_slice([ipart](int iconf)->double{return sqr(geo_Q[iconf])*chrono_topo.weight_slice[ipart][iconf];});
	  geo_Q1_rew_slice/=chrono_topo.conf_weight_slice[ipart];
	  geo_Q2_rew_slice/=chrono_topo.conf_weight_slice[ipart];
	  
	  geo_susc_rew_slice[ipart]=(geo_Q2_rew_slice)/(L*L)*1e4;
	  
	  cout<<" Geo topo susceptibility using potential of slice "<<ipart<<": "<<smart_print(geo_susc_rew_slice[ipart].ave_err())<<endl;
	}
      cout<<"Geo topo susceptibility rew propagating error on potential:         "<<smart_print(slice_ave_err(geo_susc_rew_slice))<<endl;
      
      boot_t geo_Q1_rew_sliced([](int iconf)->double{return geo_Q[iconf]*chrono_topo.weight_slice[ipart_fun(iconf)][iconf];});
      boot_t geo_Q2_rew_sliced([](int iconf)->double{return sqr(geo_Q[iconf])*chrono_topo.weight_slice[ipart_fun(iconf)][iconf];});
      boot_t conf_weight_sliced([](int iconf)->double{return chrono_topo.weight_slice[ipart_fun(iconf)][iconf];});
      geo_Q1_rew_sliced/=conf_weight_sliced;
      geo_Q2_rew_sliced/=conf_weight_sliced;
      geo_susc_rew_sliced=(geo_Q2_rew_sliced-sqr(geo_Q1_rew_sliced))/(L*L)*1e4;
      
      cout<<"Geo topo susceptibility rew weighting each part with its potential: "<<smart_print(geo_susc_rew_sliced.ave_err())<<endl;
    }
  
  for(int ilev=0;ilev<=nstout_lev;ilev++)
    {
      boot_t geo_Q_non_geo_Q([&](int iconf)->double{return geo_Q[iconf]*non_geo_Q[ilev][iconf];});
      ave_err_t reno=(geo_Q_non_geo_Q/geo_Q2).ave_err();
      Q_reno[ilev]=reno.first;
      cout<<"Renormalization constant["<<ilev<<"]: "<<smart_print(reno)<<"   "<<endl;
    }
  
  //compute noise variance
  for(int ilev=0;ilev<=nstout_lev;ilev++)
    {
      double z=(std::isnan(Q_reno[ilev])?1:Q_reno[ilev]);
      boot_t nonint([z,ilev](int iconf)->double{return non_geo_Q[ilev][iconf]-geo_Q[iconf]*z;});
      boot_t nonint2([z,ilev](int iconf)->double{return sqr(non_geo_Q[ilev][iconf]-geo_Q[iconf]*z);});
      cout<<"Noise stdev["<<ilev<<"]: "<<smart_print(sqrt(nonint2-nonint*nonint).ave_err())<<"   "<<endl;
    }
  
  //compute occupancy per topological sector
  vector<boot_t> occupancy(max_geo_Q-min_geo_Q+1);
  boot_t total_occupancy=0;
  ofstream occupancy_out("plots/occupancy_ps.xmg");
  occupancy_out<<"@type xydy"<<endl;
  for(int iQ=min_geo_Q;iQ<=max_geo_Q;iQ++)
    {
      int jQ=iQ-min_geo_Q;
      occupancy[jQ].populate([iQ](int iconf)->int{return (geo_Q[iconf]==iQ);});
      total_occupancy+=occupancy[jQ];
      auto occupancy_ave_err=occupancy[jQ].ave_err();
      if(!std::isnan(occupancy_ave_err.second)) occupancy_out<<iQ<<" "<<occupancy_ave_err<<endl;
    }
  cout<<"Check, total occupancy: "<<smart_print(total_occupancy.ave_err())<<endl;
  
  ///compute the observables sector per sector
  boot_t ene_ps;
  vector<boot_t> mag_ps(max_geo_Q-min_geo_Q+1);
  vector<boot_t> xi_g_ps(max_geo_Q-min_geo_Q+1);
  ofstream ene_ps_out("plots/ene_ps.xmg");
  ofstream mag_ps_out("plots/mag_ps.xmg");
  ofstream mag_ps_rescaled_out("plots/mag_ps_rescaled.xmg");
  ofstream xi_g_ps_out("plots/xi_g_ps.xmg");
  ofstream out_eff_ps("plots/effective_mass_ps.xmg");
  ene_ps_out<<"@type xydy"<<endl;
  mag_ps_out<<"@type xydy"<<endl;
  mag_ps_rescaled_out<<"@type xydy"<<endl;
  xi_g_ps_out<<"@type xydy"<<endl;
  out_eff_ps<<"@type xydy"<<endl;
  for(int iQ=min_geo_Q;iQ<=max_geo_Q;iQ++)
    {
      int jQ=iQ-min_geo_Q;
      int kQ=-iQ-min_geo_Q;
      
      //compute occupancy
      boot_t occ([iQ](int iconf)->int{return (geo_Q[iconf/compute_corr_each*compute_corr_each]==iQ);});
      
      //compute energy and normalize
      boot_t ene_ps([iQ](int iconf)->double{return ene[iconf]*(geo_Q[iconf]==iQ);});
      ene_ps/=occ;
      auto ene_ps_ave_err=ene_ps.ave_err();
      if(!std::isnan(ene_ps_ave_err.second)) ene_ps_out<<iQ<<" "<<ene_ps_ave_err<<endl;
      
      //compute correlation function
      vector<boot_t> correlation_ps(L/2+1);
      for(int t=0;t<L/2+1;t++)
	{
	  correlation_ps[t].populate([iQ,t](int iconf)->double{return corr[t][iconf/compute_corr_each]*
		(geo_Q[iconf/compute_corr_each*compute_corr_each]==iQ);});
	  correlation_ps[t]/=occ;
	}
      auto effective_mass_ps=effective_mass(correlation_ps);
      for(int t=0;t<L/2;t++)
	{
	  auto effective_mass_ave_err=effective_mass_ps[t].ave_err();
	  if(!std::isnan(effective_mass_ave_err.second)) out_eff_ps<<t<<" "<<effective_mass_ave_err<<endl;
	}
      out_eff_ps<<"&"<<endl;
      
      //compute magnetization and normalize
      boot_t xi_g_num([iQ](int iconf)->double{return mag[0][iconf/compute_corr_each]*
	    (geo_Q[iconf/compute_corr_each*compute_corr_each]==iQ);});
      mag_ps[jQ]=xi_g_num/occ;
      auto mag_ps_ave_err=mag_ps[jQ].ave_err();
      if(!std::isnan(mag_ps_ave_err.second)) mag_ps_out<<iQ<<" "<<mag_ps_ave_err<<endl;
      if(iQ>=0 && kQ>=0)
	{
	  auto mag_ps_symm_ave_err=(0.5*(mag_ps[jQ]+mag_ps[kQ])).ave_err();
	  if(!std::isnan(mag_ps_symm_ave_err.second)) mag_ps_rescaled_out<<pow((double)iQ/L,2)<<" "<<mag_ps_symm_ave_err<<endl;
	}
      
      //compute xi_g
      boot_t xi_g_den([iQ](int iconf)->double{return mag[1][iconf/compute_corr_each]*
	    (geo_Q[iconf/compute_corr_each*compute_corr_each]==iQ);});
      xi_g_ps[jQ]=compute_xi_g(xi_g_num,xi_g_den);
      auto xi_g_ps_ave_err=xi_g_ps[jQ].ave_err();
      if(!std::isnan(xi_g_ps_ave_err.second)) xi_g_ps_out<<iQ<<" "<<xi_g_ps_ave_err<<endl;
  }
  
  //plot the noise
  if(nsweep<1000000) plot_noise(Q_reno);
  
  //compute energy
  boot_t energy([](int iconf)->double{return ene[iconf];});
  cout<<"Energy: "<<smart_print(energy.ave_err())<<endl;
  boot_t energy_rew([](int iconf)->double{return ene[iconf]*chrono_topo.weight[iconf];});
  energy_rew/=chrono_topo.conf_weight;
  cout<<"Energy reweighted: "<<smart_print(energy_rew.ave_err())<<endl;
  
  //compute magnetization
  boot_t magnetization([](int iconf)->double{return mag[0][iconf/compute_corr_each];});
  boot_t magnetization_rew([](int iconf)->double{int icorr=iconf/compute_corr_each;return mag[0][icorr]*chrono_topo.weight[icorr*compute_corr_each];});
  
  cout<<"Magnetization: "<<smart_print(magnetization.ave_err())<<endl;
  cout<<"Magnetization_rew: "<<smart_print((magnetization_rew/chrono_topo.corr_weight).ave_err())<<endl;
  
  //compute autocorrelation time
  compute_mag_tint("plots/mag_correlation.xmg",mag[0]);
  
  //compute xi_g
  boot_t xi_g_den([](int iconf)->double{int icorr=iconf/compute_corr_each;return mag[1][icorr];});
  boot_t xi_g_den_rew([](int iconf)->double{int icorr=iconf/compute_corr_each;return mag[1][icorr]*chrono_topo.weight[icorr*compute_corr_each];});
  boot_t xi_g=compute_xi_g(magnetization,xi_g_den);
  boot_t xi_g_rew=compute_xi_g(magnetization_rew,xi_g_den_rew);
  cout<<"Xi_g: "<<smart_print(xi_g.ave_err())<<", L/Xi_g: "<<smart_print((L/xi_g).ave_err())<<", topo_susc*xi_g^2: "<<smart_print((geo_susc*sqr(xi_g)).ave_err())<<endl;
  cout<<"Xi_g_rew: "<<smart_print(xi_g_rew.ave_err())<<", L/Xi_g_rew: "<<smart_print((L/xi_g_rew).ave_err())<<", (topo_susc*xi_g^2)_rew: "<<smart_print((geo_susc_rew_sliced*sqr(xi_g_rew)).ave_err())<<endl;
  if(chrono_topo.coeff!=0)
    {
      auto out=slice_ave_err(geo_susc_rew_slice);
      auto fact=sqr(xi_g_rew).ave_err();
      //cout<<" fact: "<<smart_print(fact)<<" "<<smart_print(out)<<endl;
      out.second=sqrt(sqr(out.first*fact.second)+sqr(out.second*fact.first));
      out.first=geo_susc_rew.ave_err().first*fact.first;
      cout<<" propagating only error: "<<smart_print(out)<<endl;
    }
  
  //compute correlation function
  vector<boot_t> correlation(L/2+1);
  for(int t=0;t<L/2+1;t++) correlation[t].populate([t](int iconf)->double{return corr[t][iconf/compute_corr_each];});
  
  //print the correlation
  ofstream out_corr("plots/correlation.xmg");
  out_corr<<"@type xydy"<<endl;
  for(int t=0;t<L/2+1;t++) out_corr<<t<<" "<<correlation[t].ave_err()<<endl;
  out_corr.close();
  
  //compute and print the effective mass
  auto eff=effective_mass(correlation);
  ofstream out_eff("plots/effective_mass.xmg");
  out_eff<<"@type xydy"<<endl;
  for(int t=0;t<L/2;t++) out_eff<<t<<" "<<eff[t].ave_err()<<endl;
  out_eff.close();
  
  //compute the number of transitions
  int base_trans=0;
  int ntrans=0;
  int ntrue_trans=0;
  int prev_Q=geo_Q[nterm];
  int prev_true_Q=geo_Q[nterm];
  int same_Q=0;
  for(int iconf=nterm;iconf<nsweep;iconf++)
    {
      int cur_Q=geo_Q[iconf];
      
      //if we changed mark down
      if(prev_Q!=cur_Q)
	{
	  ntrans++;
	  same_Q=0;
	}
      
      //increase the time since which we are in the same Q
      same_Q++;
      
      //if we are in the same Q since 10, this is a state
      if(same_Q==10)
	{
	  if(prev_true_Q!=cur_Q) ntrue_trans++;
	  prev_true_Q=cur_Q;
	}
      
      //mark previous step
      prev_Q=cur_Q;
      base_trans++;
    }
  
  double freq_trans=bayes_bin_ave(ntrans,base_trans),err_trans=bayes_bin_err(ntrans,base_trans);
  double freq_true_trans=bayes_bin_ave(ntrue_trans,base_trans),err_true_trans=bayes_bin_err(ntrue_trans,base_trans);
  cout<<"Frequency of transitions: "<<freq_trans<<" "<<err_trans<<" = exp("<<log(freq_trans)<<"["<<err_trans/freq_trans<<"])"<<endl;
  cout<<"Frequency of true transitions: "<<freq_true_trans<<" "<<err_true_trans<<" = exp("<<log(freq_true_trans)<<"["<<err_true_trans/freq_true_trans<<"])"<<endl;
  
  cout<<">>Time: "<<time(0)-init_time<<endl;

  return 0;
}
