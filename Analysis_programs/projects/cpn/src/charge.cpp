#include "lib.hpp"

vector<double> charge_real;
vector<double> charge_imag;

vector<vector<double>> non_geo_Q;
vector<int> geo_Q;

void read_charge()
{
  cout<<"Reading charge"<<endl;
  
  //allocate data
  chrono_charge.setup();
  charge_real.resize(nsweep);
  charge_imag.resize(nsweep);
  
  FILE *fintop=open_file("charge","r");
  
  for(int isweep=0;isweep<nsweep;isweep++)
    {
      int iconf;
      char siconf[16],sre[20],sim[20];
      if(fscanf(fintop,"%s %s %s",siconf,sre,sim)==0) crash("reading charge at isweep %d",isweep);
      iconf=atoi(siconf);
      double re=strtod(sre,NULL);
      double im=strtod(sim,NULL);
      
      if(isweep!=iconf) crash("isweep=%d!=iconf=%d",isweep,iconf);
      
      charge_real[isweep]=re;
      charge_imag[isweep]=chrono_charge.coll[isweep]=im;
    }
  fclose(fintop);
}

//read topology
void read_topology()
{
  cout<<"Reading topology"<<endl;
  
  //allocate data
  non_geo_Q.resize(nstout_lev);
  for(int ilev=0;ilev<=nstout_lev;ilev++) non_geo_Q[ilev].resize(nsweep);
  geo_Q.resize(nsweep);
  
  FILE *fintop=open_file("topology","r");
  
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
  read_charge();
  read_topology();
  
  cout<<">>Reading time: "<<time(0)-init_read_time<<" s"<<endl;
}

//use formulazza to compute observables
vector<boot_t> compute_fully_reweighted_obs(vector<boot_t> &obs_vs_Q_fit_pars,vector<double> *pot_fit_pars,double s)
{
  boot_t a=obs_vs_Q_fit_pars[0];
  boot_t b=obs_vs_Q_fit_pars[2];
  vector<boot_t> out;
  for(int ipart=0;ipart<npart;ipart++)
    {
      boot_t e=pot_fit_pars[ipart][2];
      out.push_back(a-b*(2*e+s*s)/(4*e*e));
    }
  
  return out;
}

int main()
{
  cout.precision(16);
  
  int init_time=time(0);
  read_input(true);
  cout<<"NSweeps: "<<nsweep<<endl;
  
  //read the data and reconstruct reweigthing factor
  read_data();
  
  ofstream out_charge("plots/charge.xmg");
  for(int i=0;i<nsweep;i++) out_charge<<charge_imag[i]<<endl;
  
  //draw random indices
  prepare_boot_sample();
  
  chrono_charge.mega_reweight([](double x){return cos(sinh(charge_pot/L)/g*x);});
  
  boot_t geo_Q1([](int iconf)->double{return geo_Q[iconf];});
  boot_t geo_Q2([](int iconf)->double{return sqr(geo_Q[iconf]);});
  boot_t geo_susc=(geo_Q2-sqr(geo_Q1))/(L*L)*1e4;
  cout<<"Geo topo susceptibility: "<<smart_print(geo_susc.ave_err())<<endl;
  
  //phase
  cout<<"Phase: "<<chrono_charge.conf_weight.ave_err()<<endl;
  
  //compute energy
  boot_t energy_est([](int iconf)->double{return ene[iconf];});
  boot_t energy_rew2([](int iconf)->double{return ene[iconf]*cos(sinh(charge_pot/L)/g*charge_imag[iconf]);});
  cout<<"Energy: "<</*smart_print*/(energy_est.ave_err())<<endl;
  boot_t energy_rew([](int iconf)->double{return ene[iconf]*chrono_charge.weight[iconf];});
  energy_rew/=chrono_charge.conf_weight;
  cout<<"Energy reweighted: "<</*smart_print*/(energy_rew.ave_err())<<endl;
  
  //compute energy
  double c=cosh(charge_pot/L)/g;
  double s=sinh(charge_pot/L)/g;
  boot_t charge_rew([s,c](int iconf)->double{return
	s*charge_real[iconf]*cos(s*charge_imag[iconf])+
	c*charge_imag[iconf]*sin(s*charge_imag[iconf]);});
  boot_t charge_rew_den([s](int iconf)->double{return cos(s*charge_imag[iconf]);});
  cout<<"Charge_rew: "<</*smart_print*/((charge_rew/charge_rew_den).ave_err())<<endl;
  
  //compute magnetization
  boot_t magnetization([](int iconf)->double{return mag[0][iconf/compute_corr_each];});
  boot_t magnetization_rew([](int iconf)->double{int icorr=iconf/compute_corr_each;return mag[0][icorr]*chrono_charge.weight[icorr*compute_corr_each];});
  
  cout<<"Magnetization: "<</*smart_print*/(magnetization.ave_err())<<endl;
  cout<<"Magnetization_rew: "<</*smart_print*/((magnetization_rew/chrono_charge.corr_weight).ave_err())<<endl;
  
  //compute xi_g
  boot_t xi_g_den([](int iconf)->double{int icorr=iconf/compute_corr_each;return mag[1][icorr];});
  boot_t xi_g_den_rew([](int iconf)->double{int icorr=iconf/compute_corr_each;return mag[1][icorr]*chrono_charge.weight[icorr*compute_corr_each];});
  boot_t xi_g=compute_xi_g(magnetization,xi_g_den);
  boot_t xi_g_rew=compute_xi_g(magnetization_rew,xi_g_den_rew);
  cout<<"Xi_g: "<</*smart_print*/(xi_g.ave_err())<<", L/Xi_g: "<<smart_print((L/xi_g).ave_err())<<endl;
  cout<<"Xi_g_rew: "<</*smart_print*/(xi_g_rew.ave_err())<<", L/Xi_g_rew: "<<smart_print((L/xi_g_rew).ave_err())<<endl;
  
  //compute correlation function
  vector<boot_t> correlation(L/2+1);
  for(int t=0;t<L/2+1;t++) correlation[t].populate([t](int iconf)->double{return corr[t][iconf/compute_corr_each];});
  
  ////////////////////////////////////////// fit potential ///////////////////////////////////////////
  
  int par_fit_ord=2;
  double margin=0.45;
  double margfin=0.9;
  
  //prepare x
  vector<double> grid_in(chrono_charge.ngrid+1);
  for(int i=0;i<=chrono_charge.ngrid;i++) grid_in[i]=-chrono_charge.barr+i*chrono_charge.width;
  
  //fit
  cout<<"===============Potential fit of ord "<<par_fit_ord<<", margin "<<margin<<"=================="<<endl;
  vector<double> pot_fit_pars[npart];
  vector<double> chi2(npart);
  int ndof;
  for(int ipart=0;ipart<npart;ipart++)
    {
	auto fit_output=poly_fit(grid_in,chrono_charge.pote_slice[ipart],chrono_charge.pote_err,par_fit_ord,chrono_charge.barr*margin,+chrono_charge.barr*margfin);
	pot_fit_pars[ipart]=fit_output.coeffs;
	chi2[ipart]=fit_output.chi2;
	ndof=fit_output.ndof;
    }
  auto ave_err_chi2=ave_err_part([chi2](double ipart)->double{return chi2[ipart];});
  cout<<"Chi2 = "<<smart_print(ave_err_chi2)<<"/"<<ndof<<" = "<<ave_err_chi2.first/ndof<<endl;
  
  vector<ave_err_t> pot_fit_pars_ave_err(par_fit_ord+1);
  cout<<"Parameters"<<endl;
  for(int ipar=0;ipar<=par_fit_ord;ipar++)
    {
      pot_fit_pars_ave_err[ipar]=ave_err_part([pot_fit_pars,ipar](int ipart)->double{return pot_fit_pars[ipart][ipar];});
      cout<<ipar<<" "<<pot_fit_pars_ave_err[ipar]<<endl;
    }
  
  ofstream out_fitted_pot("fitted_potential.xmg");
  out_fitted_pot.precision(16);
  out_fitted_pot<<"@type xydy"<<endl;
  for(int i=0;i<=chrono_charge.ngrid;i++) out_fitted_pot<<grid_in[i]<<" "<<chrono_charge.pote_ave[i]-chrono_charge.pote_ave[chrono_charge.ngrid/2]<<
					    " "<<chrono_charge.pote_err[i]<<endl;
  out_fitted_pot<<"&\n@type xy"<<endl;
  for(int ord=-1;ord<=+1;ord+=2)
    for(double x=-chrono_charge.barr*3;x<+chrono_charge.barr*3;x+=0.1)
      {
	auto out=ave_err_part([x,pot_fit_pars,par_fit_ord](int ipart)->double{return poly(x,pot_fit_pars[ipart])-pot_fit_pars[ipart][0];});
	out_fitted_pot<<ord*x<<" "<<out.first+ord*out.second<<endl;
      }
  
  ///////////////////////////////////////// bin observables //////////////////////////////////////////
  
  ofstream binned_K("plots/K_vs_Q.xmg");
  ofstream binned_E("plots/E_vs_Q.xmg");
  binned_K.precision(16);
  binned_E.precision(16);
  binned_K<<"@type xydy"<<endl;
  binned_E<<"@type xydy"<<endl;
  vector<double> vs_Q;
  vector<boot_t> K_vs_Q;
  vector<boot_t> E_vs_Q;
  
  for(double b=chrono_charge.barr,w=chrono_charge.width,x=-b;x<+b+w-1e-10;x+=w)
    {
      auto is_incl=[x,w](double c){return ((c>=x-w/2 && c<x+w/2)||(c>=-x-w/2 && c<-x+w/2));};
      boot_t num_K([is_incl](int iconf){return is_incl(chrono_charge.coll[iconf])*charge_real[iconf];});
      boot_t num_E([is_incl](int iconf){return is_incl(chrono_charge.coll[iconf])*ene[iconf];});
      boot_t den([is_incl](int iconf){return is_incl(chrono_charge.coll[iconf]);});
      
      vs_Q.push_back(x);
      K_vs_Q.push_back(num_K/den);
      E_vs_Q.push_back(num_E/den);
      
      binned_K<<x<<" "<<(num_K/den).ave_err()<<endl;
      binned_E<<x<<" "<<(num_E/den).ave_err()<<endl;
    }
  
  vector<boot_t> K_vs_Q_fit_pars=poly_fit(vs_Q,K_vs_Q,2).coeffs;
  vector<boot_t> E_vs_Q_fit_pars=poly_fit(vs_Q,E_vs_Q,2).coeffs;
  
  //parameters for K
  cout<<"Fitted K vs Q:"<<endl;
  for(size_t i=0;i<K_vs_Q_fit_pars.size();i++) cout<<i<<" "<<K_vs_Q_fit_pars[i].ave_err()<<endl;
  
  //parameters for E
  cout<<"Fitted E vs Q:"<<endl;
  for(size_t i=0;i<E_vs_Q_fit_pars.size();i++) cout<<i<<" "<<E_vs_Q_fit_pars[i].ave_err()<<endl;
  
  //compute charge
  {
    boot_t a=K_vs_Q_fit_pars[0];
    boot_t b=K_vs_Q_fit_pars[2];
    cout<<"Average charge:"<<endl;
    for(int ipart=0;ipart<npart;ipart++)
      {
	boot_t e=pot_fit_pars[ipart][2];
	boot_t Qave=-(s*(2*e*(b+c-2*a*e)+b*s*s))/(4*e*e);
	cout<<ipart<<" "<<Qave.ave_err()<<endl;
      }
  }
  
  //compute energy
  vector<boot_t> Eave=compute_fully_reweighted_obs(E_vs_Q_fit_pars,pot_fit_pars,s);
  cout<<"Average energy:"<<endl;
  for(int ipart=0;ipart<npart;ipart++) cout<<ipart<<" "<<Eave[ipart].ave_err()<<endl;
  
  return 0;
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  ofstream fout("plots/reco_weight.xmg");
  vector<double> out_x,out_y[npart];
  double norm[npart];
  for(int islice=0;islice<npart;islice++)
    {
      norm[islice]=0;
      double dx=chrono_charge.width/10;
      double tmax=*max_element(chrono_charge.pote_slice[islice].begin(),chrono_charge.pote_slice[islice].end());
      for(double x=-chrono_charge.barr-3*chrono_charge.width;x<chrono_charge.barr+chrono_charge.width*3+0.001;x+=dx)
	{
	  double w=tmax-(chrono_charge.interpolate_potential(x,chrono_charge.pote_slice[islice])+chrono_charge.interpolate_potential(-x,chrono_charge.pote_slice[islice]))/2;
	  double temp=exp(-w)*(cos(sinh(charge_pot/L)*x));
	  if(islice==0) out_x.push_back(x);
	  out_y[islice].push_back(temp);
	  norm[islice]+=exp(-w)*dx;
	}
      transform(out_y[islice].begin(),out_y[islice].end(),out_y[islice].begin(),[=](double in){return in/norm[islice];});
      
      for(size_t i=0;i<out_x.size();i++) fout<<out_x[i]<<" "<<out_y[islice][i]<<endl;
      fout<<"&"<<endl;
      
      //cout<<"islice "<<islice<<" norm: "<<norm[islice]<<endl;
    }
  
  {
    vector<double> out_y_ave(out_x.size()),out_y_err(out_x.size());
    
    reconstruct_band(out_y_ave,out_y_err,out_y,npart);
    vector<double> low(out_x.size()),high(out_x.size());
    for(size_t i=0;i<out_x.size();i++)
      {
	low[i] =out_y_ave[i]-out_y_err[i];
	high[i]=out_y_ave[i]+out_y_err[i];
      }
    
    ofstream fout("plots/reco_weight_band.xmg");
    for(size_t i=0;i<out_x.size();i++) fout<<out_x[i]<<" "<<low[i]<<endl;
    for(size_t i=out_x.size()-1;i<out_x.size();i--) fout<<out_x[i]<<" "<<high[i]<<endl;
    
    auto w_fun=[](double x)
      {
	double p0=chrono_charge.interpolate_potential( 0,chrono_charge.pote_ave);
	double p1=chrono_charge.interpolate_potential(+x,chrono_charge.pote_ave);
	double p2=chrono_charge.interpolate_potential(-x,chrono_charge.pote_ave);
	double p=(p1+p2)/2-p0;
	double w=exp(p)*cos(sinh(charge_pot/L)*x);
	return w;
      };
    boot_t energy([=](int iconf)->double{return ene[iconf];});
    boot_t energy_rew([=](int iconf)->double{return ene[iconf]*w_fun(charge_imag[iconf]);});
    boot_t rew_weight([=](int iconf)->double{return w_fun(charge_imag[iconf]);});
    cout<<smart_print(energy.ave_err())<<" "<<smart_print((energy_rew/rew_weight).ave_err())<<"    "<<smart_print(rew_weight.ave_err())<<endl;
    cout<<charge_pot/L<<" "<<sinh(charge_pot/L)<<endl;
    
    ofstream out_we("/tmp/weight");
    for(int iconf=nterm;iconf<nsweep;iconf++)
      out_we<<iconf<<" "<<charge_imag[iconf]<<" "<<chrono_charge.weight[iconf]<<endl;
    
  }
  
  cout<<">>Time: "<<time(0)-init_time<<endl;
  
  return 0;
}
