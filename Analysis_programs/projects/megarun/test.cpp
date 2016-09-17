#include "include.h"

const double a=1/2.3;
const int nmom=19,nave_mom=3,nseps=3;
int T=128,TH=T/2,L=64,tsep_tab[nseps]={16,20,24};
int ncfg=33,njacks=ncfg;

//list of Q2 and momenta
int mtab[nmom][3];
int ave_tab[19],deg_tab[nave_mom];

//indipendent momenta and their degeneracy
int ex_mom_combo[nave_mom][3]={{0,0,0},{0,0,1},{0,1,1}};
int C3pts_sp_deg[nave_mom];

//3pts for Pi-Pi
jvec C3pts_Pi_Pi[nseps][nave_mom];
jvec C3pts_Pi_Pi_sp[nseps][nave_mom];

//3pts for K-K
jvec C3pts_K_K[nseps][nave_mom];
jvec C3pts_K_K_sp[nseps][nave_mom];

//3pts for D-D
jvec C3pts_D_D[nseps][nave_mom];
jvec C3pts_D_D_sp[nseps][nave_mom];

//3pts for Ds-Ds
jvec C3pts_Ds_Ds[nseps][nave_mom];
jvec C3pts_Ds_Ds_sp[nseps][nave_mom];

//3pts for D-Pi
jvec C3pts_D_Pi[nseps][nave_mom];
jvec C3pts_D_Pi_sp[nseps][nave_mom];

//3pts for Ds-K
jvec C3pts_Ds_K[nseps][nave_mom];
jvec C3pts_Ds_K_sp[nseps][nave_mom];

//2pts
jvec C2pts_Pi(T,njacks);
jvec C2pts_K(T,njacks);
jvec C2pts_D(T,njacks);
jvec C2pts_Ds(T,njacks);

//results for 2pts fit of Pi
jack M_Pi(njacks),Z_Pi(njacks);
jvec E_Pi(nave_mom,njacks);
jvec C2pts_Pi_tilde[nave_mom];

//results for 2pts fit of K
jack M_K(njacks),Z_K(njacks);
jvec E_K(nave_mom,njacks);
jvec C2pts_K_tilde[nave_mom];

//results for 2pts fit of D
jack M_D(njacks),Z_D(njacks);
jvec E_D(nave_mom,njacks);
jvec C2pts_D_tilde[nave_mom];

//results for 2pts fit of Ds
jack M_Ds(njacks),Z_Ds(njacks);
jvec E_Ds(nave_mom,njacks);
jvec C2pts_Ds_tilde[nave_mom];

//results for Zv
jvec Zv_Pi(nseps,njacks);
jvec Zv_K(nseps,njacks);
jvec Zv_D(nseps,njacks);
jvec Zv_Ds(nseps,njacks);

//fill the table of momenta
void fill_mtab()
{
  for(int imom=0;imom<nave_mom;imom++) deg_tab[imom]=0;
  
  int imom=0;
  for(int m1=-1;m1<=+1;m1++)
    for(int m2=-1;m2<=+1;m2++)
      for(int m3=-1;m3<=+1;m3++)
	{
	  int m2t=m1*m1+m2*m2+m3*m3;
	  
	  if(m2t<=2)
	    {
	      //cout<<imom<<"   "<<m1<<" "<<m2<<" "<<m3<<"   "<<m2t<<endl;
	      if(imom==nmom) crash("counting too much");
	      
	      mtab[imom][0]=m1;
	      mtab[imom][1]=m2;
	      mtab[imom][2]=m3;
	      
	      ave_tab[imom]=m2t;
	      deg_tab[m2t]++;
	      
	      imom++;
	    }
	}

  for(int imom=0;imom<nave_mom;imom++)
    C3pts_sp_deg[imom]=0;
  for(int mu=1;mu<4;mu++)
    for(int imom=0;imom<nmom;imom++)
      if(mtab[imom][3-mu]) C3pts_sp_deg[ave_tab[imom]]++;

  if(imom!=nmom) crash("not filled the entire table");
}

void load_3pts(jvec out[nseps][nave_mom],jvec out_sp[nseps][nave_mom],const char *in)
{
  FILE *fin;

  fin=open_file(in,"r");

  for(int isep=0;isep<nseps;isep++)
    for(int imom=0;imom<nave_mom;imom++)
      {
	int tsep=tsep_tab[isep];
	
	out[isep][imom]=jvec(tsep+1,njacks);
	out[isep][imom]=0;
	out_sp[isep][imom]=jvec(tsep+1,njacks);
	out_sp[isep][imom]=0;
      }

  //load 3pts
  for(int isep=0;isep<nseps;isep++)
    {
      int tsep=tsep_tab[isep];
      
      for(int mu=0;mu<4;mu++)
	for(int imom=0;imom<nmom;imom++)
	  for(int ri=0;ri<2;ri++)
	    {
	      jvec temp(T,njacks);
	      
	      //read t by t
	      for(int t=0;t<T;t++)
		{
		  int rc=fread(temp[t].data,sizeof(double),njacks,fin);
		  if(rc!=njacks) crash("reading obtained %d",rc);
		}
	      
	      //clusterize
	      temp.clusterize();
	      
	      //only real part for V0 and imaginary of VK
	      if(mu==0 && ri==0) out[isep][ave_tab[imom]]+=temp.subset(0,tsep+1);
	      if(mu!=0 && ri==1) out_sp[isep][ave_tab[imom]]+=temp.subset(0,tsep+1)*mtab[imom][3-mu];
	    }
    
      //normalize
      for(int imom=0;imom<nave_mom;imom++)
	{
	  out[isep][imom]/=deg_tab[imom];
	  out_sp[isep][imom]/=C3pts_sp_deg[imom];
	}
    }
  
  //check EOF
  double dum;
  if(fread(&dum,sizeof(double),1,fin)==1) crash("not reached EOF!");
  fclose(fin);
}

void load_2pts(jvec &out,const char *in)
{
  //read 2pts for pion
  FILE *fin=open_file(in,"r");
  //read T by T
  for(int t=0;t<T;t++)
    {
      int rc=fread(out[t].data,sizeof(double),njacks,fin);
      if(rc!=njacks) crash("reading obtained %d",rc);
    }
  out.clusterize();
  out=out.simmetrized(1);

  //check EOF
  double dum;
  if(fread(&dum,sizeof(double),1,fin)==1) crash("not reached EOF!");
  fclose(fin);
}

void load_corrs()
{
  load_2pts(C2pts_Pi,"data/corr_Z2_l64t128_Ls12_m0.000678_LL_m0.000678_LL_15_15");
  load_2pts(C2pts_K,"data/corr_Z2_l64t128_Ls12_m0.000678_LL_m0.02661_LL_15_15");
  load_2pts(C2pts_D,"data/corr_Z2_l64t128_Ls12_m0.34_LL_m0.000678_LL_15_15");
  load_2pts(C2pts_Ds,"data/corr_Z2_l64t128_Ls12_m0.34_LL_m0.02661_LL_15_15");
  load_3pts(C3pts_Pi_Pi,C3pts_Pi_Pi_sp,"data/3pt_ml_0.000678_ml_0.000678_ml_0.000678_maxmom_2");
  load_3pts(C3pts_K_K,C3pts_K_K_sp,"data/3pt_ml_0.000678_ms_0.02661_ms_0.02661_maxmom_2");
  load_3pts(C3pts_D_D,C3pts_D_D_sp,"data/3pt_ml_0.000678_mh_0.34_mh_0.34_maxmom_2");
  load_3pts(C3pts_Ds_Ds,C3pts_Ds_Ds_sp,"data/3pt_ms_0.02661_mh_0.34_mh_0.34_maxmom_2");
  load_3pts(C3pts_D_Pi,C3pts_D_Pi_sp,"data/3pt_ml_0.000678_ml_0.000678_mh_0.34_maxmom_2");
  load_3pts(C3pts_Ds_K,C3pts_Ds_K_sp,"data/3pt_ms_0.02661_ml_0.000678_mh_0.34_maxmom_2");
}

template <class T> T latt_e(T m,double p0,double p1,double p2)
{return 2*asinh(sqrt(sqr(sin(p0/2))+sqr(sin(p1/2))+sqr(sin(p2/2))+sqr(sinh(m/2))));}
template <class T> T latt_e(T m,int *p)
{return latt_e(m,2*M_PI*p[0]/L,2*M_PI*p[1]/L,2*M_PI*p[2]/L);}

double p2(double p0,double p1,double p2)
{return sqr(p0)+sqr(p1)+sqr(p2);}
double p2(int *p)
{return p2(2*M_PI*p[0]/L,2*M_PI*p[1]/L,2*M_PI*p[2]/L);}

void analysis_2pts(const char *outf,jack &M,jack &Z,jvec &E,jvec *tilde_corr,jvec C2pts,int tmin,int tmax)
{
  jack Z2(njacks);
  two_pts_fit(M,Z2,C2pts,tmin,tmax,combine("%s/eff_mass.xmg",outf).c_str());
  Z=sqrt(Z2);
  
  //reconstruct energies from lattice dispertion relation
  for(int i=0;i<nave_mom;i++)
    {
      E[i]=latt_e(M,ex_mom_combo[i]);
      cout<<"E["<<i<<"]: "<<smart_print(E[i])<<endl;
    }

  //tilded corrs
  for(int i=0;i<nave_mom;i++) tilde_corr[i]=jvec(T,njacks);
  
  //tilded at rest
  for(int t=0;t<=T/2;t++) tilde_corr[0][t]=C2pts[t]-C2pts[T/2]*exp(-M*(T/2-t))/2;
  for(int t=T/2;t<T;t++)  tilde_corr[0][t]=C2pts[T/2]*exp(-M*(t-T/2))/2;
  
  //define the other tilded
  for(int imom=0;imom<nave_mom;imom++)
    for(int t=0;t<T;t++)
      tilde_corr[imom][t]=tilde_corr[0][t]*exp(-(E[imom]-M)*t)*M/E[imom];

}

//call the 2pts analysis for Pi
void analysis_Pi()
{
  cout<<"Pi analysis"<<endl;
  cout<<"-------------"<<endl;
  analysis_2pts("plots/Pi",M_Pi,Z_Pi,E_Pi,C2pts_Pi_tilde,C2pts_Pi,19,TH);
}

//call the 2pts analysis for K
void analysis_K()
{
  cout<<"K analysis"<<endl;
  cout<<"-------------"<<endl;
  analysis_2pts("plots/K",M_K,Z_K,E_K,C2pts_K_tilde,C2pts_K,19,TH);
}

//call the 2pts analysis for D
void analysis_D()
{
  cout<<"D analysis"<<endl;
  cout<<"----------"<<endl;
  analysis_2pts("plots/D",M_D,Z_D,E_D,C2pts_D_tilde,C2pts_D,24,48);
}

//call the 2pts analysis for Ds
void analysis_Ds()
{
  cout<<"Ds analysis"<<endl;
  cout<<"----------"<<endl;
  analysis_2pts("plots/Ds",M_Ds,Z_Ds,E_Ds,C2pts_Ds_tilde,C2pts_Ds,24,48);
}

//determine Zv
void analysis_Zv(const char *outf,jvec &Zv,jvec C3pts[nseps][nave_mom],jvec C2pts_tilde[nave_mom])
{
  for(int isep=0;isep<nseps;isep++)
    {
      int tsep=tsep_tab[isep];
      
      //determine Zv
      ofstream out(combine("%s/tsep_%d/Zv_corr.xmg",outf,tsep));
      out<<"@type xydy"<<endl;
      jvec Zv_corr=C2pts_tilde[0][tsep]/C3pts[isep][0];
      for(int t=2;t<tsep-1;t++) out<<(double)t/pow(tsep,1.03)<<" "<<Zv_corr[t]<<endl;
      out<<"&"<<endl;
      
      //fill Zc
      Zv[isep]=Zv_corr[tsep/2];
    }
}

//determine Zv from Pi-Pi, K-K, D-D and Ds-Ds
void analysis_all_Zv()
{
  cout<<"Computing all Zvs"<<endl;
  cout<<"-----------------"<<endl;
  analysis_Zv("plots/Pi_Pi",Zv_Pi,C3pts_Pi_Pi,C2pts_Pi_tilde);
  analysis_Zv("plots/K_K",Zv_K,C3pts_K_K,C2pts_K_tilde);
  analysis_Zv("plots/D_D",Zv_D,C3pts_D_D,C2pts_D_tilde);
  analysis_Zv("plots/Ds_Ds",Zv_Ds,C3pts_Ds_Ds,C2pts_Ds_tilde);
  for(int isep=0;isep<nseps;isep++)
    {
      int tsep=tsep_tab[isep];
      
      cout<<"Zv_Pi["<<tsep<<"]: "<</*smart_print*/(Zv_Pi[isep])<<endl;
      cout<<"Zv_K["<<tsep<<"]: "<</*smart_print*/(Zv_K[isep])<<endl;
      cout<<"Zv_D["<<tsep<<"]: "<</*smart_print*/(Zv_D[isep])<<endl;
      cout<<"Zv_Ds["<<tsep<<"]: "<</*smart_print*/(Zv_Ds[isep])<<endl;
    }
  
    for(int isep=1;isep<nseps-1;isep++)
    {
      cout<<M_Pi.med()<<" "<</*smart_print*/(Zv_Pi[isep])<<endl;
      cout<<M_K.med()<<" "<</*smart_print*/(Zv_K[isep])<<endl;
      cout<<M_D.med()<<" "<</*smart_print*/(Zv_D[isep])<<endl;
      cout<<M_Ds.med()<<" "<</*smart_print*/(Zv_Ds[isep])<<endl;
    }

}

int band_color[3]={10,7,3};
int point_color[nave_mom]={2,10,14},symb[nave_mom]={2,4,3};
  
//output set properties
string set_set_properties(int iset,int tsep,int off)
{
  ostringstream os;
  
  os<<"@type xydy"<<endl;
  os<<"@s"<<iset+off<<" legend \"t\\ssep\\N="<<tsep<<"\""<<endl;
  os<<"@s"<<iset+off<<" line type 0"<<endl;
  os<<"@s"<<iset+off<<" symbol color "<<point_color[iset]<<endl;
  os<<"@s"<<iset+off<<" symbol "<<symb[iset]<<endl;
  os<<"@s"<<iset+off<<" symbol linewidth 1.5"<<endl;
  os<<"@s"<<iset+off<<" symbol fill pattern 1"<<endl;
  os<<"@s"<<iset+off<<" symbol fill color "<<point_color[iset]<<endl;
  os<<"@s"<<iset+off<<" errorbar linewidth 1.5"<<endl;
  os<<"@s"<<iset+off<<" errorbar color "<<point_color[iset]<<endl;
  os<<"@s"<<iset+off<<" errorbar riser linewidth 1.5"<<endl;
  
  return os.str();
}

//perform a full 3pts analysis
void analysis_3pts(const char *outf,jack &M,jack &ZM,jvec &e,jack &Ze,jvec C3pts[nseps][nave_mom],jvec C3pts_sp[nseps][nave_mom],jvec C2pts_M_tilde,jvec C2pts_m_tilde[nave_mom],jvec Zv,const char *process,double tint[3][2][2],bool deg)
{
  //compute Q2
  jvec Q2(nave_mom,njacks);
  for(int imom=0;imom<nave_mom;imom++)
    Q2[imom]=sqr(e[imom]-M)-p2(ex_mom_combo[imom]);

  //define title of fixed-momentum plots
  char mom_title[nave_mom][100];
  for(int imom=0;imom<nave_mom;imom++)
    snprintf(mom_title[imom],100,"@title \"%s  -   Momentum [%d,%d,%d]\"",process,ex_mom_combo[imom][0],ex_mom_combo[imom][1],ex_mom_combo[imom][2]);

  //string to fix axis legend
  char fix_axis_legend[]=
    "@legend 0.35612745098, 0.351470588235\n"
    "@xaxis label char size 1.250000\n"
    "@yaxis label char size 1.250000\n"
    "@xaxis label \"t/t\\ssep\\N\"";
  
  //open the removal of dt and temporal/spatial ratio
  ofstream out_3pts_remove_dt[2][nave_mom];
  ofstream out_3pts_sp_nonsp[nave_mom];
  for(int imom=0;imom<nave_mom;imom++)
    {
      //open
      for(int st=0;st<2;st++) out_3pts_remove_dt[st][imom].open(combine("%s/3pts_%sremove_dt_mom_%d_all_tseps.xmg",outf,((st==0)?"sp_":""),imom).c_str()); 
      out_3pts_sp_nonsp[imom].open(combine("%s/3pts_sp_nonsp_corr_ratio_mom_%d_all_tseps.xmg",outf,imom).c_str());

      //add the x and yaxis ticks
      ofstream *out[3]={&(out_3pts_remove_dt[0][imom]),&(out_3pts_remove_dt[1][imom]),out_3pts_sp_nonsp+imom};
      for(int i=0;i<3;i++)
	{
	  (*out[i])<<mom_title[imom]<<endl;
	  (*out[i])<<"@xaxis tick major 5"<<endl;
	  (*out[i])<<"@xaxis tick minor ticks 1"<<endl;
	  (*out[i])<<"@yaxis tick major 0.01"<<endl;
	  (*out[i])<<"@yaxis tick minor ticks 1"<<endl;
	  (*out[i])<<fix_axis_legend<<endl;
	}
      
      //print yaxis title and subtitle for remove dt
      for(int st=0;st<2;st++)
	{
	  out_3pts_remove_dt[st][imom]<<"@yaxis label \"C\\s3pts\\N(t)/{Z\\S2\\N exp[-Mt-E(t\\ssep\\N-t)]/(4 E M)}"<<((st==0)?"\\sspat\\N":"\\stime\\N")<<"\""<<endl;
	  out_3pts_remove_dt[st][imom]<<"@subtitle \"Remove time dependence analytically: should plateau where ground state dominate\""<<endl;
	}
      
      //print yaxis tile and subtitle for spatial/time ratio
      out_3pts_sp_nonsp[imom]<<"@yaxis label \"C\\s3pts\\N\\Sspat\\N(t)/C\\s3pts\\N\\Stime\\N(t)\""<<endl;
      out_3pts_sp_nonsp[imom]<<"@subtitle \"Spatial/Time ratio: should plateau to "<<(deg?"2\\xp\\0/[L(E+M)]":"an unknown constant")<<"\""<<endl;

      //write the constant value of spatial/time ratio (if known)
      if(deg)
	{
	  out_3pts_sp_nonsp[imom]<<write_line_with_error(2*M_PI/L/(e[imom]+M),M*0,0.0,1.0,2)<<endl;
	  out_3pts_sp_nonsp[imom]<<"@s0 legend \"2\\xp\\0[L(E+M)]\""<<endl;
	  out_3pts_sp_nonsp[imom]<<"@s0 linewidth 2"<<endl;
	}
    }

  //write the intro for effective masses
  ofstream out_3pts_eff_mass[nave_mom];
  ofstream out_3pts_sp_eff_mass[nave_mom];
  for(int imom=0;imom<nave_mom;imom++)
    {
      //open the files
      out_3pts_eff_mass[imom].open(combine("%s/3pts_eff_mass_mom_%d_all_tseps.xmg",outf,imom).c_str());
      out_3pts_sp_eff_mass[imom].open(combine("%s/3pts_sp_eff_mass_mom_%d_all_tseps.xmg",outf,imom).c_str());

      //define subtitle
      char subtitle[2][200];
      for(int i=0;i<2;i++)
	snprintf(subtitle[i],200,"@subtitle \"3pts %s corr. eff.mass (should plateau to E-M%s where ground state dominates)\"",(i?"spat":"time"),((deg&&(imom==0))?"=0":""));
      
      //fix title and axis
      out_3pts_eff_mass[imom]<<
	fix_axis_legend<<endl<<
	mom_title[imom]<<endl<<
	subtitle[0]<<endl<<
	"@yaxis label \"log[C\\s3pts\\N(t)/C\\s3pts\\N(t+1)]\\stime\\N\""<<endl;
      out_3pts_sp_eff_mass[imom]<<
	fix_axis_legend<<endl<<
	mom_title[imom]<<endl<<
	subtitle[1]<<endl<<
	"@yaxis label \"log[C\\s3pts\\N(t)/C\\s3pts\\N(t+1)]\\sspat\\N\""<<endl;
      
      //write the constant line
      string str=write_line_with_error(M-e[imom],M*0,0.0,1.0,2)+"\n@s0 legend \"2E-M\"\n@s0 linewidth 3";
      out_3pts_eff_mass[imom]<<str<<endl;
      out_3pts_sp_eff_mass[imom]<<str<<endl;
    }
  
  for(int isep=0;isep<nseps;isep++)
    {
      int tsep=tsep_tab[isep];
      
      //study 3pts plateaux
      jack mel[2][nave_mom];
      for(int imom=0;imom<nave_mom;imom++)
	{
	  //print effective mass from 3pts
	  jvec ap=aperiodic_effective_mass(C3pts[isep][imom]);
	  jvec ap_sp=aperiodic_effective_mass(C3pts_sp[isep][imom]);
	  
	  //handles to deal with time and spatial
	  out_3pts_eff_mass[imom]<<set_set_properties(isep,tsep,2/*offset*/);
	  out_3pts_sp_eff_mass[imom]<<set_set_properties(isep,tsep,2/*offset*/);
	  
	  //write the aperiodic effective mass
      	  for(int t=2;t<tsep-2;t++)
	    {
	      if(!std::isnan(ap[t].err()))out_3pts_eff_mass[imom]<<(double)t/tsep<<" "<<ap[t]<<endl;
	      if(imom!=0 && !std::isnan(ap_sp[t].err())) out_3pts_sp_eff_mass[imom]<<(double)t/tsep<<" "<<ap_sp[t]<<endl;
	    }
	  
	  //define dt
	  jvec dt(tsep+1,njacks);
	  for(int t=0;t<=tsep;t++)
	    dt[t]=Ze*ZM/(4*e[imom]*M)*exp(-(M*t+e[imom]*(tsep-t)));
	  
	  //fit to constant and print ratio
	  jvec c[2]={C3pts_sp[isep][imom]/dt,C3pts[isep][imom]/dt};
	  for(int st=0;st<2;st++)
	    {
	      mel[st][imom]=constant_fit(c[st],tint[imom][st][0]*tsep,tint[imom][st][1]*tsep);
	      out_3pts_remove_dt[st][imom]<<write_line_with_error
		(mel[st][imom],M*0,tint[imom][st][0],tint[imom][st][1],2)<<endl;
	      out_3pts_remove_dt[st][imom]<<"@s"<<0+3*isep<<" line color "<<band_color[isep]<<endl;
	      out_3pts_remove_dt[st][imom]<<"@s"<<0+3*isep<<" fill color "<<band_color[isep]<<endl;
	      out_3pts_remove_dt[st][imom]<<"@s"<<0+3*isep<<" fill type 1"<<endl;
	      out_3pts_remove_dt[st][imom]<<"@s"<<1+3*isep<<" line color "<<point_color[isep]<<endl;
	      out_3pts_remove_dt[st][imom]<<"@s"<<1+3*isep<<" linewidth 2"<<endl;
	      
	      out_3pts_remove_dt[st][imom]<<set_set_properties(isep,tsep,2+2*isep/*offset*/);
	      for(int t=2;t<tsep-1;t++)
		out_3pts_remove_dt[st][imom]<<(double)t/tsep<<" "<<c[st][t]<<endl;
	    }
	  
	  //break
	  out_3pts_eff_mass[imom]<<"&"<<endl;
	  out_3pts_sp_eff_mass[imom]<<"&"<<endl;
	  for(int st=0;st<2;st++) out_3pts_remove_dt[st][imom]<<"&"<<endl;
	}
          
      //define ratio
      jvec R[nave_mom],R_sp[nave_mom];
      for(int imom=0;imom<nave_mom;imom++)
	{
	  R[imom]=2*sqrt(e[imom]*M)*sqrt(C3pts[isep][imom]*C3pts[isep][imom].inverted()/(C2pts_m_tilde[imom][tsep]*C2pts_M_tilde[tsep]));
	  R_sp[imom]=2*sqrt(e[imom]*M)*sqrt(C3pts_sp[isep][imom]*C3pts_sp[isep][imom].inverted()/(C2pts_m_tilde[imom][tsep]*C2pts_M_tilde[tsep]));
	  R[imom].print_to_file("%s/tsep_%d/ratio_mom_%d.xmg",outf,tsep,imom);
	  R_sp[imom].print_to_file("%s/tsep_%d/ratio_sp_mom_%d.xmg",outf,tsep,imom);
	}

      //plot spatial-nonspatial ratios
      for(int imom=0;imom<nave_mom;imom++)
	{
	  jvec ratio_sp_nonsp=C3pts_sp[isep][imom]/C3pts[isep][imom];
	  out_3pts_sp_nonsp[imom]<<set_set_properties(isep,tsep,2*deg/*offset*/);
	  for(int t=1;t<tsep-1;t++) out_3pts_sp_nonsp[imom]<<(double)t/tsep<<" "<<ratio_sp_nonsp[t]<<endl;
	  out_3pts_sp_nonsp[imom]<<"&"<<endl;
	}

      //determine form factors
      if(deg)
	{
	  //separately determine from time
	  ofstream ffma_plot(combine("%s/tsep_%d/ffma.xmg",outf,tsep).c_str());
	  ffma_plot<<"@type xydy"<<endl;
	  jvec ffma(nave_mom,njacks);
	  for(int imom=0;imom<nave_mom;imom++)
	    {
	      jack A=R[imom][tsep/2];
	      cout<<imom<<" "<<A<<endl;
	      A=mel[1][imom];
	      cout<<imom<<" "<<A<<endl;
	      cout<<"--"<<endl;
	      ffma[imom]=A*Zv[isep]/(M+e[imom]);
	      ffma_plot<<Q2[imom].med()/a/a<<" "<<ffma[imom]<<endl;
	    }

	  //and from space
	  ofstream ffma_plot_sp(combine("%s/tsep_%d/ffma_sp.xmg",outf,tsep).c_str());
	  ffma_plot_sp<<"@type xydy"<<endl;
	  jvec ffma_sp(nave_mom,njacks);
	  for(int imom=1;imom<nave_mom;imom++)
	    {
	      jack A=R[imom][tsep/2];
	      A=mel[0][imom];
	      ffma_sp[imom]=A*Zv[isep]/(2*M_PI/L);
	      ffma_plot_sp<<Q2[imom].med()/a/a<<" "<<ffma_sp[imom]<<endl;
	    }
	}
      else
	{
	  //open the output for f+ and f0
	  ofstream fP_plot(combine("%s/tsep_%d/fP.xmg",outf,tsep).c_str());
	  ofstream f0_plot(combine("%s/tsep_%d/f0.xmg",outf,tsep).c_str());
	  fP_plot<<"@type xydy"<<endl;
	  f0_plot<<"@type xydy"<<endl;

	  //print the Q2 max
	  jack A=R[0][tsep/2];
	  A=mel[1][0];
	  jack f0_st=A*Zv[isep]/(M+e[0]);
	  f0_plot<<Q2[0].med()/a/a<<" "<<f0_st<<endl;
	  
	  for(int imom=1;imom<nave_mom;imom++)
	    {
	      //kinematic factors
	      jack P0=M+e[imom];
	      jack Q0=M-e[imom];
	      double PK=2*M_PI/L;
	      double QK=-2*M_PI/L;

	      //determine the matrix elements
	      jack V0=R[imom][tsep/2]*Zv[isep];
	      jack VK=R_sp[imom][tsep/2]*Zv[isep];
	      V0=mel[1][imom]*Zv[isep];
	      VK=mel[0][imom]*Zv[isep];
	      
	      //a global factor 3 is dropped
	      jack  delta=P0*QK-Q0*PK;
	      jack deltaP=V0*QK-Q0*VK;
	      jack deltaM=P0*VK-V0*PK;
	      
	      //determine f+ and f-
	      jack fP=deltaP/delta;
	      jack fM=deltaM/delta;
	      
	      //combine to get f0
	      jack f0=fP+fM*Q2[imom]/(M*M-e[0]*e[0]);

	      //print 
	      fP_plot<<Q2[imom].med()/a/a<<" "<<fP<<endl;
	      f0_plot<<Q2[imom].med()/a/a<<" "<<f0<<endl;
	    }
	}
    }
}

void analysis_Pi_Pi()
{
  cout<<"Pi-Pi analysis"<<endl;
  cout<<"------------------"<<endl;
  double tint[3][2][2]={{{.45,.55},{.45,.55}},
			{{.40,.50},{.32,.45}},
			{{.40,.50},{.32,.45}}};
  analysis_3pts("plots/Pi_Pi",M_Pi,Z_Pi,E_Pi,Z_Pi,C3pts_Pi_Pi,C3pts_Pi_Pi_sp,C2pts_Pi_tilde[0],C2pts_Pi_tilde,Zv_Pi,"\\xp\\0->\\xp\\0",tint,true);
}

void analysis_K_K()
{
  cout<<"K-K analysis"<<endl;
  cout<<"------------------"<<endl;
  double tint[3][2][2]={{{.45,.55},{.45,.55}},
			{{.40,.50},{.32,.45}},
			{{.40,.50},{.32,.45}}};
  analysis_3pts("plots/K_K",M_K,Z_K,E_K,Z_K,C3pts_K_K,C3pts_K_K_sp,C2pts_K_tilde[0],C2pts_K_tilde,Zv_K,"K->K",tint,true);
}

void analysis_D_D()
{
  cout<<"D-D analysis"<<endl;
  cout<<"------------"<<endl;
  double tint[3][2][2]={{{.45,.55},{.45,.55}},{{.45,.55},{.45,.55}},{{.45,.55},{.45,.55}}};
  analysis_3pts("plots/D_D",M_D,Z_D,E_D,Z_D,C3pts_D_D,C3pts_D_D_sp,C2pts_D_tilde[0],C2pts_D_tilde,Zv_D,"D->D",tint,true);
}

void analysis_Ds_Ds()
{
  cout<<"Ds-Ds analysis"<<endl;
  cout<<"------------"<<endl;
  double tint[3][2][2]={{{.45,.55},{.45,.55}},{{.45,.55},{.45,.55}},{{.45,.55},{.45,.55}}};
  analysis_3pts("plots/Ds_Ds",M_Ds,Z_Ds,E_Ds,Z_Ds,C3pts_Ds_Ds,C3pts_Ds_Ds_sp,C2pts_Ds_tilde[0],C2pts_Ds_tilde,Zv_Ds,"D\\ss\\N->D\\ss\\N",tint,true);
}

void analysis_D_Pi()
{
  cout<<"D-Pi analysis"<<endl;
  cout<<"---------------"<<endl;
  double tint[3][2][2]={{{.45,.55},{.52,.63}},
			{{.69,.82},{.48,.68}},
			{{.76,.89},{.61,.81}}};
  analysis_3pts("plots/D_Pi",M_D,Z_D,E_Pi,Z_Pi,C3pts_D_Pi,C3pts_D_Pi_sp,C2pts_D_tilde[0],C2pts_Pi_tilde,sqrt(Zv_D*Zv_Pi),"D->\\xp\\0",tint,false);
}

void analysis_Ds_K()
{
  cout<<"Ds-K analysis"<<endl;
  cout<<"---------------"<<endl;
  double tint[3][2][2]={{{.45,.55},{.52,.63}},
			{{.69,.82},{.48,.68}},
			{{.76,.89},{.61,.81}}};
  analysis_3pts("plots/Ds_K",M_Ds,Z_Ds,E_K,Z_K,C3pts_Ds_K,C3pts_Ds_K_sp,C2pts_Ds_tilde[0],C2pts_K_tilde,sqrt(Zv_Ds*Zv_K),"D\\ss\\N->K",tint,false);
}

int main()
{
  debug_fit=0;
  
  fill_mtab();
  load_corrs();
  
  
  analysis_Pi();
  analysis_K();
  analysis_D();
  analysis_Ds();
  
  
  analysis_all_Zv();
  
  
  analysis_Pi_Pi();
  analysis_K_K();
  analysis_D_Pi();
  analysis_Ds_K();
  analysis_D_D();
  analysis_Ds_Ds();
  
  return 0;
}
