#include "convert.cpp"

jvec M(nsemi*2,njacks);
jvec Z(nsemi*2,njacks);

int p[3][3]={{0,0,0},
	     {1,0,0},
	     {1,1,0}};

//time interval for 2pts
int tmin_fit_2pts[nmes];
int tmax_fit_2pts[nmes];
void set_tint()
{
  tmin_fit_2pts[PI]=14;
  tmin_fit_2pts[KA]=14;
  tmin_fit_2pts[D1]=14;
  tmin_fit_2pts[D2]=14;
  tmin_fit_2pts[D3]=14;
  tmin_fit_2pts[D4]=14;
  tmin_fit_2pts[DS1]=17;
  tmin_fit_2pts[DS4]=17;
  
  tmax_fit_2pts[PI]=48;
  tmax_fit_2pts[KA]=48;
  tmax_fit_2pts[D1]=19;
  tmax_fit_2pts[D2]=19;
  tmax_fit_2pts[D3]=19;
  tmax_fit_2pts[D4]=19;
  tmax_fit_2pts[DS1]=30;
  tmax_fit_2pts[DS4]=30;
}

//write the header
string print_header(string name,int itsep,int isemi,int imom,int icur)
{
  ostringstream out;
  out<<"@xaxis  label \"t/tsep("<<tsep_list[itsep]<<")\""<<endl;
  out<<"@title \""<<name<<" "<<semi[isemi].get_full_name()<<"   -   "<<curr_name[icur]<<"    -   ["<<
    p[imom][0]<<","<<p[imom][1]<<","<<p[imom][2]<<"]\""<<endl;
  out<<"@type xydy"<<endl;
  out<<"@s"<<0<<""<<" line type 0"<<endl;
  out<<"@xaxis  label char size 1.600000"<<endl;
  
  return out.str();
}

int main()
{
  set_mesons();
  set_tint();
  
  read_all_corrs();
  
  //fit all mesons where asked
  for(int isemi=0;isemi<nsemi;isemi++)
    for(int i12=1;i12<=2;i12++)
    {
      semi_t s=semi[isemi];
      int imes=s.get_imes(i12);
      jack Z2(njacks);
      int icorr=icorr_2pts(isemi,i12);
      corr_2pts[icorr].print_to_file(combine("plots/corr_2pts_%s_%s.xmg",
					     s.get_name().c_str(),s.get_mes_name(i12).c_str()));
      two_pts_fit(M[icorr],Z2,corr_2pts[icorr],tmin_fit_2pts[imes],tmax_fit_2pts[imes],combine("plots/eff_mass_mes_%s_%s.xmg",
											       s.get_name().c_str(),s.get_mes_name(i12).c_str()));
      Z[icorr]=sqrt(Z2);
    }
  
  for(int isemi=0;isemi<nsemi;isemi++)
    for(int itsep=0;itsep<ntseps;itsep++)
      {
	int tsep=tsep_list[itsep];
	for(int icur=0;icur<ncurrs;icur++)
	  for(int imom=0;imom<nind_mom;imom++)
	    {
	      int ind=icorr_3pts(isemi,itsep,icur,imom);
	      jvec c3=-corr_3pts[ind];
	      
	      ofstream out_rat(combine("plots/ratio_3pts_semi%d_cur%d_tsep%d_imom%d.xmg",isemi,icur,tsep,imom));
	      int imes1=icorr_2pts(isemi,1),imes2=icorr_2pts(isemi,2);
	      jack E1=latt_en(M[imes1],p[imom][0],p[imom][1],p[imom][2]);
	      jack E2=M[imes2];
	      out_rat<<print_header("ratio", itsep,isemi,imom,icur);
	      out_rat<<"@type xydy"<<endl;
	      out_rat<<"@target G0.S"<<ind*3<<endl;
	      out_rat<<"@s"<<ind*3<<" line type 0"<<endl;
	      for(int t=1;t<tsep;t++)
		{
		  //cout<<" imes1 "<<imes1<<" imes2 "<<imes2<<" Z1 "<<Z[imes1]<<" Z2 "<<Z[imes2]<<" E1 "<<E1<<" E2 "<<E2<<" tsep "<<tsep<<" t "<<t<<endl;
		  jack tdep=Z[imes1]*Z[imes2]*exp(-E2*t)*exp(-E1*(tsep-t))/(4*E1*E2);
		  out_rat<<(double)t/tsep<<" "<<c3[t]/tdep<<endl;
		}
	      
	      //////
	      
	      jvec a=aperiodic_effective_mass(c3);
	      ofstream out_eff(combine("plots/eff_mass_3pts_semi%d_cur%d_tsep%d_imom%d.xmg",isemi,icur,tsep,imom));
	      out_eff<<print_header("eff_mass", itsep,isemi,imom,icur);
	      out_eff<<"@target G0.S"<<ind*3<<endl;
	      out_eff<<"@s"<<ind*3<<""<<" line type 0"<<endl;
	      for(int t=1;t<tsep-2;t++) out_eff<<(double)t/tsep<<" "<<a[t]<<endl;
	      out_eff<<"&"<<endl;
	      if(itsep==0) out_eff<<write_constant_with_error(E2-E1,0,1)<<endl;
	    }
      }
  
  
  for(int itsep=0;itsep<ntseps;itsep++)
    {
      int isemi=PI_to_PI;
      int icorr=icorr_2pts(isemi,1);
    //cout<<"Zv: "<<(corr_2pts[icorr_2pts(PI_to_PI,2)][tsep_list[itsep]])<<" "<<
      cout<<(sqr(Z[icorr])*exp(-M[icorr]*tsep_list[itsep])/(2*M[icorr])/(corr_3pts[icorr_3pts(isemi,itsep,V0,0)]))[tsep_list[itsep]/2]<<endl;
    }
  // //prepare te
  // for(int isemi=0;isemi<nsemi;isemi++)
  //   {
  //     semi_t &s=semi[isemi];
  //     int imes1=s.get_imes(1);
  //     int imes2=s.get_imes(2);
      
  //     for(int itsep=0;itsep<ntseps;itsep++)
  // 	{
  // 	  int tsep=tsep_list[itsep];
  // 	  for(int imom=0;imom<nind_mom;imom++)
  // 	    {
  // 	      jvec tdep(tsep+1,njacks);
  // 	      //for(int t=0;t<=tsep;t++)
  // 	      //tdep[t]=Z[imes1]*Z[imes2]*exp()
  // 	    }
  // 	}
  //   }
  // semi_t semi=semi[DS1_to_K];
  
  // jvec mes1=load_meson(semi,1).simmetrized(1);
  // effective_mass(mes1).print_to_file("/tmp/mes1.xmg");
  // jvec mes2=load_meson(semi,2).simmetrized(1);
  // effective_mass(mes2).print_to_file("/tmp/mes2.xmg");
  
  // (effective_mass(mes2)-effective_mass(mes1)).print_to_file("/tmp/mes_diff.xmg");
  
  // for(int itsep=0;itsep<ntseps;itsep++)
  //   {
  //     int tsep=tsep_list[itsep];
      
  //     vector<jvec> cS0=load_0(semi,tsep,raw_S0);
  //     vector<jvec> cV0=load_0(semi,tsep,raw_V0);
  //     vector<jvec> cVK=load_K(semi,tsep,raw_VX,raw_VY,raw_VZ);
  //     vector<jvec> cTK=load_K(semi,tsep,raw_TX,raw_TY,raw_TZ);
      
  //     for(int imom=0;imom<nind_mom;imom++)
  // 	{
  // 	  ofstream out("/tmp/VK_mom_"+to_string(imom)+"_tsep_"+to_string(tsep)+".xmg");
  // 	  out<<"@type xydy"<<endl;
	  
  // 	  jvec y=/*cS0[imom]/*/cV0[imom];
  // 	  for(int t=0;t<tsep;t++)
  // 	    out<<(double)t/tsep<<" "<<y[t]<<endl;
  // 	  out<<endl;
  // 	}
  //   }
  
  return 0;
}


