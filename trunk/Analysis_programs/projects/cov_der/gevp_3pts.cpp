#include "common.cpp"

int main(int narg,char **arg)
{
  if(narg<2) crash("use %s input" ,arg[0]);
  scan_input(arg[1]);
  
  //load 3pts
  jvec *three_pts=new jvec[nlevls];
  int three_pts_par[3]={+1,-1,+1};
  ofstream three_pts_plot("three_pts_corr.xmg");
  three_pts_plot<<"@type xydy"<<endl;
  for(int iop=0;iop<nlevls;iop++)
    {
      int jop=op_map[iop];
      three_pts[iop]=jvec_load(infile,T,njacks,11+jop)/2;
      if(tsep==TH) three_pts[iop]=three_pts[iop].simmetrized(three_pts_par[jop]);
      else         three_pts[iop]=three_pts[iop].subset(0,tsep+1);
      
      //add Zv
      three_pts[iop]*=ZV_mom[ibeta];
      
      three_pts_plot<<three_pts[iop]<<"&"<<endl;
    }
  
  //////////////////////////////////// solve gevp ///////////////////////////////////
  
  //resolve gevp
  g->gevp();
  g->check_orthogonality();
  
  //reorder and convert to full eig
  if(reorder) g->reorder_eig();
  g->convert_to_full_eig_ve();

  //print
  if(DEBUG_GEVP)
    {
      cout<<endl;
      cout<<"Checking eigenvalues"<<endl;
      for(int t=0;t<TH;t++)
	{
	  for(int ilev=0;ilev<nlevls;ilev++) cout<<smart_print(g->eig_va[ilev][t])<<" ";
	  cout<<endl;
	}
      cout<<endl;
    }
  
  {
    ofstream out_gevp("gevp.xmg");
    out_gevp<<"@type xydy"<<endl;
    for(int ilev=0;ilev<nlevls;ilev++)
      {
	jvec m=effective_mass(g->eig_va[ilev]);
	for(int t=0;t<TH;t++)
	  if(fabs(m[t].med())<10)
	    out_gevp<<t<<" "<<m[t]<<endl;
	out_gevp<<"&"<<endl;
      }
  }
  
  ////////////////////////////////// tests and fits ///////////////////////////////////
  
  //fit M, E, Z for vector
  cout<<endl;
  jack M_V(njacks),E_V_fit(njacks),Z2(njacks);
  two_pts_fit(M_V,Z2,VKVK_st,tfit_VV_min,tfit_VV_max,"VKVK_st.xmg");
  jack E_V=latt_en(M_V,th[ibeta]);
  jack Z_Vs_st=sqrt(Z2);
  two_pts_fit(E_V,Z2,VKVK_mv,tfit_VV_min,tfit_VV_max,"VKVK_mv.xmg");  
  jack Z_Vs_mv=sqrt(Z2);
  cout<<"Moving V effect: "<<smart_print((Z_Vs_mv/Z_Vs_st-1)*100)<<"%"<<endl;
  
  //test dispertion relation for vector
  cout<<endl;
  double qi=M_PI*th[ibeta]/TH;
  cout<<"Checking dispertion relation"<<
    ", EV(fit)="<<smart_print(E_V_fit)<<
    ", EV(reco)="<<smart_print(E_V)<<endl;
  
  //determine M_Ps and Z_P for non impr operators
  cout<<endl;
  jvec M_P(nlevls,njacks),Z_P(nlevls,njacks);
  ofstream PP_plots("PP.xmg");
  for(int iop=0;iop<nlevls;iop++)
    {
      jack Z2(njacks);
      two_pts_fit(M_P[iop],Z2,g->data[iop*nlevls+iop],tfit_PP_min,tfit_PP_max);
      Z_P[iop]=sqrt(Z2);
      
      //write plots
      PP_plots<<write_constant_fit_plot(effective_mass(g->data[iop*nlevls+iop]),M_P[iop],tfit_PP_min,tfit_PP_max,iop*3);
      cout<<"Pseudoscalar mass: "<<smart_print(M_P[iop])<<endl;
    }
  
  //writing q2 max
  jack q2_max=sqr(M_V-M_P[1]);
  cout<<endl;
  cout<<"Q2_max: "<<smart_print(q2_max/sqr(lat_med[ibeta]))<<" GeV^2"<<endl;
  
  //print R(J/Psi)
  cout<<endl;
  cout<<"R(J/Psi): "<<smart_print(M_V/M_P[1])<<endl;
  
  //fit P
  cout<<endl;
  jvec M_P_impr(nlevls,njacks),Z_P_impr(nlevls,njacks);
  ofstream PP_impr_plots("PP_impr.xmg");
  for(int iopt_op=0;iopt_op<nlevls;iopt_op++)
    {
      jack Z2(njacks);
      two_pts_fit(M_P_impr[iopt_op],Z2,g->eig_va[iopt_op],tfit_impr_min,tfit_impr_max);
      Z_P_impr[iopt_op]=sqrt(Z2);
      
      //write plots
      PP_impr_plots<<write_constant_fit_plot(effective_mass(g->eig_va[iopt_op]),
					     M_P_impr[iopt_op],tfit_impr_min,tfit_impr_max,iopt_op*3);
      
      cout<<"Pseudoscalar state "<<iopt_op<<" mass: "<<smart_print(M_P_impr[iopt_op])<<endl;
    }
  
  //plot coeff
  for(int ilev=0;ilev<nlevls;ilev++)
    {
      ofstream out_coeff_plot(combine("out_coeff_%d_plot.xmg",ilev).c_str());
      out_coeff_plot<<"@type xydy"<<endl;
      for(int iop=0;iop<nlevls;iop++)
        {
          for(int t=1;t<=TH;t++) out_coeff_plot<<g->eig_ve[iop*nlevls+ilev][t]<<endl;
          out_coeff_plot<<"&"<<endl;
        }
      out_coeff_plot.close();
    }

  //compute Q2
  cout<<endl;
  for(int ilev=0;ilev<nlevls;ilev++)
    {
      double q2=3*sqr(qi);
      jack a2Q2=sqr(M_P_impr[ilev]-E_V)-q2;
      cout<<"Q2["<<ilev<<"]: "<<smart_print(a2Q2/sqr(lat_med[ibeta]))<<" GeV^2, "
	  <<smart_print(a2Q2)<<"in latt units"<<endl;
    }
  
  //cout ratio of massess
  cout<<endl;
  for(int iopt_op=1;iopt_op<nlevls;iopt_op++)
    cout<<"M"<<iopt_op<<"/M0: "<<smart_print(M_P_impr[iopt_op]/M_P_impr[0])<<endl;

  //check if slope in three pts is the good one
  cout<<endl;
  ofstream three_pts_slope_plots("three_pts_slope.xmg");
  jvec three_pts_slope(nlevls,njacks);
  for(int iop=0;iop<nlevls;iop++)
    {
      int tslope_min=3,tslope_max=5;
      jvec corr=numerical_derivative(-log(three_pts[iop]));
      three_pts_slope[iop]=constant_fit(corr,tslope_min,tslope_max);
      three_pts_slope_plots<<write_constant_fit_plot(corr,three_pts_slope[iop],tslope_min,tslope_max,iop*3);
      cout<<"Slope determined from three pts: "<<smart_print(three_pts_slope[iop])
	  <<", expected from two pts: "<<smart_print(E_V-M_P[iop])<<endl;
    }
  
  //check the V-to-P_ground_state matrix element
  jack glb_coeff=2*M_V/(M_P[0]+M_V)*qi;
  ofstream ground_state_me_plots("ground_state_me_plots.xmg");
  ofstream ground_state_dt_plots("ground_state_dt_plots.xmg");
  ground_state_me_plots<<"@type xydy"<<endl;
  ground_state_dt_plots<<"@type xydy"<<endl;
  jvec *ground_dT=new jvec[nlevls];
  
  for(int iop=0;iop<nlevls;iop++)
    {
      ground_dT[iop]=jvec(tsep+1,njacks);
      for(int t=0;t<=tsep;t++)
	{
	  ground_dT[iop][t]=(Z_P[iop]*Z_Vs_mv)/(2*E_V*2*M_P[iop])*exp(-E_V*t)*exp(-M_P[iop]*(tsep-t));
	  //ground_dT[iop][t]=VKVK_mv[t]*g->data[iop*nlevls+iop][tsep-t]/(Z_P[iop]*Z_Vs_mv);
	}
      jvec corr=three_pts[iop]/ground_dT[iop]/glb_coeff;
      ground_state_me_plots<<corr<<"&"<<endl;
      ground_state_dt_plots<<ground_dT[iop]<<"&"<<endl;
    }
  
  //write operator contents
  cout<<endl;
  for(int ilev=0;ilev<nlevls;ilev++)
    {
      cout<<"Content of lev: "<<ilev<<": ";
      jack n2(njacks);
      n2=0;
      for(int iop=0;iop<nlevls;iop++) n2+=sqr(g->eig_ve[iop*nlevls+ilev][tfit_op]);
      for(int iop=0;iop<nlevls;iop++) cout<<smart_print(sqr(g->eig_ve[iop*nlevls+ilev][tfit_op])/n2)<<" ";
      cout<<endl;
    }

  //try the improved one
  jvec *impr_corr=new jvec[nlevls];
  ofstream opt_state_me_plots("opt_state_me_plots.xmg");
  ofstream test("/tmp/test.xmg");
  opt_state_me_plots<<"@type xydy"<<endl;
  test<<"@type xydy"<<endl;
  for(int ilev=0;ilev<nlevls;ilev++)
    {
      jvec impr_three_pts(tsep+1,njacks),impr_dT(tsep+1,njacks);
      for(int t=0;t<=tsep;t++)
	impr_dT[t]=(Z_P_impr[ilev]*Z_Vs_mv)/(2*E_V*2*M_P_impr[ilev])*
	  exp(-E_V*t)*exp(-M_P_impr[ilev]*(tsep-t));
      
      impr_three_pts=0;
      for(int iop=0;iop<nlevls;iop++)
	{
	  jack w=g->full_eig_ve[iop*nlevls+ilev][tfit_op];
	  impr_three_pts+=three_pts[iop]*w;
	}
      
      impr_corr[ilev]=-impr_three_pts/impr_dT/glb_coeff;
      for(int t=0;t<=tsep;t++)
	if(fabs(impr_corr[ilev][t].med())>impr_corr[ilev][t].err() && fabs(impr_corr[ilev][t].med())<3)
	  opt_state_me_plots<<t<<" "<<impr_corr[ilev][t]<<endl;
      opt_state_me_plots<<"&"<<endl;
      test<<fabs(impr_dT)<<fabs(impr_three_pts)<<"&"<<endl;
      
      int tslope_min=3,tslope_max=5;
      jvec corr_slope=numerical_derivative(-log(fabs(impr_three_pts)));
      three_pts_slope[ilev]=constant_fit(corr_slope,tslope_min,tslope_max);
      cout<<"Slope determined from improved three pts: "<<smart_print(three_pts_slope[ilev])
          <<", expected from two pts: "<<smart_print(E_V-M_P_impr[ilev])<<endl;
    }
  
  jack fmel=constant_fit(impr_corr[1],tfit_3pts_fexc_min,tfit_3pts_fexc_max,"first_state_mel.xmg");
  cout<<"First excited state matrix element: "<<fmel<<endl;
  
  return 0;
}
