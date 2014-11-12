#include "common.cpp"

int main(int narg,char **arg)
{
  if(narg<2) crash("use %s input" ,arg[0]);
  scan_input(arg[1]);
  
  //print the 0 0 two pts
  ofstream two_pts_plot("two_pts_corr.xmg");
  two_pts_plot<<"@type xydy"<<endl;
  two_pts_plot<<g->data[0*nlevls+0]<<endl;

  //load 3pts
  jvec *three_pts=new jvec[nlevls];
  int three_pts_par[3]={+1,-1,+1};
  ofstream three_pts_plot("three_pts_corr.xmg");
  three_pts_plot<<"@type xydy"<<endl;
  for(int iop=0;iop<nlevls;iop++)
    {
      int jop=op_map[iop];
      cout<<"iop: "<<iop<<", jop: "<<jop<<endl;
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
  two_pts_fit(E_V_fit,Z2,VKVK_mv,tfit_VV_min,tfit_VV_max,"VKVK_mv.xmg");  
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
  
  //check if slope in three pts is the good one
  cout<<endl;
  ofstream three_pts_slope_plots("three_pts_slope.xmg");
  jvec three_pts_slope(nlevls,njacks);
  for(int iop=0;iop<nlevls;iop++)
    {
      int tslope_min=5,tslope_max=10;
      jvec corr=numerical_derivative(-log(three_pts[iop]));
      three_pts_slope[iop]=constant_fit(corr,tslope_min,tslope_max);
      three_pts_slope_plots<<write_constant_fit_plot(corr,three_pts_slope[iop],tslope_min,tslope_max,iop*3);
      cout<<"Slope determined from three pts: "<<smart_print(three_pts_slope[iop])
	  <<", expected from two pts: "<<smart_print(E_V-M_P[iop])<<endl;
    }
  
  //check the V-to-P_ground_state matrix element
  //jvec *ground_dT=new jvec[nlevls];
  //for(int iop=0;iop<nlevls;iop++)
  //{
  //ground_dT[iop]=jvec(tsep+1,njacks);
  //for(int t=0;t<=tsep;t++)
  //{
  //ground_dT[iop][t]=(Z_P[iop]*Z_Vs_mv)/(2*E_V*2*M_P[iop])*exp(-E_V*t)*exp(-M_P[iop]*(tsep-t));
  ////ground_dT[iop][t]=VKVK_mv[t]*g->data[iop*nlevls+iop][tsep-t]/(Z_P[iop]*Z_Vs_mv);
  //}
  //}

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
  jvec Q2(nlevls,njacks);
  for(int ilev=0;ilev<nlevls;ilev++)
    {
      jack a2Q2=sqr(M_P_impr[ilev]-E_V)-3*sqr(qi);
      jack a2Q2_max=sqr(M_P_impr[ilev]-M_V);
      Q2[ilev]=a2Q2/sqr(lat_med[ibeta]);
      cout<<"Q2["<<ilev<<"]: "<<smart_print(Q2[ilev])<<" GeV^2="
	  <<smart_print(a2Q2)<<" (latt units), Q2["<<ilev<<"] max: "
	  <<smart_print(a2Q2_max/sqr(lat_med[ibeta]))<<" GeV^2"<<endl;
    }
  
  //cout ratio of massess
  cout<<endl;
  for(int iopt_op=1;iopt_op<nlevls;iopt_op++)
    cout<<"M"<<iopt_op<<"/M0: "<<smart_print(M_P_impr[iopt_op]/M_P_impr[0])<<endl;

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
  
  //compute coefficient for converting mel in ff
  jack glb_coeff_gr=2*M_P[0]/(M_P_impr[0]+M_V)*qi;
  jack glb_coeff_ex=2*M_P_impr[1]/(M_P_impr[1]+M_V)*qi;

  //try the improved one
  jvec *impr_corr=new jvec[nlevls];
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
      
      //divide by correct time dependance and put coeff
      impr_corr[ilev]=-impr_three_pts/impr_dT/((ilev==0)?glb_coeff_gr:glb_coeff_ex);
      
      int tslope_min=3,tslope_max=5;
      jvec three_pts_slope(nlevls,njacks);
      jvec corr_slope=numerical_derivative(-log(fabs(impr_three_pts)));
      three_pts_slope[ilev]=constant_fit(corr_slope,tslope_min,tslope_max);
      cout<<"Slope determined from improved three pts: "<<smart_print(three_pts_slope[ilev])
          <<", expected from two pts: "<<smart_print(E_V-M_P_impr[ilev])<<endl;
    }
  
  //fit the ground state of the simple 3pts
  jvec dT(tsep+1,njacks);
  int ii=1;
  for(int t=0;t<=tsep;t++) dT[t]=(Z_P[ii]*Z_Vs_mv)/(2*E_V*2*M_P[ii])*exp(-E_V*t)*exp(-M_P[ii]*(tsep-t));
  dT.print_to_file("dT.xmg");
  jack fmel_gr_simple=constant_fit(three_pts[ii]/dT/glb_coeff_gr,tfit_3pts_fexc_min,tfit_3pts_fexc_max,"ground_state_mel_simple.xmg");
  cout<<"Ground state matrix element simple: "<<fmel_gr_simple<<endl;
  
  jack fmel_gr=constant_fit(impr_corr[0],tfit_3pts_fexc_min,tfit_3pts_fexc_max,"ground_state_mel.xmg");
  cout<<"Ground state matrix element: "<<fmel_gr<<endl;
  cout<<"Ground state matrix element corrected to Q2=0: "<<fmel_gr*exp(-Q2[0]/(16*sqr(0.54)))<<endl;
  
  jack fmel_ex=constant_fit(impr_corr[1],tfit_3pts_fexc_min,tfit_3pts_fexc_max,"first_state_mel.xmg");
  cout<<"First excited state matrix element: "<<fmel_ex<<endl;
  
  return 0;
}
