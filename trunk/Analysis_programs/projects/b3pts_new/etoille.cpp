#include "../b3pts_damir/common.cpp"

int auto_find_3pts_interval;
int ic_fix=2;
int tmin_3pts=0,tmax_3pts=0;

double Zt_med[4]={0.73,0.74,0.78,0.82};
double Za_med[4]={0.746,0.746,0.772,0.780};
FILE *fin; 

const int i2_P5P5=0;
const int i2_VKVK=1;

struct info_3pts_new
{
  int id;
  PARITY pa;
public:
  info_3pts_new(int id,PARITY pa) : id(id),pa(pa) {}
private:
  info_3pts_new();
};

info_3pts_new iS0P5(0,EVN);
info_3pts_new iS0VK(1,ODD);
info_3pts_new iVKP5(2,EVN);
info_3pts_new iV0P5(3,ODD);
info_3pts_new iTKP5(4,ODD);
info_3pts_new iVJVK(5,ODD);
info_3pts_new iVKVJ(6,ODD);
info_3pts_new iP5VK(7,EVN);
info_3pts_new iAKVK(8,EVN);
info_3pts_new iAJVK(9,EVN); //Horr
info_3pts_new iA0VK(10,ODD);
info_3pts_new iAKV0(11,ODD);
info_3pts_new iA0V0(12,ODD); //Horr
info_3pts_new iTJVK(13,ODD); //Horr
info_3pts_new iTKVJ(14,ODD); //Horr
info_3pts_new iBKVK(15,ODD);
info_3pts_new iBJVK(16,ODD);
info_3pts_new iBKVJ(17,ODD);
info_3pts_new iBKV0(18,EVN);

const int njacks=16;
int tint_Detoille[2][2][2];
int tint_H[2][10][2][2];

int icombo_2pts(int iel,int im1,int ith2,int im2,int r=0)
{return iel+2*(im2+12*(r+2*(im1+12*ith2)));}

int icombo_3pts(int iel,int ith1,int im1,int ith2,int im2)
{return iel+19*(im1+3*(ith1+9*(im2+10*ith2)));}

jvec load_2pts_VKVK(int sme,int im1,int im2,int ith=4)
{
  jvec a=jvec_load(combine("../DATA/2pts_30_%02d",sme).c_str(),T,njacks,icombo_2pts(i2_VKVK,im1,ith,im2,0));
  jvec b=jvec_load(combine("../DATA/2pts_30_%02d",sme).c_str(),T,njacks,icombo_2pts(i2_VKVK,im1,ith,im2,1));
  jvec c=jvec_load(combine("../DATA/2pts_30_%02d",sme).c_str(),T,njacks,icombo_2pts(i2_VKVK,im1,8-ith,im2,0));
  jvec d=jvec_load(combine("../DATA/2pts_30_%02d",sme).c_str(),T,njacks,icombo_2pts(i2_VKVK,im1,8-ith,im2,1));
  
  return (a+b+c+d).simmetrized(1)/4;
}

jvec load_2pts_P5P5(int sme,int im1,int im2,int ith=4)
{
  jvec a=jvec_load(combine("../DATA/2pts_30_%02d",sme).c_str(),T,njacks,icombo_2pts(i2_P5P5,im1,ith,im2,0));
  jvec b=jvec_load(combine("../DATA/2pts_30_%02d",sme).c_str(),T,njacks,icombo_2pts(i2_P5P5,im1,ith,im2,1));
  jvec c=jvec_load(combine("../DATA/2pts_30_%02d",sme).c_str(),T,njacks,icombo_2pts(i2_P5P5,im1,8-ith,im2,0));
  jvec d=jvec_load(combine("../DATA/2pts_30_%02d",sme).c_str(),T,njacks,icombo_2pts(i2_P5P5,im1,8-ith,im2,1));
  
  return (a+b+c+d).simmetrized(1)/4;
}

jvec load_3pts(info_3pts_new &info,int im_spec,int ith1,int im1,int ith2,int im2)
{
  jvec a=jvec_load(combine("../DATA/3pts_sp%d_30_30",im_spec).c_str(),T,njacks,icombo_3pts(info.id,ith1,im1,ith2,im2));
  int pa_map[2]={1,-1};
  if(tsep==T/2 && info.pa!=UNC)return a.simmetrized(pa_map[info.pa]);
  else return a.subset(0,tsep+1);
}

void read_fit_tint(const char *path)
{
  fin=open_file(path,"r");
  
  read_formatted_from_file_expecting((char*)&(tint_Detoille[0][0][0]),fin,"%d","Tint_Detoille_SL");
  read_formatted_from_file((char*)&(tint_Detoille[0][0][1]),fin,"%d","Tint_Detoille_SL");

  read_formatted_from_file_expecting((char*)&(tint_Detoille[0][1][0]),fin,"%d","Tint_Detoille_SS");
  read_formatted_from_file((char*)&(tint_Detoille[0][1][1]),fin,"%d","Tint_Detoille_SS");

  read_formatted_from_file_expecting((char*)&(tint_Detoille[1][0][0]),fin,"%d","Tint_Dsetoille_SL");
  read_formatted_from_file((char*)&(tint_Detoille[1][0][1]),fin,"%d","Tint_Dsetoille_SL");
  
  read_formatted_from_file_expecting((char*)&(tint_Detoille[1][1][0]),fin,"%d","Tint_Dsetoille_SS");
  read_formatted_from_file((char*)&(tint_Detoille[1][1][1]),fin,"%d","Tint_Dsetoille_SS");
  
  
  read_formatted_from_file_expecting((char*)&(tint_H[0][0][0][0]),fin,"%d","Tint_H0_SL");
  read_formatted_from_file((char*)&(tint_H[0][0][0][1]),fin,"%d","Tint_H0_SL");
    
  read_formatted_from_file_expecting((char*)&(tint_H[0][0][1][0]),fin,"%d","Tint_H0_SS");
  read_formatted_from_file((char*)&(tint_H[0][0][1][1]),fin,"%d","Tint_H0_SS");
    
  read_formatted_from_file_expecting((char*)&(tint_H[0][9][0][0]),fin,"%d","Tint_H9_SL");
  read_formatted_from_file((char*)&(tint_H[0][9][0][1]),fin,"%d","Tint_H9_SL");
    
  read_formatted_from_file_expecting((char*)&(tint_H[0][9][1][0]),fin,"%d","Tint_H9_SS");
  read_formatted_from_file((char*)&(tint_H[0][9][1][1]),fin,"%d","Tint_H9_SS");
  
  
  read_formatted_from_file_expecting((char*)&(tint_H[1][0][0][0]),fin,"%d","Tint_Hs0_SL");
  read_formatted_from_file((char*)&(tint_H[1][0][0][1]),fin,"%d","Tint_Hs0_SL");
    
  read_formatted_from_file_expecting((char*)&(tint_H[1][0][1][0]),fin,"%d","Tint_Hs0_SS");
  read_formatted_from_file((char*)&(tint_H[1][0][1][1]),fin,"%d","Tint_Hs0_SS");
    
  read_formatted_from_file_expecting((char*)&(tint_H[1][9][0][0]),fin,"%d","Tint_Hs9_SL");
  read_formatted_from_file((char*)&(tint_H[1][9][0][1]),fin,"%d","Tint_Hs9_SL");
    
  read_formatted_from_file_expecting((char*)&(tint_H[1][9][1][0]),fin,"%d","Tint_Hs9_SS");
  read_formatted_from_file((char*)&(tint_H[1][9][1][1]),fin,"%d","Tint_Hs9_SS");
  
  for(int isp=0;isp<2;isp++)
    for(int ih=1;ih<=8;ih++)
      for(int imm=0;imm<2;imm++)
	for(int iSL=0;iSL<2;iSL++)
	  tint_H[isp][ih][iSL][imm]=
	    (int)(tint_H[isp][0][iSL][imm]+(tint_H[isp][9][iSL][imm]-tint_H[isp][0][iSL][imm])/9.0*ih+0.5);
  
  
  read_formatted_from_file_expecting((char*)&(auto_find_3pts_interval),fin,"%d","AutoFind3PtsInterval");
}

int main(int narg,char **arg)
{
  debug_load=debug_fit=0;
  
  if(narg!=2) crash("Use: %s path",arg[0]);
  
  read_set_pars("../data_pars");
  read_fit_tint("analysis_pars");
  
  //open output graph and flush output file
  ofstream mel_file("plots/mel.xmg");
  mel_file<<"@type xydy"<<endl;
  {ofstream mel_dat("mel.dat");}
  
  jack M[2][5];
  const char tag_isp[3]="ls";
  for(int isp=0;isp<2;isp++)
    for(int ith=4;ith<9;ith++)
      {
	//load 2pts for D
	jvec Detoille_corr_SL=load_2pts_VKVK( 0,isp,ic_fix,ith);
	jvec Detoille_corr_SS=load_2pts_VKVK(30,isp,ic_fix,ith);
	stamp(combine("%s/D%c_etoille_SL_th%d",arg[1],tag_isp[isp],ith-4),Detoille_corr_SL);
	stamp(combine("%s/D%c_etoille_SS_th%d",arg[1],tag_isp[isp],ith-4),Detoille_corr_SS);
	
	//fit the D*
	jack MD(njacks),ZDL(njacks),ZDS(njacks);
	two_pts_SL_fit(MD,ZDL,ZDS,Detoille_corr_SL,Detoille_corr_SS,
		       tint_Detoille[isp][0][0],tint_Detoille[isp][0][1],
		       tint_Detoille[isp][1][0],tint_Detoille[isp][1][1],
		       combine("plots/D%c_th%d_etoille_effmass_SL.xmg",tag_isp[isp],ith-4).c_str(),
		       combine("plots/D%c_th%d_etoille_effmass_SS.xmg",tag_isp[isp],ith-4).c_str());
	cout<<"MD"<<tag_isp[isp]<<"*: "<<smart_print(MD)<<endl;
	cout<<"ZD*sm: "<<smart_print(ZDS)<<endl;
	M[isp][ith-4]=MD;
	
	for(int ih=0;ih<10;ih++)
	  {
	    //load 2pts for H
	    jvec H_corr_SL=load_2pts_P5P5( 0,isp,ic_fix+ih,4);
	    jvec H_corr_SS=load_2pts_P5P5(30,isp,ic_fix+ih,4);
	    stamp(combine("%s/B%c_m%d_SL",arg[1],tag_isp[isp],ih),H_corr_SL);
	    stamp(combine("%s/B%c_m%d_SS",arg[1],tag_isp[isp],ih),H_corr_SS);
	    
	    //fit H
	    jack MH(njacks),ZHL(njacks),ZHS(njacks);
	    two_pts_SL_fit(MH,ZHL,ZHS,H_corr_SL,H_corr_SS,
			   tint_H[isp][ih][0][0],tint_H[isp][ih][0][1],
			   tint_H[isp][ih][1][0],tint_H[isp][ih][1][1],
			   combine("plots/H%c_m%d_effmass_SL.xmg",tag_isp[isp],ih).c_str(),
			   combine("plots/H%c_m%d_effmass_SS.xmg",tag_isp[isp],ih).c_str());
	    cout<<"MH["<<ih<<"]: "<<smart_print(MH)<<endl;
	    cout<<"ZH["<<ih<<"]sm: "<<smart_print(ZHS)<<endl;
	    
	    //define the time-dependance
	    jvec dt_an(tsep+1,njacks);
	    for(int t=0;t<=tsep;t++) dt_an[t]=ZDS*ZHS*exp(-MD*t)*exp(-MH*(tsep-t))/(4*MD*MH);
	    jvec dt_nu=Detoille_corr_SL.subset(0,tsep+1)*H_corr_SL.subset(0,tsep+1).inverted()/(ZHL*ZDL);
	    
	    //load 3pts
	    jvec HD_AK=(load_3pts(iAKVK,isp,ith,ic_fix,4,ih)+load_3pts(iAKVK,isp,8-ith,ic_fix,4,ih))/2;
	    stamp(combine("%s/B%c_m%d_AK_D%c_th%d",arg[1],tag_isp[isp],ih,tag_isp[isp],ith-4),HD_AK);
	    (HD_AK/dt_nu).print_to_file(combine("/tmp/B%c_m%d_AK_D%c_th%d.xmg",tag_isp[isp],ih,tag_isp[isp],ith-4).c_str());
	    jvec HD_A0=(load_3pts(iA0VK,isp,ith,ic_fix,4,ih)-load_3pts(iA0VK,isp,8-ith,ic_fix,4,ih))/2;
	    stamp(combine("%s/B%c_m%d_A0_D%c_th%d",arg[1],tag_isp[isp],ih,tag_isp[isp],ith-4),HD_A0);
	    (HD_A0/dt_nu).print_to_file(combine("/tmp/B%c_m%d_A0_D%c_th%d.xmg",tag_isp[isp],ih,tag_isp[isp],ith-4).c_str());
	    jvec HD_P5=(load_3pts(iP5VK,isp,ith,ic_fix,4,ih)-load_3pts(iP5VK,isp,8-ith,ic_fix,4,ih))/2;
	    stamp(combine("%s/B%c_m%d_P5_D%c_th%d",arg[1],tag_isp[isp],ih,tag_isp[isp],ith-4),HD_P5);
	    (HD_P5/dt_nu).print_to_file(combine("/tmp/B%c_m%d_P5_D%c_th%d.xmg",tag_isp[isp],ih,tag_isp[isp],ith-4).c_str());
	    jvec HD_AJ=(load_3pts(iAJVK,isp,ith,ic_fix,4,ih)-load_3pts(iAJVK,isp,8-ith,ic_fix,4,ih))/2;
	    stamp(combine("%s/B%c_m%d_AJ_D%c_th%d",arg[1],tag_isp[isp],ih,tag_isp[isp],ith-4),HD_AJ);
	    (HD_AJ/dt_nu).print_to_file(combine("/tmp/B%c_m%d_AJ_D%c_th%d.xmg",tag_isp[isp],ih,tag_isp[isp],ith-4).c_str());
	    jvec HD_BK=(load_3pts(iBKVK,isp,ith,ic_fix,4,ih)+load_3pts(iBKVK,isp,8-ith,ic_fix,4,ih))/2;
	    stamp(combine("%s/B%c_m%d_BK_D%c_th%d",arg[1],tag_isp[isp],ih,tag_isp[isp],ith-4),HD_BK);
	    (HD_BK/dt_nu).print_to_file(combine("/tmp/B%c_m%d_BK_D%c_th%d.xmg",tag_isp[isp],ih,tag_isp[isp],ith-4).c_str());
	    jvec HD_BJ=load_3pts(iBJVK,isp,ith,ic_fix,4,ih); //unknown how to combine with the others
	    stamp(combine("%s/B%c_m%d_BJ_D%c_th%d",arg[1],tag_isp[isp],ih,tag_isp[isp],ith-4),HD_BJ);
	    (HD_BJ/dt_nu).print_to_file(combine("/tmp/B%c_m%d_BJ_D%c_th%d.xmg",tag_isp[isp],ih,tag_isp[isp],ith-4).c_str());
	    jvec HD_TJ=load_3pts(iTJVK,isp,ith,ic_fix,4,ih); //unknown how to combine with the others
	    stamp(combine("%s/B%c_m%d_TJ_D%c_th%d",arg[1],tag_isp[isp],ih,tag_isp[isp],ith-4),HD_TJ);
	    (HD_TJ/dt_nu).print_to_file(combine("/tmp/B%c_m%d_TJ_D%c_th%d.xmg",tag_isp[isp],ih,tag_isp[isp],ith-4).c_str());
	    jvec HD_VJ=(load_3pts(iVJVK,isp,ith,ic_fix,4,ih)-load_3pts(iVJVK,isp,8-ith,ic_fix,4,ih)-
			load_3pts(iVKVJ,isp,ith,ic_fix,4,ih)+load_3pts(iVKVJ,isp,8-ith,ic_fix,4,ih))/4;
	    stamp(combine("%s/B%c_m%d_VJ_D%c_th%d",arg[1],tag_isp[isp],ih,tag_isp[isp],ith-4),HD_VJ);
	    (HD_VJ/dt_nu).print_to_file(combine("/tmp/B%c_m%d_VJ_D%c_th%d.xmg",tag_isp[isp],ih,tag_isp[isp],ith-4).c_str());
	    
	    //compute the matrix-element correlator
	    jvec HD_AK_mel_corr_an=HD_AK/dt_an/(2*sqrt(MD*MH))*Za_med[ibeta];
	    jvec HD_AK_mel_corr_nu=HD_AK/dt_nu/(2*sqrt(MD*MH))*Za_med[ibeta];
	    //cout<<"dt: "<<1/dt_an/(2*sqrt(MD*MH))*Za_med[ibeta]<<endl;
	    
	    //fit matrix element
	    jack mel(njacks);
	    auto_find_3pts_interval=1;
	    if(auto_find_3pts_interval)
	      {
		//find pivot point for comparing error
		double err_piv=1e300;
		for(int t=1;t<tsep;t++)
		  err_piv=min(HD_AK_mel_corr_an[t].err(),min(HD_AK_mel_corr_nu[t].err(),err_piv));
		//cout<<"ERR MAX: "<<err_piv<<endl;
		
		//make a fit to constant mixing nu and an
		double min_chi2=1e300;
		
		//merge the two derivatives
		jvec tmp=interleave(HD_AK_mel_corr_an,HD_AK_mel_corr_nu);
		jvec tmp_der=interleave(simmetric_derivative(HD_AK_mel_corr_an),simmetric_derivative(HD_AK_mel_corr_nu));
		for(int imin=1;imin<tsep-3;imin++)
		  for(int imax=imin+3;imax<tsep-1;imax++)
		    {
		      //fit to const
		      jack mel_fit=constant_fit(tmp,2*imin,2*imax);
		      
		      int ndof=-1;
		      double tot_chi2=0;
		      for(int iel=2*imin;iel<2*imax;iel++)
			{
			  //fit corr to const
			  double contr=1.3*sqr((tmp[iel].med()-mel_fit.med())/tmp[iel].err());
			  //fit err to its min
			  double contr_err=1.4*sqr(tmp[iel].err()/err_piv);
			  //fit der to zero
			  double contr_der=sqr(tmp_der[iel].med()/tmp_der[iel].err());
			  tot_chi2+=contr+contr_der+contr_err;
			  ndof+=3;
			}
		      tot_chi2/=pow(ndof,1);
		      
		      //cout<<imin<<" "<<imax<<"  "<<ndof<<" "<<tot/ndof<<endl;
		      //if found new best, mark it
		      if(tot_chi2<min_chi2)
			{
			  min_chi2=tot_chi2;
			  tmin_3pts=imin;
			  tmax_3pts=imax;
			  mel=mel_fit;
			}
		    }
		cerr<<tmin_3pts<<" "<<tmax_3pts<<endl;
	      }
	    else
	      {
		read_formatted_from_file((char*)&(tmin_3pts),fin,"%d","Tint_3Pts");
		read_formatted_from_file((char*)&(tmax_3pts),fin,"%d","Tint_3Pts");
	      }
	    
	    //refit the numerical which looks more stable
	    mel=constant_fit(HD_AK_mel_corr_nu,tmin_3pts,tmax_3pts);
	    cout<<"Matrix element: "<<mel<<endl;
	    
	    //cout<<"Min chi2: "<<ymin<<" "<<ymax<<", "<<min_chi2<<endl;
	    ofstream cfit(combine("plots/HD%c_m%d_th%d_AK_mel.xmg",tag_isp[isp],ih,ith-4).c_str());
	    jvec tmp_an=HD_AK_mel_corr_an;
	    jvec tmp_nu=HD_AK_mel_corr_nu;
	    tmp_an[0]=tmp_an[1];
	    tmp_nu[0]=tmp_nu[1];
	    cfit<<write_constant_fit_plot(tmp_an,mel,tmin_3pts,tmax_3pts,0)<<"&"<<endl;
	    cfit<<write_constant_fit_plot(tmp_nu,mel,tmin_3pts,tmax_3pts,3)<<endl;
	    
	    //write the numerical derivative relative error
	    numerical_derivative(HD_AK_mel_corr_nu).print_rel_err_to_file
	      (combine("plots/HD%c_m%d_th%d_AK_der_nu.xmg",tag_isp[isp],ih,ith-4).c_str());
	    
	    //write to the graph and save in jack format
	    mel_file<<mel<<endl;
	    mel.append_to_binfile("mel.dat");
	    
	    //consider also the case of BK mel
	    jvec HD_BK_mel_corr_an=HD_BK/dt_an/(2*sqrt(MD*MH))*Zt_med[ibeta];
	    jvec HD_BK_mel_corr_nu=HD_BK/dt_nu/(2*sqrt(MD*MH))*Zt_med[ibeta];
	    
	    HD_BK_mel_corr_an.print_to_file(combine("plots/HD%c_m%d_th%d_BK_mel_an.xmg",tag_isp[isp],ih,ith-4).c_str());
	    HD_BK_mel_corr_nu.print_to_file(combine("plots/HD%c_m%d_th%d_BK_mel_nu.xmg",tag_isp[isp],ih,ith-4).c_str());
	    
	    /*
	    (HD_BK/HD_AK*Zt_med[ibeta]/Za_med[ibeta]).print_to_file(combine("plots/HD%c_m%d_th%d_BAK_ratio.xmg",tag_isp[isp],ih,ith-4).c_str());
	    */
	  }
    }
  
  for(int isp=0;isp<1;isp++)
    for(int ith=0;ith<5;ith++)
      {
	auto le=latt_en(M[isp][0],th_S0[ith]);
	auto fe=sqrt(sqr(M[isp][0])+3*sqr(momentum(th_S0[ith])));
	
	//cerr<<th_S0[ith]<<" "<<sqr(momentum(th_S0[ith]))<<" "<<smart_print(M[isp][ith])<<" "<<smart_print(le)<<" "<<smart_print(fe)<<endl;
	jack th_c=sqrt((sqr(M[isp][ith])-sqr(M[isp][0]))/3)*TH/M_PI;
	jack th_l=2*asin(1/sqrt(3)*sqrt(sqr(sinh(M[isp][ith]/2))-sqr(sinh(M[isp][0]/2))))*TH/M_PI;
	cerr<<th_S0[ith]<<" "<<smart_print(th_c)<<" "<<smart_print(th_l)<<endl;
      }
  cerr<<"&"<<endl;
  for(int isp=0;isp<1;isp++)
    for(int ith=0;ith<5;ith++)
      cerr<<sqr(momentum(th_S0[ith]))<<" "<<latt_en(M[isp][0],th_S0[ith])<<endl;
  cerr<<"&"<<endl;
  for(int isp=0;isp<1;isp++)
    for(int ith=0;ith<5;ith++)
      cerr<<sqr(momentum(th_S0[ith]))<<" "<<cont_en(M[isp][0],th_S0[ith])<<endl;
  
  return 0;
}
