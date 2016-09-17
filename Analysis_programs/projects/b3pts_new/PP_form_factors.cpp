#include "../b3pts_damir/common.cpp"

#define DECIDE 2

#if DECIDE >= 2
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TBox.h>
#include <termios.h>
#endif

enum disp_rel_t{FITTEN,LATTEN,CONTEN}; //CONTEN only to define q2
const disp_rel_t use_disp_rel=FITTEN;

const int njacks=16;
const int i3_VKP5=2;
const int i3_V0P5=3;

const int i2_P5P5=0;

const int add_minus=1;

int im_spec_2pts,im1_2pts,im2_2pts;
int im_spec_3pts,im1_3pts,im2_3pts;

int tmin_K,tmax_K;
int tmin_H,tmax_H;
int tmin_V0,tmax_V0;
int tmin_VK,tmax_VK;
int tmin_TK,tmax_TK;

struct f_results_t
{
  int ith_D;
  int ith_K;
  jack Q2;
  jack fP;
  jack f0;
  double PK,PD;
  bool fP_def;
};

std::vector<f_results_t> f_results;

int icombo_2pts(int iel,int im1,int ith2,int im2,int r=0)
{return iel+2*(im2+12*(r+2*(im1+12*ith2)));}

int icombo_3pts(int iel,int ith1,int im1,int ith2,int im2)
{return iel+19*(im1+3*(ith1+9*(im2+10*ith2)));}

jvec load_2pts_P5P5(int sme,int im1,int ith2,int im2)
{
  jvec a=jvec_load(combine("../DATA/2pts_30_%02d",sme).c_str(),T,njacks,icombo_2pts(i2_P5P5,im1,4+ith2,im2,0));
  jvec b=jvec_load(combine("../DATA/2pts_30_%02d",sme).c_str(),T,njacks,icombo_2pts(i2_P5P5,im1,4+ith2,im2,1));
  jvec c=jvec_load(combine("../DATA/2pts_30_%02d",sme).c_str(),T,njacks,icombo_2pts(i2_P5P5,im1,4-ith2,im2,0));
  jvec d=jvec_load(combine("../DATA/2pts_30_%02d",sme).c_str(),T,njacks,icombo_2pts(i2_P5P5,im1,4-ith2,im2,1));
  
  return (a+b+(c+d)*add_minus).simmetrized(1)/(2+2*add_minus);
}

jvec load_3pts_V0P5(int im_spec,int ith1,int im1,int ith2,int im2)
{
 jvec a=jvec_load(combine("../DATA/3pts_sp%d_30_30",im_spec).c_str(),T,njacks,icombo_3pts(i3_V0P5,4+ith1,im1,4+ith2,im2));
 jvec b=jvec_load(combine("../DATA/3pts_sp%d_30_30",im_spec).c_str(),T,njacks,icombo_3pts(i3_V0P5,4-ith1,im1,4-ith2,im2));
 jvec temp=-(a+add_minus*b)/(1+add_minus);
 if(tsep==T/2) return temp.simmetrized(-1);
 else return temp.subset(0,tsep+1);
}

jvec load_3pts_VKP5(int im_spec,int ith1,int im1,int ith2,int im2)
{
 jvec a=jvec_load(combine("../DATA/3pts_sp%d_30_30",im_spec).c_str(),T,njacks,icombo_3pts(i3_VKP5,4+ith1,im1,4+ith2,im2));
 jvec b=jvec_load(combine("../DATA/3pts_sp%d_30_30",im_spec).c_str(),T,njacks,icombo_3pts(i3_VKP5,4-ith1,im1,4-ith2,im2));
 jvec temp=-(a-add_minus*b)/(1+add_minus);
 if(tsep==T/2) return temp.simmetrized(1);
 else return temp.subset(0,tsep+1);
}

int sign(int i)
{return (i>=0)?1:-1;}

void compute_momenta(jack &P2,jack &P0,double &PK,jack &Q2,jack &Q0,double &QK,jack &m,int sith_m,double &Pm,jack e,jack &M,int sith_M,double &PM,jack E)
{
  switch(use_disp_rel)
    {
    case FITTEN:
      break;
    case LATTEN:
      e=latt_en(m,th_S0[abs(sith_m)]);
      E=latt_en(M,th_S0[abs(sith_M)]);
      break;
    case CONTEN:
      e=cont_en(m,th_S0[abs(sith_m)]);
      E=cont_en(M,th_S0[abs(sith_M)]);
      break;
    }

  P0=E+e;
  Q0=E-e;

  //we are describing the process B->D
  PM=sign(sith_M)*momentum(th_S0[abs(sith_M)]);
  Pm=sign(sith_m)*momentum(th_S0[abs(sith_m)]); //we do not need to put a minus
  //because B->D(th) = D(th)->B
  
  PK=PM+Pm;
  QK=PM-Pm;
  
  P2=P0*P0-3*PK*PK;
  Q2=Q0*Q0-3*QK*QK;
}

//fP,f0,V0,VK,P2,P0,PK,Q2,Q0,QK,EK[0],ED[0]
void solve_fP_f0(jack &fP,jack &f0,jack &V0,jack &VK,jack &P2,jack &P0,double PK,jack &Q2,jack &Q0,double QK,jack &m,jack &M)
{
  //a global factor 3 is dropped
  jack  delta=P0*QK-Q0*PK;
  jack deltaP=V0*QK-Q0*VK;
  jack deltaM=P0*VK-V0*PK;
  
  fP=deltaP/delta;
  jack fM=deltaM/delta;
  
  f0=fP+fM*Q2/(M*M-m*m);
}

char getsinglekey()
{char c;scanf("%c",&c);return c;}
  
char getkey()
{
  int comboL[6]={27,91,49,59,50,68};
  int comboR[6]={27,91,49,59,50,67};
  int combol[3]={27,91,68};
  int combor[3]={27,91,67};
  
  char k=getsinglekey();
  if(k==comboR[0])
    {
      if(getsinglekey()==comboR[1])
	{
	  char t=getsinglekey();
	  if(t==combol[2]) return 'l';
	  if(t==combor[2]) return 'r';
	  if(t==comboR[2])
	    {
	      if(getsinglekey()==comboR[3])
		if(getsinglekey()==comboR[4])
		  {
		    char r=getsinglekey();
		    if(r==comboR[5]) return 'R';
		    if(r==comboL[5]) return 'L';
		  }
	    }
	}
    }

  return k;
}

#if DECIDE >= 2
TCanvas *tela;
TBox *box,*tip;
TGraphErrors *graph_corr,*graph_help;

void decide_tint(jvec corr,jvec help,int &tmin,int &tmax,const char *title)
{
  static int init_tela=true;

  if(init_tela)
    {
      init_tela=false;
      gSystem->ProcessEvents();
      tela=new TCanvas;
      box=new TBox;
      tip=new TBox;
      graph_help=new TGraphErrors;
      graph_corr=new TGraphErrors;

      tela->SetWindowPosition(100,100);
      tela->SetWindowSize(600,600);
      //tela->Divide(2,1);
  
      //draw the current fitting
      //tela->cd(1);
      graph_corr->SetMarkerStyle(20);
      graph_corr->Draw("AP");
      
      //graph_corr->SetColor(3);
      graph_help->SetLineStyle(2);
      graph_help->Draw("SAME");
  
      box->SetFillColor(2);
      box->Draw("SAME");
      
      tip->SetFillColor(3);
      tip->Draw("SAME");
  
      graph_corr->SetName("corr");
      graph_corr->Draw("P");
    }
  
  graph_corr->SetTitle(title);
  int ip=0;
  for(int i=1;i<corr.nel-1;i++)
    {
      graph_corr->SetPoint(ip,i,corr[i].med());
      graph_corr->SetPointError(ip,0,corr[i].err());

      graph_help->SetPoint(ip,i+0.3,help[i].med());
      graph_help->SetPointError(ip,0,help[i].err());

      ip++;
    }

  struct termios info;
  tcgetattr(0, &info);          /* get current terminal attirbutes; 0 is the file descriptor for stdin */
  info.c_lflag &= ~ICANON;      /* disable canonical mode */
  info.c_cc[VMIN] = 1;          /* wait until at least one keystroke available */
  info.c_cc[VTIME] = 0;         /* no timeout */
  tcsetattr(0, TCSANOW, &info); /* set immediately */
  
  char k;
  do
    {
      jack res=constant_fit(corr,tmin,tmax);
      
      box->SetX1(tmin-0.5);
      box->SetX2(tmax+0.5);
      box->SetY1(res.med()-res.err());
      box->SetY2(res.med()+res.err());
      
      tip->SetX1(tmax+0.2);
      tip->SetX2(tmax+0.5);
      tip->SetY1(res.med()-res.err());
      tip->SetY2(res.med()+res.err());
      
      tela->Modified();
      tela->Update();
   
      k=getkey();
      
      switch(k)
	{
	case 'l':tmin--;tmax--;break;
	case 'L':tmax--;break;
	case 'r':tmin++;tmax++;break;
	case 'R':tmax++;break;
	}
    }
  while(k!='q');

}
#endif

void read_input(const char *path)
{
  FILE *input_file=open_file(path,"r");
  
  read_formatted_from_file_expecting((char*)&im_spec_2pts,input_file,"%d","im_spec");
  im_spec_3pts=im_spec_2pts;
  read_formatted_from_file_expecting((char*)&im1_3pts,input_file,"%d","im_S0");
  im1_2pts=im1_3pts;
  read_formatted_from_file_expecting((char*)&im2_3pts,input_file,"%d","im_S1");
  im2_2pts=im2_3pts+2;
  
  cout<<"Considering "<<ori_meson[im_spec_2pts][im1_3pts]<<" -> "<<prod_meson[im_spec_2pts][min(im2_3pts,1)]<<" transition"<<endl;
  
  read_formatted_from_file_expecting((char*)&tmin_K,input_file,"%d","tint_K");
  read_formatted_from_file((char*)&tmax_K,input_file,"%d","tint_K");

  read_formatted_from_file_expecting((char*)&tmin_H,input_file,"%d","tint_H");
  read_formatted_from_file((char*)&tmax_H,input_file,"%d","tint_K");
  
  read_formatted_from_file_expecting((char*)&tmin_V0,input_file,"%d","tint_V0");
  read_formatted_from_file((char*)&tmax_V0,input_file,"%d","tint_V0");
  
  read_formatted_from_file_expecting((char*)&tmin_VK,input_file,"%d","tint_VK");
  read_formatted_from_file((char*)&tmax_VK,input_file,"%d","tint_VK");
  
  read_formatted_from_file_expecting((char*)&tmin_TK,input_file,"%d","tint_TK");
  read_formatted_from_file((char*)&tmax_TK,input_file,"%d","tint_TK");
  
  fclose(input_file);
}

int main()
{
#if DECIDE >= 2
  TApplication myapp("app",NULL,NULL);
#endif

#if DECIDE >= 1
  ifstream decide_file("decided");
#if DECIDE >= 2
  ofstream new_decide_file("new_decided");
#endif
#endif

  read_set_pars("../data_pars");
  read_input("analysis_pars");
 
  debug_load=0;
  
  //fitting K
  jvec ZK(5,njacks),ZK_loc(5,njacks),EK(5,njacks);
  jvec P5P5_K_00[5],P5P5_K_30[5];
  for(int ith=0;ith<5;ith++)
    {
      P5P5_K_00[ith]=load_2pts_P5P5(00,im_spec_2pts,ith,im1_2pts);
      P5P5_K_30[ith]=load_2pts_P5P5(30,im_spec_2pts,ith,im1_2pts);
      
#if DECIDE >= 1
      //if(!(decide_file>>tmin_K>>tmax_K))
      decide_file>>tmin_K>>tmax_K;
#if DECIDE >=2
      decide_tint(effective_mass(P5P5_K_00[ith]),effective_mass(P5P5_K_30[ith]),
		  tmin_K,tmax_K,combine("K_%0d",ith).c_str());
      new_decide_file<<tmin_K<<" "<<tmax_K<<endl;
#endif
#endif

      two_pts_SL_fit(EK[ith],ZK_loc[ith],ZK[ith],P5P5_K_00[ith],P5P5_K_30[ith],tmin_K,tmax_K,tmin_K,tmax_K,
		     combine("plots/Pi_mass_fit/Sm_Lc_effmass_%d.xmg",ith).c_str(),
		     combine("plots/Pi_mass_fit/Sm_Sm_effmass_%d.xmg",ith).c_str(),
		     combine("plots/Pi_mass_fit/chi2_%d.txt",ith).c_str());
      
      cout<<"EK["<<ith<<"]: "<<smart_print(EK[ith])<<endl;
    }
    
  //fitting D
  jvec ZD(5,njacks),ZD_loc(5,njacks),ED(5,njacks);
  jvec P5P5_D_00[5],P5P5_D_30[5];
  for(int ith=0;ith<5;ith++)
    {
      P5P5_D_00[ith]=load_2pts_P5P5(00,im_spec_2pts,ith,im2_2pts);
      P5P5_D_30[ith]=load_2pts_P5P5(30,im_spec_2pts,ith,im2_2pts);
      
#if DECIDE >= 1
      //if(!(decide_file>>tmin_H>>tmax_H))
      decide_file>>tmin_H>>tmax_H;
#if DECIDE >= 2
      decide_tint(effective_mass(P5P5_D_00[ith]),effective_mass(P5P5_D_30[ith]),
		  tmin_H,tmax_H,combine("D_%0d",ith).c_str());      
      new_decide_file<<tmin_H<<" "<<tmax_H<<endl;
#endif
#endif

      two_pts_SL_fit(ED[ith],ZD_loc[ith],ZD[ith],P5P5_D_00[ith],P5P5_D_30[ith],tmin_H,tmax_H,tmin_H,tmax_H,
		     combine("plots/D_mass_fit/Sm_Lc_%d.xmg",ith).c_str(),
		     combine("plots/D_mass_fit/Sm_Sm_%d.xmg",ith).c_str(),
		     combine("plots/D_mass_fit/chi2_%d.txt",ith).c_str());
      cout<<"ED["<<ith<<"]: "<<smart_print(ED[ith])<<endl;
    }
  
  //write down the dispertion relation
  int subt_stand=0;
  ofstream disp_out("plots/Dispertion_relation.xmg");
  disp_out<<"@type xydy\n@s0 legend \"Data\""<<endl;
  for(int ith=0;ith<5;ith++) disp_out<<3*sqr(momentum(th_S0[ith]))<<" "<<sqr(ED[ith])-sqr(ED[0])*subt_stand<<endl;
  disp_out<<"&\n@s1 legend \"Continuum disp rel\""<<endl;
  for(int ith=0;ith<5;ith++) disp_out<<3*sqr(momentum(th_S0[ith]))<<" "<<
			       sqr(cont_en(ED[0],th_S0[ith]))-sqr(cont_en(ED[0],th_S0[0]))*subt_stand<<endl;
  disp_out<<"&\n@s2 legend \"Lattice disp rel\""<<endl;
  for(int ith=0;ith<5;ith++) disp_out<<3*sqr(momentum(th_S0[ith]))<<" "
				     <<sqr(latt_en(ED[0],th_S0[ith]))-sqr(latt_en(ED[0],th_S0[0]))*subt_stand<<endl;
  disp_out<<"&"<<endl;
  
  //print tables
  ZD.print_to_file("plots/ZD_SM_table.txt");
  ZK.print_to_file("plots/ZPi_SM_table.txt");
  
  //prepare table of momenta header
  ofstream momenta_table("plots/momenta_table.txt");
  momenta_table<<"ith1\tith2\t   p_Pi\t\t   p_D\t\t   P2\t\t   Q2"<<endl;
  momenta_table.precision(5);
  momenta_table<<std::fixed<<std::showpos;
  
  //prepare the matrix element table
  ofstream mel_table("plots/matrix_element_table.txt");
  mel_table<<"ith1\tith2\t   V0\t\t   VK"<<endl;
  
  for(int sith1=-4;sith1<=0;sith1++)
    for(int sith2=-4;sith2<5;sith2++)
      //if(sith1<=0 || sith1!=sith2) //eliminate repeated: -1,-1 -> 1,1
      {
	int ith1=fabs(sith1);
	int ith2=fabs(sith2);
	
	jvec V0P5=load_3pts_V0P5(im_spec_3pts,sith1,im1_3pts,sith2,im2_3pts);
	V0P5.print_to_file("plots/V0P5_correlators/%+d_%+d.xmg",sith1,sith2);
	jvec VKP5=load_3pts_VKP5(im_spec_3pts,sith1,im1_3pts,sith2,im2_3pts);
	VKP5.print_to_file("plots/VKP5_correlators/%+d_%+d.xmg",sith1,sith2);
	jvec dT(tsep+1,njacks);
	jvec dTn(tsep+1,njacks);
	jack EK_l;
	jack ED_l;
	
	switch(use_disp_rel)
	  {
	  case FITTEN:
	    EK_l=EK[ith1]; //smallest error
	    ED_l=ED[ith2];
	    break;
	  case LATTEN:
	  case CONTEN:
	    EK_l=latt_en(EK[0],th_S0[ith1]);
	    ED_l=latt_en(ED[0],th_S0[ith2]);
	    break;
	    /* //not sensible
	      case 2:
	      EK_l=cont_en(EK[0],th_S0[ith1]);
	      ED_l=cont_en(ED[0],th_S0[ith2]);
	      break;
	    */
	  }
	
	//write the analytical behavior recosntructed from 2pts
	for(int t=0;t<=tsep;t++)
	  {
	    //dT[t]=(ZK[ith1]*ZD[ith2])/(2*EK[ith1]*2*ED[ith2])*exp(-EK[ith1]*t)*exp(-ED[ith2]*(tsep-t));
	    dT[t]=(ZK[ith1]*ZD[ith2])/(2*EK_l*2*ED_l)*exp(-EK_l*t)*exp(-ED_l*(tsep-t));
	    //dT[t]=(ZK[0]*ZD[0])/(2*EK_l*2*ED_l)*exp(-EK_l*t)*exp(-ED_l*(tsep-t));
	    
	    dTn[t]=P5P5_K_00[ith1][t]*P5P5_D_00[ith2][tsep-t]/ZK_loc[ith1]/ZD_loc[ith2];
	  }
	dT.print_to_file("plots/time_dependances/%+d_%+d.xmg",sith1,sith2);
	dTn.print_to_file("plots/time_dependances/%+d_%+d_num.xmg",sith1,sith2);
	
#if DECIDE >= 1
	
	//if(!(decide_file>>tmin_V0>>tmax_V0>>tmin_VK>>tmax_VK))
	decide_file>>tmin_V0>>tmax_V0>>tmin_VK>>tmax_VK;
	cout<<ith1<<" "<<ith2<<" "<<tmin_V0<<" "<<tmax_V0<<" "<<tmin_VK<<" "<<tmax_VK<<endl;
#if DECIDE >= 2
	  {
	    string title=combine("D(%+02d) -  #Pi (%+02d)",sith2,sith1);
	    decide_tint(V0P5/dT,V0P5/dTn,tmin_V0,tmax_V0,title.c_str());
	    decide_tint(VKP5/dT,VKP5/dTn,tmin_VK,tmax_VK,title.c_str());
	  }
	new_decide_file<<tmin_V0<<" "<<tmax_V0<<" "<<tmin_VK<<" "<<tmax_VK<<endl;
#endif
#endif

	{
	  ofstream out(combine("plots/V0P5_correlators/eff_mass_%+d_%+d.xmg",sith1,sith2).c_str());
	  out<<"@type xydy"<<endl;
	  out<<aperiodic_effective_mass(V0P5)<<endl;
	  out<<write_constant_with_error(EK_l-ED_l,tmin_V0,tmax_V0);
	}

	{
	  ofstream out(combine("plots/VKP5_correlators/eff_mass_%+d_%+d.xmg",sith1,sith2).c_str());
	  out<<"@type xydy"<<endl;
	  out<<aperiodic_effective_mass(VKP5)<<endl;
	  out<<write_constant_with_error(EK_l-ED_l,tmin_VK,tmax_VK);
	}

	//compute the matrix element
	jack V0=/*Zv[ibeta]*/-constant_fit(-(Zv[ibeta]*V0P5/dT).subset(1,23),tmin_V0,tmax_V0,
				       combine("plots/V0_matr_el/%+d_%+d.xmg",sith1,sith2).c_str(),
				       combine("plots/V0_matr_el/%+d_%+d_chi2.xmg",sith1,sith2).c_str());
	jack VK=/*Zv[ibeta]*/-constant_fit(-(Zv[ibeta]*VKP5/dT).subset(1,23),tmin_VK,tmax_VK,
				       combine("plots/VK_matr_el/%+d_%+d.xmg",sith1,sith2).c_str(),
				       combine("plots/VK_matr_el/%+d_%+d_chi2.xmg",sith1,sith2).c_str());
	mel_table<<sith1<<"\t"<<sith2<<"\t"<<smart_print(V0)<<"\t"<<smart_print(VK)<<endl;	
	
	//compute kinematic factors
	jack Q2(njacks);
	jack P2,P0,Q0;
	jack f0(njacks),fP(njacks);
	double PK,QK;
	double PM,Pm;
	compute_momenta(P2,P0,PK, Q2,Q0,QK, EK[0],sith1,Pm,EK_l, ED[0],sith2,PM,ED_l);
	momenta_table<<sith1<<"\t"<<sith2<<"\t"<<Pm<<"\t"<<PM<<"\t"<<P2.med()<<"\t"<<Q2.med()<<endl;
	
	//solve
	if(ith1==0&&ith2==0) f0=(Q0*V0-3*QK*VK)/(sqr(ED[0])-sqr(EK[0]));
	else solve_fP_f0(fP,f0,V0,VK,P2,P0,PK,Q2,Q0,QK,EK[0],ED[0]);
	
	//push the results
	f_results_t res;
	res.ith_K=ith1;
	res.ith_D=ith2;
	res.Q2=Q2;
	res.fP=fP;
	res.f0=f0;
	res.PD=PM;
	res.PK=Pm;
	res.fP_def=!(ith1==0&&ith2==0);
	
	f_results.push_back(res);
      }
  

  //print
  ofstream out_fP("plots/fP.xmg");
  ofstream out_f0("plots/f0.xmg");
  ofstream *out[2]={&out_fP,&out_f0};
  out_fP<<"@type xydy"<<endl;
  out_f0<<"@type xydy"<<endl;
  //const char legend[4][20]={"\\xp\\0 stand","D stand","Breit frame","Else"};
  const char legend[5][20]={"p=0","p=1","p=2","p=3","p=4"};
  //int colors[4]={1,2,4,15};
  int colors[5]={1,2,3,4,5};
  //for(int ifilt=0;ifilt<4;ifilt++)
  for(int ifilt=0;ifilt<5;ifilt++)
    {
      //write legend
      for(int i=0;i<2;i++)
	{
	  (*(out[i]))<<"@s"<<ifilt<<" line type 0"<<endl;
	  (*(out[i]))<<"@s"<<ifilt<<" symbol "<<ifilt+1<<endl;
	  (*(out[i]))<<"@s"<<ifilt<<" symbol size 0.5"<<endl;
	  (*(out[i]))<<"@s"<<ifilt<<" symbol color "<<colors[ifilt]<<endl;
	  (*(out[i]))<<"@s"<<ifilt<<" symbol linewidth 1.5"<<endl;
	  (*(out[i]))<<"@s"<<ifilt<<" errorbar linewidth 1.5"<<endl;
	  (*(out[i]))<<"@s"<<ifilt<<" errorbar riser linewidth 1.5"<<endl;
	  (*(out[i]))<<"@s"<<ifilt<<" errorbar color "<<colors[ifilt]<<endl;
	  (*(out[i]))<<"@s"<<ifilt<<" legend\""<<legend[ifilt]<<"\""<<endl;
	}
	  
      for(std::vector<f_results_t>::iterator it=f_results.begin();it!=f_results.end();it++)
	{
	  bool incl=false;
	  
	  bool K_stand=(fabs(it->PK)<1e-8);
	  bool D_stand=(fabs(it->PD)<1e-8);
	  bool breit_frame=(fabs(it->PK+it->PD)<1e-8);
	  bool else_frame=(!K_stand)&&(!D_stand)&&(!breit_frame);
	  switch(ifilt)
	    {
	    case 0:incl=K_stand;break;
	    case 1:incl=D_stand;break;
	    case 2:incl=breit_frame;break;
	    case 3:incl=else_frame;break;
	    }
	  
	  incl=(ifilt==(it->ith_K));
	  
	  double lat_med[4]={0.486508,0.422773,0.335339,0.268402};
	  if(incl) out_f0<<it->Q2.med()/*/sqr(lat_med[ibeta])*/<<" "<<it->f0<<endl;
	  if(incl && it->fP_def) out_fP<<it->Q2.med()/*/sqr(lat_med[ibeta])*/<<" "<<it->fP<<endl;
	}
      out_fP<<"&"<<endl;
      out_f0<<"&"<<endl;      
    }

  //print
  {
    ofstream out_fP("plots/fP_pimot.xmg");
    ofstream out_f0("plots/f0_pimot.xmg");
    ofstream *out[2]={&out_fP,&out_f0};
    out_fP<<"@type xydy"<<endl;
    out_f0<<"@type xydy"<<endl;
    const char legend[5][20]={"\\xp\\0 stand","\\xp\\0 move1","\\xp\\0 move2","\\xp\\0 move3","\\xp\\0 move4"};
    int colors[5]={1,2,4,15,14};
    for(int ifilt=0;ifilt<5;ifilt++)
      {
	//write legend
	for(int i=0;i<2;i++)
	  {
	    (*(out[i]))<<"@s"<<ifilt<<" line type 0"<<endl;
	    (*(out[i]))<<"@s"<<ifilt<<" symbol "<<ifilt+1<<endl;
	    (*(out[i]))<<"@s"<<ifilt<<" symbol size 0.5"<<endl;
	    (*(out[i]))<<"@s"<<ifilt<<" symbol color "<<colors[ifilt]<<endl;
	    (*(out[i]))<<"@s"<<ifilt<<" symbol linewidth 1.5"<<endl;
	    (*(out[i]))<<"@s"<<ifilt<<" errorbar linewidth 1.5"<<endl;
	    (*(out[i]))<<"@s"<<ifilt<<" errorbar riser linewidth 1.5"<<endl;
	    (*(out[i]))<<"@s"<<ifilt<<" errorbar color "<<colors[ifilt]<<endl;
	    (*(out[i]))<<"@s"<<ifilt<<" legend\""<<legend[ifilt]<<"\""<<endl;
	  }
	
	for(std::vector<f_results_t>::iterator it=f_results.begin();it!=f_results.end();it++)
	  {
	    bool incl=((it->ith_K)==ifilt);

	    double lat_med[4]={0.486508,0.422773,0.335339,0.268402};
	    if(incl) out_f0<<it->Q2.med()/*/sqr(lat_med[ibeta])*/<<" "<<it->f0<<endl;
	    if(incl && it->fP_def) out_fP<<it->Q2.med()/*/sqr(lat_med[ibeta])*/<<" "<<it->fP<<endl;
	  }
	out_fP<<"&"<<endl;
	out_f0<<"&"<<endl;      
      }
  }
  
#if DECIDE >= 1
  decide_file.close();
#endif

#if DECIDE >= 2
  new_decide_file.close();
  system("mv new_decided decided");
#endif

  return 0;
}
