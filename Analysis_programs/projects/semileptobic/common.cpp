#include "include.h"

//basic parameters
const int T=96,L=48;
const int njacks=60;
const int nind_mom=3;

//list of time separations
const int ntseps=8;
int tsep_list[ntseps]={22,26,28,30,32,34,36,38};

//quark name and short name
const int nquarks=6;
string qname[nquarks]={"l","s","h_0.51","h_0.57","h_0.63","h_0.69"};
enum Q_T{LI,ST,C1,C2,C3,C4};

//meson in term of constituents
const int nmes=8;
typedef pair<Q_T,Q_T> mes_cont_t;
const char full_mes_name[nmes][114]={"\\xp\\0","K","D\\S1\\N","D\\S2\\N","D\\S3\\N","D\\S4\\N","D\\S4\\N\\ss\\N","D\\S4\\N\\ss\\N"};
enum MES_T{PI,KA,D1,D2,D3,D4,DS1,DS4};
vector<mes_cont_t> mes(nmes);

//get the index of the meson from the quark content
int get_imes_from_quark(Q_T q1,Q_T q2)
{
  int imes=0;
  while(imes<nmes && make_pair(max(q1,q2),min(q1,q2))!=mes[imes]) imes++;
  if(imes>=nmes) crash("imes %d should be < %d asking for %d %d",imes,nmes,q1,q2);
  return imes;
}

//3pts function identifier
struct semi_t
{
  Q_T qdec,qpro,qspec;
  int get_imes(int imes){return get_imes_from_quark(imes==1?qpro:qdec,qspec);}
  // Lsh means D-->k but the K has the momentum
  string get_name()
  {return qname[qspec]+qname[qpro]+qname[qdec];}
  string get_full_name()
  {return (string)full_mes_name[get_imes(2)]+" -> "+full_mes_name[get_imes(1)];}
  string get_mes_name(int i12)
  {return qname[qspec]+(i12==1?qname[qpro]:qname[qdec]);}
  void init (MES_T mdec,MES_T mpro,Q_T spec)
  {
    Q_T q_dec_1=mes[mdec].first;
    Q_T q_dec_2=mes[mdec].second;
    Q_T q_pro_1=mes[mpro].first;
    Q_T q_pro_2=mes[mpro].second;
    
    if(q_dec_1==q_pro_1 && q_dec_1==spec)
      {
	qspec=q_dec_1;
	qdec=q_dec_2;
	qpro=q_pro_2;
      }
    else
      if(q_dec_1==q_pro_2 && q_dec_1==spec)
	{
	  qspec=q_dec_1;
	  qdec=q_dec_2;
	  qpro=q_pro_1;
	}
      else
	if(q_dec_2==q_pro_1 && q_dec_2==spec)
	  {
	    qspec=q_dec_2;
	    qdec=q_dec_1;
	    qpro=q_pro_2;
	  }
	else
	  if(q_dec_2==q_pro_2 && q_dec_2==spec)
	    {
	      qspec=q_dec_2;
	      qdec=q_dec_1;
	      qpro=q_pro_1;
	    }
	  else
	    crash("incompatible mesons (%s,%s) and (%s,%s)",
		  qname[q_dec_1].c_str(),qname[q_dec_2].c_str(),
		  qname[q_pro_2].c_str(),qname[q_pro_2].c_str());
    
    // cout<<"Decripting mes "<<mdec<<" -> "<<mpro<<" qspec: "<<spec<<endl;
    // cout<<" "<<get_name()<<" "<<mdec<<" "<<mpro<<endl;
  }
};

//all three points known
int nsemi_known=10,nsemi=nsemi_known;
enum SEMI_CORR_ID{PI_to_PI,KA_to_PI,KA_to_KA,D1_to_PI,D4_to_PI,DS1_to_KA,DS4_to_KA,D2_to_D1,D3_to_D1,D4_to_D1};
vector<semi_t> semi(nsemi_known);

//set the meson content
void set_mesons()
{
  //describe mesons
  mes[PI]=make_pair(LI,LI);
  mes[KA]=make_pair(ST,LI);
  mes[D1]=make_pair(C1,LI);
  mes[D2]=make_pair(C2,LI);
  mes[C3]=make_pair(C3,LI);
  mes[D4]=make_pair(C4,LI);
  mes[DS1]=make_pair(C1,ST);
  mes[DS4]=make_pair(C4,ST);
  
  //describe 3pts
  semi[PI_to_PI].init(PI,PI,LI);
  semi[KA_to_KA].init(KA,KA,LI);
  semi[KA_to_PI].init(KA,PI,LI);
  semi[D1_to_PI].init(D1,PI,LI);
  semi[D4_to_PI].init(D4,PI,LI);
  semi[DS1_to_KA].init(DS1,KA,ST);
  semi[DS4_to_KA].init(DS4,KA,ST);
  semi[D2_to_D1].init(D2,D1,LI);
  semi[D3_to_D1].init(D3,D1,LI);
  semi[D4_to_D1].init(D4,D1,LI);
}

//currents
int ncurrs=4;
enum CURR_T{S0,VK,V0,TK};
char curr_name[4][10]={"S0","Vi","V0","Ti"};

//2pts function
int icorr_2pts(int isemi,int i12)
{return isemi*2+(i12-1);}
vector<jvec> corr_2pts(icorr_2pts(nsemi-1,3-1)+1);

//return the index of a 3pts function
int icorr_3pts(int isemi,int itsep,int icur,int imom)
{return imom+nind_mom*(icur+ncurrs*(itsep+ntseps*isemi));}
vector<jvec> corr_3pts(icorr_3pts(nsemi-1,ntseps-1,ncurrs-1,nind_mom-1)+1);

//write and read all correlators
const char corr_data_path[]="compressed_corrs";
void write_all_corrs()
{
  FILE *fout=open_file(corr_data_path,"w");
  for(auto &c : corr_3pts) c.write_to_binfile(fout);
  for(auto &c : corr_2pts) c.write_to_binfile(fout);
}
void convert_all_corrs();
void read_all_corrs()
{
  if(file_exists(corr_data_path))
    {
      FILE *fin=open_file(corr_data_path,"r");
      //init all 3pts
      for(int isemi=0;isemi<nsemi;isemi++)
	for(int itsep=0;itsep<ntseps;itsep++)
	  for(int icur=0;icur<ncurrs;icur++)
	    for(int imom=0;imom<nind_mom;imom++)
	      {
		int icorr=icorr_3pts(isemi,itsep,icur,imom);
		corr_3pts[icorr]=jvec(tsep_list[itsep]+1,njacks);
		corr_3pts[icorr].load(fin);
	      }
      //init all 2pts
      for(int isemi=0;isemi<nsemi;isemi++)
	for(int i12=1;i12<=2;i12++)
	  {
	    corr_2pts[icorr_2pts(isemi,i12)]=jvec(T/2+1,njacks);
	    corr_2pts[icorr_2pts(isemi,i12)].load(fin);
	  }
    }
  else convert_all_corrs();
}

template <class T> T latt_en(T m,int ipx,int ipy,int ipz)
{
  double px=2*M_PI/L*ipx;
  double py=2*M_PI/L*ipy;
  double pz=2*M_PI/L*ipz;
  return 2*asinh(sqrt(sqr(sin(px/2))+sqr(sin(py/2))+sqr(sin(pz/2))+sqr(sinh(m/2))));
}
