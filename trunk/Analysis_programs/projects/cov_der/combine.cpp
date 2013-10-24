#include <cmath>
#include <map>
#include <vector>
#include <string>
#include <string.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <tr1/array>
#include <algorithm>

using namespace std;

//define gamma
const int nga=19;
const char lga[nga][3]={"S0","V1","V2","V3","V0","P5","A1","A2","A3","A0","T1","T2","T3","B1","B2","B3","C1","C2","C3"};
#define V1_ID 1
#define V2_ID 2
#define V3_ID 3
#define P5_ID 5
#define A0_ID 9
#define A0_ID 9
#define T1_ID 10
#define T2_ID 11
#define T3_ID 12
#define B1_ID 13
#define B2_ID 14
#define B3_ID 15
#define C1_ID 16
#define C2_ID 17
#define C3_ID 18

//define the operator
struct op_contr_t
{
  int mu1,mu2,iga;
  double weight;
  op_contr_t(int mu1,int mu2,int iga,double weight) : mu1(mu1),mu2(mu2),iga(iga),weight(weight) {}
private:
  op_contr_t();
};
typedef vector<op_contr_t> op_t;

//correlation specifications
const int NT=48;
template <class T> class corr_t : public tr1::array<T,NT>
{
public:
  corr_t() {}
  corr_t(double a){for(auto it=this->begin();it!=this->end();it++) (*it)=a;}
  ostream &operator<<(ostream &out)
  {
    for(auto it=this->begin();it!=this->end();it++) out<<*it<<endl;
    
    return out;
  }
};
map<int,corr_t<double>> corr_list[1024];

//endianess
void doubles_to_doubles_changing_endianess(double *data,int ndoubles)
{
  for(int idouble=0;idouble<ndoubles;idouble++)
    {
      char *ptr=(char*)(data+idouble);
      
      swap(ptr[7],ptr[0]);
      swap(ptr[6],ptr[1]);
      swap(ptr[5],ptr[2]);
      swap(ptr[4],ptr[3]);
    }
}

//jacknife
const int njacks=16;
class jack : public tr1::array<double,njacks+1>
{
public:
  void clusterize(int nconfs)
  {
    double fact=(double)njacks/(njacks-1)/nconfs;
    for(int iclust=0;iclust<njacks;iclust++) (*this)[njacks]+=(*this)[iclust];
    for(int iclust=0;iclust<njacks;iclust++) (*this)[iclust]=((*this)[njacks]-(*this)[iclust])*fact;
    (*this)[njacks]/=nconfs;
  }
  void change_endianess() {doubles_to_doubles_changing_endianess(&((*this)[0]),njacks+1);}
};

template <class T> T operator+(T a,T b)
{T c;transform(a.begin(),a.end(),b.begin(),c.begin(),plus<double>());return c;}
template <class T> T operator-(T a,T b)
{T c;transform(a.begin(),a.end(),b.begin(),c.begin(),minus<double>());return c;}
template <class T> T operator*(T a,T b)
{T c;transform(a.begin(),a.end(),b.begin(),c.begin(),multiplies<double>());return c;}
template <class T> T operator/(T a,T b)
{T c;transform(a.begin(),a.end(),b.begin(),c.begin(),divides<double>());return c;}

template <class T> T operator+(T a,double b)
{T c;transform(a.begin(),a.end(),c.begin(),bind2nd(plus<double>(),b));return c;}
template <class T> T operator-(T a,double b)
{T c;transform(a.begin(),a.end(),c.begin(),bind2nd(minus<double>(),b));return c;}
template <class T> T operator*(T a,double b)
{T c;transform(a.begin(),a.end(),c.begin(),bind2nd(multiplies<double>(),b));return c;}
template <class T> T operator/(T a,double b)
{T c;transform(a.begin(),a.end(),c.begin(),bind2nd(divides<double>(),b));return c;}

template <class T> T operator+=(T &a,T b){return a=a+b;}
template <class T> T operator-=(T &a,T b){return a=a-b;}
template <class T> T operator*=(T &a,T b){return a=a*b;}
template <class T> T operator/=(T &a,T b){return a=a/b;}

//return the id composing spatial der
inline int get_spat_id(int sm,int mu1,int mu2) {return mu2+4*(mu1+4*sm);}

//find the correct gamma in the list
int ga_unmap(const char *tag)
{
  int iga=0;
  while(strcmp(tag,lga[iga])) iga++;
  if(iga==nga)
    {
      fprintf(stderr,"op %s not present in list\n",tag);
      exit(1);
    }
  return iga;
}

//scan a file
void scan(const char *path)
{
  ifstream fin(path);
  
  string hash;
  while(fin>>hash)
    {
      int iga_so=0,iga_si=0;
      int ispat_so=0,ispat_si=0;
      corr_t<double> temp;
      if(hash=="#")
	{
	  string name;
	  fin>>name;
	  //read all the parameters
	  int ism_so,ism_si,mu1_so,mu2_so,mu1_si,mu2_si;
	  if(name=="sm_so")
	    fin>>ism_so>>hash>>mu1_so>>hash>>mu2_so>>hash>>ism_si>>hash>>mu1_si>>hash>>mu2_si>>hash>>name;
	  iga_si=ga_unmap(name.substr(0,2).c_str());
	  iga_so=ga_unmap(name.substr(2,2).c_str());
	  ispat_si=get_spat_id(ism_si,mu1_si,mu2_si);
	  ispat_so=get_spat_id(ism_so,mu1_so,mu2_so);
	  for(int t=0;t<NT;t++) fin>>temp[t];
	}
      int ispat=ispat_si+32*ispat_so;
      int iga=iga_si+nga*iga_so;
      corr_list[ispat][iga]=temp;
    }
  
  fin.close();
}

//search and return
corr_t<double> &corr(int ism_so,op_contr_t op_so,int ism_si,op_contr_t op_si)
{
  int ispat_so=get_spat_id(ism_so,op_so.mu1,op_so.mu2);
  int ispat_si=get_spat_id(ism_si,op_si.mu1,op_si.mu2);
  int ispat=ispat_si+32*ispat_so;
  int iga=op_si.iga+nga*op_so.iga;
  
  if(corr_list[ispat].find(iga)==corr_list[ispat].end())
    {
      cerr<<"Gamma: "<<iga<<" (si: "<<op_si.iga<<", so: "<<op_so.iga<<") not found in "<<ispat
	  <<" [so: "<<ispat_so<<" (sm: "<<ism_so<<", mu1: "<<op_so.mu1<<", mu2: "<<op_so.mu2<<"), "
	  <<" [si: "<<ispat_si<<" (sm: "<<ism_si<<", mu1: "<<op_si.mu1<<", mu2: "<<op_si.mu2<<")"<<endl;
      exit(1);
    }
  
  return corr_list[ispat][iga];
}

//list
void init_list_op(vector<op_t> &list_op)
{
  op_t P5,A0,CI,VI,TI;
  P5.push_back(op_contr_t(0,0, 5, +1));
  A0.push_back(op_contr_t(0,0, 9, +1));
  for(int i1=0;i1<3;i1++)
    {
      int i2=(i1+1)%3,i3=(i1+2)%3;
      CI.push_back(op_contr_t(   0,1+i1, 16+i1, +1.0));
      VI.push_back(op_contr_t(1+i2,1+i3,  1+i1, +1.0));
      TI.push_back(op_contr_t(1+i2,1+i3, 10+i1, +1.0));
      VI.push_back(op_contr_t(1+i3,1+i2,  1+i1, -1.0));
      TI.push_back(op_contr_t(1+i3,1+i2, 10+i1, -1.0));
    }
  list_op.push_back(P5);
  list_op.push_back(A0);
  list_op.push_back(CI);
  list_op.push_back(VI);
  list_op.push_back(TI);
}

//combine operators
corr_t<double> combine(int ism_so,op_t &op_so,int ism_si,op_t &op_si)
{
  corr_t<double> total(0);
  for(uint icontr_so=0;icontr_so<op_so.size();icontr_so++)
    for(uint icontr_si=0;icontr_si<op_si.size();icontr_si++)
      total+=corr(ism_so,op_so[icontr_so],ism_si,op_si[icontr_si])*op_so[icontr_so].weight*op_si[icontr_si].weight;
  
  return total;
}

//add the various contributions
void put_all(corr_t<jack> *all,int ism_so,int ism_si,int ijack,vector<op_t> &list_op)
{
  int nop=list_op.size();
  for(int iop_so=0;iop_so<nop;iop_so++)
    for(int iop_si=0;iop_si<nop;iop_si++)
      {
	auto contr=combine(ism_so,list_op[iop_so],ism_si,list_op[iop_si]);
	for(int t=0;t<NT;t++) all[iop_si+nop*(ism_si+2*(iop_so+nop*ism_so))][t][ijack]+=contr[t];
      }
}

int main()
{
  vector<op_t> list_op;
  init_list_op(list_op);

  corr_t<jack> *data=new corr_t<jack>[2*2*5*5];
  const int nconfs=240,clust_size=nconfs/njacks;
  //#pragma omp parallel for
  for(int iclust=0;iclust<njacks;iclust++)
    for(int iconf=100+clust_size*iclust;iconf<100+clust_size*(iclust+1);iconf++)
      {
	//read
	char path[100];
	sprintf(path,"out/%3d_0/2pts_corr",iconf);
	scan(path);
	
	if(0)
	for(int i=0;i<6;i++)
	  for(int j=0;j<6;j++)
	    {
	      corr_t<double> t=corr(0,list_op[3][i],0,list_op[3][j]);
	      cout.precision(16);
	      cout<<i<<" "<<list_op[3][i].mu1<<" "<<list_op[3][i].mu2<<" "<<list_op[3][i].iga<<" "<<list_op[3][i].weight
		  <<"  "<<j
		  <<" "<<list_op[3][j].mu1<<" "<<list_op[3][j].mu2<<" "<<list_op[3][j].iga<<" "<<list_op[3][j].weight
		  <<"  "<<t[1]<<endl;
	    }
	
	//put all
	for(int ism_so=0;ism_so<2;ism_so++)
	  for(int ism_si=0;ism_si<2;ism_si++)
	    put_all(data,ism_so,ism_si,iclust,list_op);
	
	cout<<iconf<<" "<<iclust<<endl;
      }
  
  //clusterize
  for(int i=0;i<100;i++)
    for(int t=0;t<NT;t++)
      {
	data[i][t].clusterize(nconfs);
	data[i][t].change_endianess();
      }
  ofstream fout("data");
  fout.write((char*)data,2*2*5*5*NT*17*sizeof(double));
  fout.close();
  
  return 0;
}

