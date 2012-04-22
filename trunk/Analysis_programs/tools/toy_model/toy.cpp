#include "include.h"

int nel=10;
int nestr=25000;
int njack=sqrt(nestr);

ran_gen extr(0);

jack temp_fill_gauss(double med,double sigma)
{
  jack a(njack);
  a=0;
  int clust_size=nestr/njack;
  for(int i=0;i<njack;i++) for(int ies=0;ies<clust_size;ies++) a.data[i]+=extr.get_gauss(med,sigma*sqrt(nestr-1));
  for(int i=0;i<njack;i++) a.data[njack]+=a.data[i];
  for(int i=0;i<njack;i++) a.data[i]=(a.data[njack]-a.data[i])/(nestr-clust_size);
  a.data[njack]/=nestr;
  
  return a;
}

int main()
{
  jvec f(nel,njack);
  
  int nfill=20000;
  int nm1=0,nm2=0,nq1=0,nq2=0;
  double em1=0,em2=0,eq1=0,eq2=0;
  for(int ifill=0;ifill<nfill;ifill++)
    {
      double reg=0.5;
      for(int iel=0;iel<nel;iel++)
	{
	  f[iel]=temp_fill_gauss(0,0.1);
	  if(iel>0) f[iel]=f[iel]*(1-reg)+f[iel-1]*reg;
	}
      for(int iel=0;iel<nel;iel++) f[iel]+=iel+1.3;
      
      //cout<<f<<endl;
      
      jack m(njack),q(njack);
      linear_fit(m,q,f,0,nel);
      em1+=m.err();
      eq1+=q.err();
      if(m.err()>fabs((m.med()-1))) nm1++;
      if(q.err()>fabs((q.med()-1.3))) nq1++;
      
      jvec e(nel-1,njack);
      for(int iel=0;iel<nel-1;iel++) e[iel]=f[iel+1]-f[iel];
      m=constant_fit(e,0,nel-1);
      em2+=m.err();
      if(m.err()>fabs((m.med()-1))) nm2++;
      for(int iel=0;iel<nel;iel++) f[iel]-=m*iel;
      q=constant_fit(f,0,nel);
      eq2+=q.err();
      if(q.err()>fabs((q.med()-1.3))) nq2++;
    }
  
  double pm1=(double)nm1/nfill;
  double pq1=(double)nq1/nfill;
  double pm2=(double)nm2/nfill;
  double pq2=(double)nq2/nfill;
  double epm1=sqrt(pm1*(1-pm1)/nfill);
  double epq1=sqrt(pq1*(1-pq1)/nfill);
  double epm2=sqrt(pm2*(1-pm2)/nfill);
  double epq2=sqrt(pq2*(1-pq2)/nfill);
  
  cout<<pm1<<" "<<epm1<<"\t"<<pq1<<" "<<epq1<<endl;
  cout<<pm2<<" "<<epm2<<"\t"<<pq2<<" "<<epq2<<endl;
  
  cout<<em1/nfill<<" "<<eq1/nfill<<endl;
  cout<<em2/nfill<<" "<<eq2/nfill<<endl;
  
  return 0;
}
