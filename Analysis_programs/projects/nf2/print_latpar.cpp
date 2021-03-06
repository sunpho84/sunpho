#include "common.cpp"

int main()
{
  init_latpars();
  
  cout<<"ml_phys: "<<ml_phys<<endl;
  cout<<"ms_phys: "<<ms_phys<<endl;
  cout<<"mc_phys: "<<mc_phys<<endl;
  cout<<endl;
  cout<<"Zp380: "<<Zp[0]<<endl;
  cout<<"a380: "<<lat[0]<<endl;
  cout<<"amc_phys380: "<<amc_phys[0]<<endl;
  cout<<"1/a380: "<<1/lat[0]<<endl<<endl;
  
  cout<<"Zp390: "<<Zp[1]<<endl;
  cout<<"a390: "<<lat[1]<<endl;
  cout<<"1/a390: "<<1/lat[1]<<endl<<endl;
  
  cout<<"Zp405: "<<Zp[3]<<endl;
  cout<<"a405: "<<lat[2]<<endl;
  cout<<"1/a405: "<<1/lat[2]<<endl<<endl;
  
  cout<<"Zp420: "<<Zp[3]<<endl;
  cout<<"a420: "<<lat[3]<<endl;
  cout<<"amc_phys420: "<<amc_phys[3]<<endl;
  cout<<"1/a420: "<<1/lat[3]<<endl<<endl;

  cout<<"mc: "<<mc_phys<<endl<<endl;
  cout<<"r0: "<<smart_print(r0)<<endl<<endl;
  
  double ghat=0.45;
  //double Amed=-1.200,Aerr=0.5;
  double Amed=-2.000,Aerr=0.9;
  boot A(nboot,njack);
  A.fill_gauss(Amed,Aerr,28732);
  
  cout<<A<<endl;
  cout<<"db0: "<<db0<<endl;
  //cout<<"r0: "<<r0<<endl;
  cout<<"f0:"<<f0<<endl;
  
  boot arg=4*M_PI*f0;
  double c=3.0/4*(1+3*sqr(ghat));
  boot den=arg*arg;

  boot corr1=A/db0;
  boot corr2=c/den*log(den);
  
  cout<<corr1<<endl;
  cout<<corr2<<endl;
  cout<<corr1+corr2<<endl;
  
  ofstream t("/tmp/data_cecilia");
  t<<"iboot\t\t2B0\t\tf0"<<endl;
  for(int iboot=0;iboot<nboot;iboot++)
    t<<iboot+1<<"\t\t"<<db0[iboot]<<"\t\t"<<f0[iboot]<<endl;
  t<<"med\t\t"<<db0.med()<<"\t\t"<<f0.med()<<endl;
  t<<"err\t\t"<<db0.err()<<"\t"<<f0.err()<<endl;
  
  cout<<db0<<endl;
  cout<<f0<<endl;
  cout<<"ccond: ["<<smart_print(pow(db0*f0*f0/4,1.0/3)*1000)<<" MeV]^3"<<endl;
  
  return 0;
}
