#include "find_theta_internal.cpp"

double m1=1.49356,m2=1.55528;

int main()
{
  double bp1=find_p1_bf(m1,m2);
  double bp2=-bp1;
  
  cout<<bp1<<" "<<Q2_fun(m1,bp1,m2,bp2)<<endl;
  
  bp1=find_p1_bf(m1,m2,latt_e);
  bp2=-bp1;
  
  cout<<bp1<<" "<<Q2_fun(m1,bp1,m2,bp2,latt_e)<<endl;
  
  return 0;
}
