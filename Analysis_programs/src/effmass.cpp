#pragma once

double effective_mass(double a,double b,int t,int TH,double m=1,double e=1.e-10,int par=1,int dt=1)
{
  double targ=b/a;
  
  double yl;
  double yr;
  
  //increment the range up to reaching opposite sign
  int q;
  do
    {
      if(par==1)
	{
	  yl=cosh((m-e)*(TH-(t+dt)))/cosh((m-e)*(TH-t))-targ;
	  yr=cosh((m+e)*(TH-(t+dt)))/cosh((m+e)*(TH-t))-targ;
	}
      else
	{
	  yl=sinh((m-e)*(TH-(t+dt)))/sinh((m-e)*(TH-t))-targ;
	  yr=sinh((m+e)*(TH-(t+dt)))/sinh((m+e)*(TH-t))-targ;
	}
      q=((yl<0 && yr<0) || (yl>=0 && yr>=0));
      if(q)
	{
	  e*=2;
	  if(m<=e) m+=(e-m);
	}
    }
  while(q);
  
  //bisect
  double xl=m-e,xr=m+e,ym;
  do
    {
      m=(xl+xr)/2;
      if(par==1) ym=cosh(m*(TH-(t+dt)))/cosh(m*(TH-t))-targ;
      else       ym=sinh(m*(TH-(t+dt)))/sinh(m*(TH-t))-targ;
      if((yl<0 && ym<0) || (yl>0 && ym>0))
	xl=m;
      else
	xr=m;
    }
  while(fabs(ym)>1.e-14);
  
  return m;
}
