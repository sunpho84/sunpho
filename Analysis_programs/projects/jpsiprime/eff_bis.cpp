jvec effective_mass(double *x,jvec y,int TH,int par=1)
{
  int njack=y.njack;
  int nx=y.nel;
  jvec b(nx-1,njack);
  
  for(int ix=0;ix<nx-1;ix++)
    {
      jack temp=-log(y[ix+1]/y[ix]);
      double miniz=temp.med();
      double einiz=temp.err();
      if(einiz==0) einiz=1.e-10;

      for(int ijack=0;ijack<=njack;ijack++)
        {
          double m=miniz;
          double e=einiz;

          double targ=y[ix+1][ijack]/y[ix][ijack];

          double yl;
          double yr;

          int q;
          do
            {
	      double tp1=x[ix+1],t=x[ix];
	      if(par==1)
		{
		  yl=cosh((m-e)*(TH-tp1))/cosh((m-e)*(TH-t))-targ;
		  yr=cosh((m+e)*(TH-tp1))/cosh((m+e)*(TH-t))-targ;
		}
	      else
		{
		  yl=sinh((m-e)*(TH-tp1))/sinh((m-e)*(TH-t))-targ;
		  yr=sinh((m+e)*(TH-tp1))/sinh((m+e)*(TH-t))-targ;
		}

              q=((yl<0 && yr<0) || (yl>=0 && yr>=0));
              //cout<<t<<" "<<ijack<<" "<<yl<<" "<<yr<<" "<<e<<endl;                                                                                                                                                                                                        
              if(q)
                {
                  e*=2;
                  if(m<=e) m+=(e-m);
                }
            }
          while(q);
          double xl=m-e,xr=m+e,ym;
          do
            {
	      double tp1=x[ix+1],t=x[ix];
              m=(xl+xr)/2;
              ym=cosh(m*(TH-(tp1)))/cosh(m*(TH-t))-targ;
              if((yl<0 && ym<0) || (yl>0 && ym>0))
                xl=m;
              else
                xr=m;
            }
          while(fabs(ym)>1.e-14);
          b.data[ix].data[ijack]=m;
        }
    }
  return b;
}
