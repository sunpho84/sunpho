string write_polygon(double *x,VTYPE y)
{
  ostringstream out;

  out<<"@type xy"<<endl;
  
  //lower line
  for(int i=0;i<y.nel;i++)    out<<x[i]<<" "<<y.data[i].med()-y[i].err()<<endl;
  
  //upper line
  for(int i=y.nel-1;i>=0;i--) out<<x[i]<<" "<<y.data[i].med()+y[i].err()<<endl;

  out<<"&"<<endl;
  
  return out.str();
}

string write_ave_line(double *x,VTYPE y)
{
  ostringstream out;

  out<<"@type xy"<<endl;
  
  //line
  for(int i=0;i<y.nel;i++)    out<<x[i]<<" "<<y.data[i].med()<<endl;
  
  out<<"&"<<endl;
  
  return out.str();
}

string write_poly_with_error(VTYPE p,double xmin,double xmax,int npoints=100)
{
  double x[npoints];
  
#ifdef BOOT
  VTYPE y(npoints,p.nboot,p.njack);
#else
  VTYPE y(npoints,p.njack);
#endif
  for(int ip=0;ip<npoints;ip++)
    {  
      x[ip]=xmin+(xmax-xmin)/(npoints-1)*ip;
      y[ip]=p[0];
      
      double R=x[ip];
      for(int ipow=1;ipow<p.nel;ipow++)
	{
	  y[ip]+=p[ipow]*R;
	  R*=x[ip];
	}
    }
  
  ostringstream out;
  out<<write_polygon(x,y);
  out<<write_ave_line(x,y);
  
  return out.str();
}

string write_line_with_error(TYPE q,TYPE m,double xmin,double xmax,int npoints=100)
{
#ifdef BOOT
  VTYPE p(2,m.nboot,m.njack);
#else
  VTYPE p(2,m.njack);
#endif
  p[0]=q;
  p[1]=m;
  
  return write_poly_with_error(p,xmin,xmax,npoints);
}

string write_constant_with_error(TYPE c,double xmin,double xmax)
{
  TYPE q=c;
  TYPE m=c*0;
  
  return write_line_with_error(q,m,xmin,xmax,2);
}
