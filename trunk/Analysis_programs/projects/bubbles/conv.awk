{
    if($1!=$1+0)
    {
	n=0;
	if($1=="#")
	{
	    printf(" ");
	    for(i=0;i<NF/2;i++) printf("%s ",$(i+1));
	}
	printf("\n")
    }
    else
    {
	a[n]=$1;
	b[n]=$2;
	c[n]=$3;
	d[n]=$4;
	n++;
	if(n==48)
	{
	    n=0;
	    for(dt=0;dt<48;dt++)
	    {
		f=g=0;
		t1=0;
		for(t1=0;t1<48;t1++)
		{
		    t2=(t1+dt)%48;
		    f+=a[t1]*c[t2];
		    g+=b[t1]*d[t2];
		}
		printf("%+016.16lg %+016.16lg\n",f/48,g/48);
	    }
	}
    }
}