#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <cmath>
#include <cstdarg>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <omp.h>

#include "random.hpp"

using namespace std;

const double log_eps=-numeric_limits<double>::infinity();
inline double sqr(double x){return x*x;}
rnd_gen_t glb_rnd_gen;

//crash promptin error message
#define crash(...) internal_crash(__LINE__,__FILE__,__VA_ARGS__)
void internal_crash(int line,const char *file,const char *temp,...)
{
  char buffer[1024];
  va_list args;
  
  va_start(args,temp);
  vsprintf(buffer,temp,args);
  va_end(args);
  
  cerr<<"ERROR at line "<<line<<" of file "<<file<<": "<<buffer<<endl;
  exit(1);
}

//compute the logarithm of the sum of the exponentials
inline void log_sum(double &sum,double add)
//{sum=add+log1p(exp(sum-add));}
{
  double diff=sum-add;
  if(diff>0) sum+=log1p(exp(-diff));
  else sum=add+log1p(exp(diff));
}

//compute the logarithm of an observable (it may be Z itself) from the info on simulated ensembles
void calc_log_z_powoss(double *out_log_z,int noss,double ***out_log_oss,int nout_points,double *out_points,double *in_log_z,int nin_setup,double *in_setup,double **azione,double ***in_log_oss,double *powo,int npowo,int *nconf,double *log_nconf)
{
  for(int iout=0;iout<nout_points;iout++) out_log_z[iout]=log_eps;
  for(int iout=0;iout<nout_points;iout++)
    for(int ioss=0;ioss<noss;ioss++)
      for(int ipowo=0;ipowo<npowo;ipowo++)
	out_log_oss[iout][ioss][ipowo]=log_eps;
  
#pragma omp parallel for
  for(int iout=0;iout<nout_points;iout++)
    for(int is=0;is<nin_setup;is++)
      for(int ic=0;ic<nconf[is];ic++)
	{
	  double log_den=log_eps;
	  for(int isk=0;isk<nin_setup;isk++)
	    {
	      double expo=(out_points[iout]-in_setup[isk])*azione[is][ic];
	      log_sum(log_den,log_nconf[isk]-in_log_z[isk]+expo);
	    }
	  log_sum(out_log_z[iout],-log_den);
	  
	  for(int ioss=0;ioss<noss;ioss++)
	    for(int ipowo=0;ipowo<npowo;ipowo++)
	      log_sum(out_log_oss[iout][ioss][ipowo],-(log_den-powo[ipowo]*in_log_oss[is][ic][ioss]));
	}
}

//Wrapper for the calculation of z
void calc_log_z(double *out_log_z,int nout_points,double *out_points,double *in_log_z,int nin_setup,double *in_setup,double **azione,int *nconf,double *log_nconf)
{
  int noss=0;
  double ***in_log_oss=NULL;
  double ***out_log_oss=NULL;
  double *powo=NULL;
  int npowo=0;
  
  calc_log_z_powoss(out_log_z,noss,out_log_oss,nout_points,out_points,in_log_z,nin_setup,in_setup,azione,in_log_oss,powo,npowo,nconf,log_nconf);
}

//Calculate new z starting from prior
void multi_z_rew(int nsetup,double *setup,double **azione,int *nconf,double *log_nconf,double *log_z,const char *tag,int ord_rec)
{
  double *new_log_z=new double[nsetup];
  double res;
  double shift;
  
  do
    {
      //Calculate new log_z
      calc_log_z(new_log_z,nsetup,setup,log_z,nsetup,setup,azione,nconf,log_nconf);
      
      //Shift and calculate the residual
      res=0;
      shift=new_log_z[nsetup/2];
      for(int is=0;is<nsetup;is++)
        {
          //Calculate stopping condition
          res+=sqr(expm1(new_log_z[is]-log_z[is]));
          log_z[is]=new_log_z[is]-shift; //Save new values of log_z
        }
      for(int iord=0;iord<ord_rec;iord++) cout<<" ";
      cout<<tag<<res<<" "<<nsetup<<endl;
    }
  while(res>1.e-7);
  
  delete [] new_log_z;
}

// Implementazione MultiLivello ricorsiva
void rec_multi_z_rew(int nsetup,double *setup,double **azione,int *nconf,double *log_nconf,double *log_z,bool guess,const char *tag="",int ord_rec=0)
{
  //Spacca in due blocchi solo se ci sono più di 3 punti
  if(nsetup>3 && guess==false) 
    { 
      int nsetup_half=nsetup/2;
      double *log_z1=new double[nsetup_half+1];
      double *log_z2=new double[nsetup-nsetup_half];
      
      //Crea due sotto blocchi di z
      for(int is=0;is<nsetup;is++)
        {
          if(is<=nsetup_half) log_z1[is]=log_z[is];
          if(is>=nsetup_half) log_z2[is-nsetup_half]=log_z[is];
        }
      
      //Ripesa i due blocchi separatamente
      int nsetup1=nsetup_half+1,nsetup2=nsetup-nsetup_half;
      double *setup1=setup,*setup2=setup+nsetup_half;
      double **azione1=azione,**azione2=azione+nsetup_half;
      int *nconf1=nconf,*nconf2=nconf+nsetup_half;
      double *log_nconf1=log_nconf,*log_nconf2=log_nconf+nsetup_half;
      rec_multi_z_rew(nsetup1,setup1,azione1,nconf1,log_nconf1,log_z1,guess,">a>",ord_rec+1);
      rec_multi_z_rew(nsetup2,setup2,azione2,nconf2,log_nconf2,log_z2,guess,"<b<",ord_rec+1);
      
      //Unisce i due blocchi
      for(int is=0;is<nsetup;is++)
        if(is<=nsetup_half) log_z[is]=log_z1[is];
        else 
          log_z[is]=log_z[is-1]+log_z2[is-nsetup_half]-log_z2[is-nsetup_half-1];
      
      delete [] log_z1;
      delete [] log_z2;
    }
  
  //Ripesamento finale
  multi_z_rew(nsetup,setup,azione,nconf,log_nconf,log_z,tag,ord_rec);
}

double boot_err(double *v,int n)
{
  double sum=0;
  double sum2=0;
  
  for(int i=1;i<n;i++)
    {
      sum+=v[i];
      sum2+=v[i]*v[i];
    }
  sum/=n-1;
  sum2/=n-1;
  
  return sqrt(sum2-sum*sum);
}

void crea_osservabili(double *out_oss,double *x)
{
  double x1=x[0],x2=x[1],x3=x[2],x4=x[3];
  out_oss[0]=x1;
  out_oss[1]=x2-x1*x1;
  out_oss[2]=1-x4/(3*x2*x2);
  //out_oss[3]=(x4-4*x3*x1+6*x2*x1*x1-3*x1*x1*x1*x1)/sqr(x2-x1*x1);
  out_oss[3]=sqr(x2)/x4;
}

int contaelementi(string perc)
{
  int i=0;
  ifstream file(perc.c_str());
  if(!file.good()) crash("error opening file: %s",perc.c_str());
  string dum;
  
  while(file>>dum) i++;
  
  return i;
}

int contarighe(string perc)
{
  int i=0;
  ifstream file(perc.c_str());
  if(!file.good()) crash("error opening file: %s",perc.c_str());
  char dum[10000];
  
  while(file.getline(dum,10000)) i++;
  
  return i;
}

void compute_partial_weights(int nout_points,double *out_points,double *log_z,int nraw_setup,double *raw_par,double **raw_azione,int *nconf,double *log_nconf)
{
  double part_z[nout_points][nraw_setup];
  double tot_z[nout_points];
#pragma omp parallel for
  for(int iout=0;iout<nout_points;iout++)
    {
      tot_z[iout]=log_eps;
      for(int is=0;is<nraw_setup;is++)
	{
	  part_z[iout][is]=log_eps;
	  for(int ic=0;ic<nconf[is];ic++)
	    {
	      double log_den=log_eps;
	      for(int isk=0;isk<nraw_setup;isk++)
		{
		  double expo=(out_points[iout]-raw_par[isk])*raw_azione[is][ic];
		  log_sum(log_den,log_nconf[isk]-log_z[isk]+expo);
		}
	      log_sum(part_z[iout][is],-log_den);
	    }
	  log_sum(tot_z[iout],part_z[iout][is]);
	}
    }
  
  //write down
  ofstream part_z_out("part_weights");
  for(int iout=0;iout<nout_points;iout++)
    {
      part_z_out<<out_points[iout]<<" ";
      for(int is=0;is<nraw_setup;is++) part_z_out<<exp(part_z[iout][is]-tot_z[iout])<<" ";
      part_z_out<<endl;
    }
  part_z_out.close();
}

int main(int narg,char **arg)
{
  cout.precision(16);
  cout<<scientific;
  
  cout<<"Using "<<omp_get_max_threads()<<" threads"<<endl;
  glb_rnd_gen.init(10032);
  
  //Conta le righe
  string conf_path="configfile";
  int nraw_setup=contarighe(conf_path.c_str());
  ifstream conf_file(conf_path.c_str());
  if(!conf_file.good()||nraw_setup==0) crash("check %s",conf_path.c_str());
  
  //Definisce gli array per i dati
  int nall_conf[nraw_setup];          //Numero di tutte le conf. file
  int nconf[nraw_setup];              //Numero di conf effettivo
  double log_nconf[nraw_setup];       //Logaritmo di nconf effettivo
  double raw_par[nraw_setup];         //Beta,massa,etc
  int aut[nraw_setup];                //Autocorrelazione
  string *input_path=new string[nraw_setup];   //Percorsi input
  
  //Carica file configurazione, e ne conta le righe
  string aut_path("autoc");
  ifstream aut_file(aut_path.c_str());
  for(int is=0;is<nraw_setup;is++)
    {
      double temp;
      if(!(conf_file>>raw_par[is])) crash("loading raw_par[%d]",is);
      if(!(conf_file>>input_path[is])) crash("reading path %d",is);
      nall_conf[is]=contarighe(input_path[is]);
      if(nall_conf[is]==0) crash("file %s has 0 lines",input_path[is].c_str());
      if(!(aut_file>>temp)) crash("check file: %s",aut_path.c_str());
      temp=2*temp+1;
      aut[is]=int(temp+0.5);
    }
  
  //Conta le osservabili e definisce i vettori con i dati raw
  int ndati=contaelementi(input_path[0].c_str())/nall_conf[0];
  int noss=ndati-1;
  double *raw_azione[nraw_setup];
  double **log_raw_oss[nraw_setup];
  double **raw_oss[nraw_setup];
  cout<<"N. observables: "<<noss<<endl;
  cout<<"N. raw setups: "<<nraw_setup<<endl;
  
  for(int is=0;is<nraw_setup;is++)
    {
      raw_oss[is]=new double*[nall_conf[is]];
      log_raw_oss[is]=new double*[nall_conf[is]];
      raw_azione[is]=new double[nall_conf[is]];
      for(int ic=0;ic<nall_conf[is];ic++)
	{
          log_raw_oss[is][ic]=new double[noss];
          raw_oss[is][ic]=new double[noss];
        }
    }
  
  //Carica i dati raw
  cout<<"Loading data"<<endl;
#pragma omp parallel for
  for(int is=0;is<nraw_setup;is++)
    {
      ifstream input_file(input_path[is].c_str());
      if(!input_file) crash("unable to open file: %s",input_path[is].c_str());
      
      //Carica nel posto giusto o nel temporaneo
      for(int ic=0;ic<nall_conf[is];ic++)
        {
          input_file>>raw_azione[is][ic];
          for(int io=0;io<noss;io++)
            {
              input_file>>raw_oss[is][ic][io];
              log_raw_oss[is][ic][io]=raw_oss[is][ic][io];
            }
	  
          if(!input_file.good())
            {
              cerr<<"Erorre in file: "<<input_path[is]<<endl;
              exit(1);
            }
        }
      input_file.close();
    }
  
  //Sposta in su del minimo l'osservabile
  double minoss[noss];
  for(int io=0;io<noss;io++)
    {
      minoss[io]=log_raw_oss[0][0][io];
      for(int is=0;is<nraw_setup;is++)
        for(int ic=0;ic<nall_conf[is];ic++)
          minoss[io]=min(minoss[io],log_raw_oss[is][ic][io]);
    }
  for(int io=0;io<noss;io++)
    {
      if(minoss[io]<0) cout<<"Osservabile "<<io<<" spostata"<<endl;
      for(int is=0;is<nraw_setup;is++)
        for(int ic=0;ic<nall_conf[is];ic++)
          {
            if(minoss[io]<0) log_raw_oss[is][ic][io]-=2*minoss[io];
            log_raw_oss[is][ic][io]=log(log_raw_oss[is][ic][io]);
          }
    }
  
  ///////////////////// inizio reweighting ///////////////////////
  
  //Parsa i parmetri di input
  const int nboot=100;
  const int npowo=4;
  double powo[npowo]={1,2,3,4};
  
  int *iswap[nraw_setup];
  for(int is=0;is<nraw_setup;is++) iswap[is]=new int[nall_conf[is]/aut[is]];
  
  double log_z[nraw_setup];
  
  string z_path="z_save";
  fstream *z_file=new fstream;
  
  bool z_good=true;
  z_file->open(z_path.c_str(),std::ios::in);
  
  int nout_points=contarighe("base_out");
  if(nout_points==0)
    {
      cerr<<"File base_out assente!"<<endl;
      exit(1);
    }
  
  double *out_points=new double[nout_points];
  
  double *out_log_z=new double[nout_points];
  double ***out_log_oss=new double**[nout_points];
  for(int ip=0;ip<nout_points;ip++)
    {
      out_log_oss[ip]=new double*[noss];
      for(int io=0;io<noss;io++) out_log_oss[ip][io]=new double[npowo];
    }
  
  //Vettori per osservabili ripesate
  double out_rew_oss[nout_points][noss][npowo][nboot];
  //Vettori per osservabili non ripesate
  double out_oss[nraw_setup][noss][npowo][nboot];
  
  ifstream file_out_points("base_out");
  for(int ip=0;ip<nout_points;ip++) 
    if(!(file_out_points>>out_points[ip]))
      {
	cerr<<"Errore in file base_out"<<endl;
	exit(1);
      }
  
  //Parte principale
  for(int iboot=0;iboot<nboot;iboot++)
    {
      cout<<"Boot "<<iboot+1<<"/"<<nboot<<endl;
      
      for(int is=0;is<nraw_setup;is++)
        {
          int tau;
          if(iboot==0) tau=1;
          else tau=aut[is];
	  
          nconf[is]=nall_conf[is]/tau;
          log_nconf[is]=log(nconf[is]); //Log N
          
          //piglia 1 per blocco
          if(iboot!=0) 
            for(int ic=0;ic<nconf[is];ic++)
              {
                int iw=ic*tau+glb_rnd_gen.get_unif(0,tau);
                iswap[is][ic]=iw;
		
                for(int io=0;io<noss;io++)
                  swap(log_raw_oss[is][ic][io],log_raw_oss[is][iw][io]);
		swap(raw_azione[is][ic],raw_azione[is][iw]);
              }
        }
      
      //Calcola le osservabili sui soli file originali
      for(int is=0;is<nraw_setup;is++)
        for(int io=0;io<noss;io++)
          for(int ipowo=0;ipowo<npowo;ipowo++)
            {
              out_oss[is][io][ipowo][iboot]=0;
              for(int ic=0;ic<nconf[is];ic++)
                {
                  double temp_oss;
                  if(iboot==0) temp_oss=raw_oss[is][ic][io];
                  else temp_oss=raw_oss[is][iswap[is][ic]][io];
                  out_oss[is][io][ipowo][iboot]+=pow(temp_oss,powo[ipowo]);
                }
              out_oss[is][io][ipowo][iboot]/=nconf[is];
            }
      
      //Copia le osservabili non ripesate nel mucchio
      for(int is=0;is<nraw_setup;is++) 
        for(int io=0;io<noss;io++)
          {
            double temp_oss[npowo],temp_in[npowo];
            for(int ipowo=0;ipowo<npowo;ipowo++)
              temp_in[ipowo]=out_oss[is][io][ipowo][iboot];
            crea_osservabili(temp_oss,temp_in);
            for(int ipowo=0;ipowo<npowo;ipowo++)
                out_oss[is][io][ipowo][iboot]=temp_oss[ipowo];
          }
            
      //Prova a caricare z
      if(z_good)
        {
	  cout<<"Reading Z"<<endl;
          for(int is=0;is<nraw_setup;is++)z_good=z_good and ((*z_file)>>log_z[is]);
          if(!z_good)
            {
              cout<<"Opening "<<z_path<<" in writing mode"<<endl;
              delete z_file;
              z_file=new fstream;
              (*z_file)<<scientific;
              z_file->precision(16);
              z_file->open(z_path.c_str(),std::ios::out | std::ios::app);
              if(!z_file->good()) {cerr<<"Boh"<<endl;exit(3);}
            }
        }
      
      if(!z_good)
        {
          //Reweighting
	  cout<<"Computing z"<<endl;
          if(iboot==0) for(int is=0;is<nraw_setup;is++) log_z[is]=0;
          rec_multi_z_rew(nraw_setup,raw_par,raw_azione,nconf,log_nconf,log_z,iboot);
          
          //Output
          for(int is=0;is<nraw_setup;is++) (*z_file)<<log_z[is]<<endl;
        }
      
      if(iboot==0)
	{
	  cout<<" Computing partial weights"<<endl;
	  compute_partial_weights(nout_points,out_points,log_z,nraw_setup,raw_par,raw_azione,nconf,log_nconf);
	}
      
      cout<<" Reweighting observables"<<endl;
      calc_log_z_powoss(out_log_z,noss,out_log_oss,nout_points,out_points,log_z,nraw_setup,raw_par,raw_azione,log_raw_oss,powo,npowo,nconf,log_nconf);
      
      //Copia le osservabili ripesate nel mucchio
      cout<<" Create observables"<<endl;
      for(int ip=0;ip<nout_points;ip++)
        for(int io=0;io<noss;io++)
          {
            double temp_oss[npowo],temp_in[npowo];
	    
            for(int ipowo=0;ipowo<npowo;ipowo++)
              temp_in[ipowo]=exp(out_log_oss[ip][io][ipowo]-out_log_z[ip]);
	    
            crea_osservabili(temp_oss,temp_in);
            for(int ipowo=0;ipowo<npowo;ipowo++)
              {
                out_rew_oss[ip][io][ipowo][iboot]=temp_oss[ipowo];
                if(ipowo==0 && minoss[io]<0)
                  out_rew_oss[ip][io][ipowo][iboot]+=2*minoss[io];
              }
          }
      
      //rimette a posto
      if(iboot!=0)
        for(int is=0;is<nraw_setup;is++)
          for(int ic=0;ic<nconf[is];ic++)
            {
              int iw=iswap[is][ic];
              
              for(int io=noss-1;io>=0;io--)
                swap(log_raw_oss[is][ic][io],log_raw_oss[is][iw][io]);
               swap(raw_azione[is][ic],raw_azione[is][iw]);
            }
    }
  
  string fout_path[npowo]={"osservabili","suscettivita","binder","binder_malvagio"};  
  
  //Stampa le osservabili non ripesate
  for(int ipowo=0;ipowo<npowo;ipowo++)
    {
      ofstream fout(fout_path[ipowo].c_str());
      fout.precision(16);
      fout<<scientific;
      for(int is=0;is<nraw_setup;is++) 
	{
	  fout<<raw_par[is]<<"\t";
	  for(int io=0;io<noss;io++)
	    fout<<out_oss[is][io][ipowo][0]<<" "<<boot_err(out_oss[is][io][ipowo],nboot)<<"\t";
	  fout<<endl;
	}
      fout.close();
    }
  
  //Stampa le osservabili ripesate
  for(int ipowo=0;ipowo<npowo;ipowo++)
    {
      ofstream fout((fout_path[ipowo]+"_rew").c_str());
      fout.precision(16);
      fout<<scientific;
      for(int ip=0;ip<nout_points;ip++) 
	{
	  fout<<out_points[ip]<<"\t";
	  for(int io=0;io<noss;io++)
	    fout<<out_rew_oss[ip][io][ipowo][0]<<" "<<boot_err(out_rew_oss[ip][io][ipowo],nboot)<<"\t";
	  fout<<endl;
	}
      fout.close();
    }
  
  //Salva i dati grezzi ripesati
  FILE *f=fopen("data.raw","w");
  fwrite(&nraw_setup,sizeof(int),1,f);
  fwrite(&nout_points,sizeof(int),1,f);
  fwrite(&noss,sizeof(int),1,f);
  fwrite(&npowo,sizeof(int),1,f);
  fwrite(&nboot,sizeof(int),1,f);
  for(int is=0;is<nraw_setup;is++) fwrite(&(raw_par[is]),sizeof(double),1,f);
  for(int is=0;is<nraw_setup;is++)
    for(int io=0;io<noss;io++)
      for(int ipowo=0;ipowo<npowo;ipowo++)
	{
	  fwrite(out_oss[is][io][ipowo]+1,sizeof(double),nboot-1,f);
	  fwrite(out_oss[is][io][ipowo],sizeof(double),1,f);
	}
  for(int ip=0;ip<nout_points;ip++) fwrite(&(out_points[ip]),sizeof(double),1,f);
  for(int ip=0;ip<nout_points;ip++) 
    for(int io=0;io<noss;io++)
      for(int ipowo=0;ipowo<npowo;ipowo++)
	{       
	  fwrite(out_rew_oss[ip][io][ipowo]+1,sizeof(double),nboot-1,f);
	  fwrite(out_rew_oss[ip][io][ipowo],sizeof(double),1,f);
	}
  fclose(f);

  return 0;
}
