#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TSystem.h>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

const double chrono_topo_coeff=0.3;
const double chrono_topo_width=0.3;
double Q_min=+1000000;
double Q_max=-1000000;

using namespace std;

vector<double> chrono_topo_past_values;

TApplication myapp("app",NULL,NULL);
TCanvas tela;
TGraph graph;

//compute the topodynamical potential using past history
double compute_theta_pot(double Q)
{
  double topotential=0;
  for(vector<double>::iterator it=chrono_topo_past_values.begin();it!=chrono_topo_past_values.end();it++)
    {
      double q=*it;
      double diff=Q-q,f=diff/chrono_topo_width;
      if(fabs(f)<5)
        {
          double cont=exp(-f*f/2);
          topotential+=cont;
        }
    }
  topotential*=chrono_topo_coeff;
    
  return topotential;
}

//draw the chronological topological potential
void draw_chrono_topo_potential()
{
  double Q_diff=Q_max-Q_min;
  int n=ceil(Q_diff/chrono_topo_width*20);
  if(n==0) n=1;
  double dQ=Q_diff/n;
  
  //compute
  ofstream out("/tmp/out.xmg");
  for(int i=0;i<=n;i++)
    {
      double x=Q_min+i*dQ;
      double y=compute_theta_pot(x);
      graph.SetPoint(i,x,y);
      out<<x<<" "<<y<<endl;
    }
  tela.Modified();
  tela.Update();
}

int main()
{
  tela.SetWindowPosition(100,100);
  tela.SetWindowSize(1200,1200);
  graph.Draw("APL");

  gSystem->ProcessEvents();

  std::ifstream fin("/tmp/twith");
  
  double x,y;
  int i=0,n=10;
  while(fin>>x>>y)
    {
      chrono_topo_past_values.push_back(+y);
      chrono_topo_past_values.push_back(-y);
      
      Q_max=max(Q_max,y);
      Q_min=min(Q_min,y);

      Q_max=max(Q_max,-y);
      Q_min=min(Q_min,-y);

      //if((i%n)==(n-1)) draw_chrono_topo_potential();
      i++;
    }
  
  draw_chrono_topo_potential();
  myapp.Run();
  
  return 0;
}
