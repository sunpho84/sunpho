#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "pos.hpp"
#include "tile.hpp"

#include <TApplication.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraph.h>
#include <TPolyLine.h>
#include <TArc.h>
#include <TCrown.h>
#include <TText.h>

TApplication myapp("App",NULL,NULL);
TCanvas tela("Esagoni","Esagoni",0,0,800,800);
TPad pad("pad","",0,0,1,1);
TPolyLine poly[19];
TArc red_arc,white_arc;
TPolyLine red_rect[19],white_rect[19];
TText labels[19];
TCrown red_crown,white_crown;

const double R=sqrt(3);
const double corner_coords[7][2]={{-1,-R},{+1,-R},{+2,0},{+1,+R},{-1,+R},{-2,0},{-1,-R}};

//initialize the graphic
void init_graph()
{
  pad.Range(-4,-4,+4,+4);
  pad.Draw();
  gPad=&pad;
  
  white_arc.SetLineColor(WHI);
  white_arc.SetFillColor(WHI);
  white_arc.Draw("SAME");

  red_arc.SetLineColor(RED);
  red_arc.SetFillColor(RED);
  red_arc.Draw("SAME");
  
  for(int i=0;i<19;i++)
    {
      white_rect[i].SetLineColor(WHI);
      white_rect[i].SetFillColor(WHI);
      
      red_rect[i].SetLineColor(RED);
      red_rect[i].SetFillColor(RED);
    }
  
  white_crown.SetLineColor(WHI);
  white_crown.SetFillColor(WHI);
  white_crown.Draw("SAME");
  
  red_crown.SetLineColor(RED);
  red_crown.SetFillColor(RED);
  red_crown.Draw("SAME");
}

//draw an exagon
void draw_exa(int ele,double x,double y,double l,bool col)
{
  for(int icor=0;icor<=6;icor++) poly[ele].SetPoint(icor,x+l*corner_coords[icor][0],y+l*corner_coords[icor][1]);
  poly[ele].SetFillColor(col?RED:WHI);
  poly[ele].Draw("f");
  poly[ele].Draw();
}

//add an arc
void add_arc(TArc &arc,double x,double y,double r,double amin,double amax)
{arc.DrawArc(x,y,r,amin,amax);}
void add_arc(double x,double y,double l,int icor,int col)
{add_arc(col?red_arc:white_arc,x+l*corner_coords[icor][0],y+l*corner_coords[icor][1],l,60*icor,120+60*icor);}

//draw a crown
void add_crown(TCrown &crown,double x,double y,double l,int irot)
{
  double ang=60*irot;
  double c=cos(M_PI*(ang+30)/180),s=sin(M_PI*(ang+30)/180);
  crown.DrawCrown(x-R*2*l*c,y-R*2*l*s,2*l,2*l*1.5,ang,ang+60);
}

//draw a rectangle
void add_rectangle(TPolyLine &rect,double *c)
{
  for(int i=0;i<4;i++) rect.SetPoint(i,c[2*i+0],c[2*i+1]);
  rect.Draw("f");
}
void add_rectangle(TPolyLine &rect,double x,double y,double l,int itile,int icor)
{
  double c[8];
  for(int ic=0;ic<2;ic++)
    {
      c[2*0+ic]=l*corner_coords[(icor+1)%6][ic];
      c[2*1+ic]=l*corner_coords[(icor+5)%6][ic];//(c[2*0+ic]+l*corner_coords[(icor+2)%6][ic])/2;
      c[2*2+ic]=l*(corner_coords[(icor+4)%6][ic]+corner_coords[(icor+5)%6][ic])/2;
      c[2*3+ic]=l*(corner_coords[(icor+1)%6][ic]+corner_coords[(icor+2)%6][ic])/2;
    }
  for(int i=0;i<4;i++)
    {
      c[2*i+0]+=x;
      c[2*i+1]+=y;
    }

  add_rectangle(rect,c);
}

//draw a whole tile in a certain position
void draw_tile(double x,double y,double l,int itile,int orie=0)
{
  draw_exa(itile,x,y,l,tiles[itile].color);
  
  //add crown
  for(int icor=0;icor<6;icor++)
    if(tiles[itile].fc[icor]==tiles[itile].fc[(icor+5)%6])
      {
	int icorie=(icor+orie)%6;
	add_crown(tiles[itile].fc[icor]?red_crown:white_crown,x,y,l,icorie);
	add_arc(tiles[itile].fc[icor]?red_arc:white_arc,
		x+l*(corner_coords[icorie][0]+corner_coords[(icorie+5)%6][0])/2,
		y+l*(corner_coords[icorie][1]+corner_coords[(icorie+5)%6][1])/2,
		l/2,icorie*60-60,icorie*60+120);
      }
  //add arcs
  for(int icor=0;icor<6;icor++) add_arc(x,y,l,(icor+orie)%6,tiles[itile].fc[icor]);
  
  //add rectangles
  for(int icor=0;icor<6;icor++)
    {
      int a=tiles[itile].fc[(icor+5)%6];
      int b=tiles[itile].fc[(icor)];
      int c=tiles[itile].fc[(icor+1)%6];
      if((a==b)&&(a==c)) add_rectangle(b?red_rect[itile]:white_rect[itile],x,y,l,itile,(icor+orie)%6);
    }
  
  //add label
  labels[itile].SetTextAlign(22);
  char text[10];
  sprintf(text,"%d,%d",itile,orie);
  labels[itile].DrawText(x,y,text);
}

//draw a tile by position
void draw_tile(int ipos,int itile,int orie=0)
{
  double xmin=0,ymin=0;
  const double l=0.4;
  draw_tile(xmin+l*grid[ipos].shift[0],ymin+l*grid[ipos].shift[1]*R,l,itile,orie);
}

//draw the whole proposed configuration
void draw_conf(top_t top,int *orie)
{
  pad.Clear();
	  
  for(int ipos=0;ipos<19;ipos++) draw_tile(ipos,top[ipos],orie[ipos]);
  
  tela.Modified();
  tela.Update();      
}
