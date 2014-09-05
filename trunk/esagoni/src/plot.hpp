#ifndef _PLOT_HPP
#define _PLOT_HPP

#include <TApplication.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraph.h>
#include <TPolyLine.h>
#include <TArc.h>
#include <TCrown.h>
#include <TText.h>

#include "tile.hpp"
#include "types.hpp"

extern TApplication myapp;
extern TCanvas tela;
extern TPad pad;
extern TPolyLine poly[19];
extern TArc red_arc,white_arc;
extern TPolyLine red_rect[19],white_rect[19];
extern TText labels[19];
extern TCrown red_crown,white_crown;

void init_graph();
void draw_conf(top_t top,int *orie);

#endif
