#ifndef _UTILS_HPP
#define _UTILS_HPP

#include "global_variables.hpp"

typedef int coords_t[nmu];
typedef int neighs_t[2*nmu];

void crash(const char *templ,...);
FILE *open_file(const char* path,const char* mod);
void coords_of_site(coords_t coords,int site);
int site_of_coords(coords_t x);

#endif
