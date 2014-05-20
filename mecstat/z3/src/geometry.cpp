#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_GEOMETRY

#include "geometry.hpp"

//return the coordinate of a site
int site_of_coords(coords c)
{
  int site=0;
  for(int mu=0;mu<NDIMS;mu++) site=site*L[mu]+c[mu];
  return site;
}

//get coords of site
void coords_of_site(coords c,int site)
{
  for(int mu=NDIMS-1;mu>=0;mu--)
    {
      c[mu]=site%L[mu];
      site/=L[mu];
    }
}
