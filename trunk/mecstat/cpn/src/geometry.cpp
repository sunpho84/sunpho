#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "geometry.hpp"

//get coords of site
void geometry_t::coords_of_site(coords_t coords,int site,coords_t sizes)
{
  for(int mu=0;mu<NMU;mu++)
    {
      coords[mu]=site%sizes[mu];
      site/=sizes[mu];
    }
}

//return site of coords
int geometry_t::site_of_coords(coords_t coords,coords_t sizes)
{
  int site=0;
  for(int mu=0;mu<NMU;mu++) site=site*sizes[mu]+coords[mu];
  return site;
}

//initialize the geometry
void geometry_t::start(coords_t ext_glb_sizes)
{
  //copy glb sizes
  for(int mu=0;mu<NMU;mu++) glb_sizes[mu]=ext_glb_sizes[mu];
}
