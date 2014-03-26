#ifndef _PER_SITE_NEIGHS_HPP
#define _PER_SITE_NEIGHS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <vector>

#include "geometry.hpp"

//holds the neighbors
class per_site_neighs_t : public std::vector<coords_t>
{
public:
  void add_neighbor(coords_t disp){this->push_back(disp);}
};

#endif
