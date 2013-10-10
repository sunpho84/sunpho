#ifndef _NEW_TYPES_H
#define _NEW_TYPES_H

#include <vector>
#include <string>
#include <map>

using namespace std;

//double vectors
typedef vector<double> double_vect_t;

//structure holding vars
struct var_t
{
  void reset();
  int type;
  string name;
  double double_numb;
  double_vect_t* double_vect;
  var_t();
  double_vect_t operator=(double_vect_t &in);
  double operator=(double in);
};

//list of var
typedef map<string,var_t> var_map_t;

//driver for parser
class driver_t
{
public:
  void *scanner;
  FILE *fin;
  var_map_t var_map;
  driver_t(const char *path);
  virtual ~driver_t(){destroy_scanner();}
protected:
  void init_scanner();
  void destroy_scanner();
};

//to be moved
double average(double_vect_t in);

#endif
