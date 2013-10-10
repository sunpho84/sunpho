#ifndef _NEW_TYPES_H
#define _NEW_TYPES_H

#include <vector>
#include <string>
#include <map>

using namespace std;

typedef double (*single_arg_function_t) (double);
typedef double (*double_arg_function_t) (double,double);

//double vectors
struct double_vect_t : vector<double>
{
  double average();
  double naive_error();
  double_vect_t operator+(double_vect_t in);
  double_vect_t operator-(double_vect_t in);
  double_vect_t operator*(double_vect_t in);
  double_vect_t operator+(double in);
  double_vect_t operator-(double in);
  double_vect_t operator*(double in);
  double_vect_t operator/(double in);
  double_vect_t operator/(double_vect_t in);
  double_vect_t operator+=(double_vect_t in) {return (*this)=(*this)+in;};
  double_vect_t operator-=(double_vect_t in) {return (*this)=(*this)-in;};
  double_vect_t operator*=(double_vect_t in) {return (*this)=(*this)*in;};
  double_vect_t operator/=(double_vect_t in) {return (*this)=(*this)/in;};
  double_vect_t operator+=(double in) {return (*this)=(*this)+in;};
  double_vect_t operator-=(double in) {return (*this)=(*this)-in;};
  double_vect_t operator*=(double in) {return (*this)=(*this)*in;};
  double_vect_t operator/=(double in) {return (*this)=(*this)/in;};
  double_vect_t operator-();
};

double_vect_t operator+(double in1,double_vect_t in2);
double_vect_t operator-(double in1,double_vect_t in2);
double_vect_t operator*(double in1,double_vect_t in2);
double_vect_t operator/(double in1,double_vect_t in2);
double_vect_t pow(double_vect_t in1,double in2);

//structure holding vars
struct var_t
{
  void reset();
  int type;
  string name;
  double double_numb;
  double_vect_t* double_vect;
  single_arg_function_t single_arg_function;
  double_arg_function_t double_arg_function;
  var_t();
  var_t(single_arg_function_t in);
  var_t(double_arg_function_t in);
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
  driver_t();
  void init_scanner();
  void destroy_scanner();
};

#endif
