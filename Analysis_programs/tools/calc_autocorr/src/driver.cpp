#include <cmath>

#include "debug.hpp"
#include "new_types.hpp"
#include "parse.hpp"
#include "redefine_yy.hpp"

//extern
void parser_set_extra(driver_t*,void*);
int parser_lex_init(void**);
int parser_lex_destroy(void*);

//initializator
driver_t::driver_t(const char *path)
{
  //add function to the list
  var_map["pow"]=(double(*)(double,double))pow;
  var_map["atan2"]=(double(*)(double,double))atan2;
  var_map["sqrt"]=sqrt;
  var_map["sin"]=sin;
  var_map["cos"]=cos;
  var_map["tan"]=tan;
  var_map["log"]=log;
  var_map["asin"]=asin;
  var_map["acos"]=acos;
  var_map["atan"]=atan;
  var_map["exp"]=exp;
  
  //add doubles constant to the list
  var_map["M_PI"]=M_PI;
  
  //open file                                                                                                            
  fin=fopen(path,"r");
  if(fin==NULL) CRASH("opening %s",path);

  init_scanner();
}

//initialize the scanner
void driver_t::init_scanner()
{
  parser_lex_init(&scanner);
  parser_set_extra(this,scanner);
}

//destroy the scanner
void driver_t::destroy_scanner()
{
  parser_lex_destroy(scanner);
}
