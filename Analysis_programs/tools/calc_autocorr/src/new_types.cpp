#include <algorithm>
#include <functional>
#include <iostream>
#include <cmath>

#include "debug.hpp"
#include "new_types.hpp"
#include "parse.hpp"

using namespace std;
using namespace std::placeholders;

//creator for var_t
var_t::var_t() : type(UNINIT_VAR) {}

//creator for single and double argument function
var_t::var_t(single_arg_function_t in) : type(SINGLE_ARG_FUNCTION),single_arg_function(in) {}
var_t::var_t(double_arg_function_t in) : type(DOUBLE_ARG_FUNCTION),double_arg_function(in) {}

void var_t::reset()
{
  //cout<<"Resetting var "<<name<<endl;
  if(type==DOUBLE_VECT_VAR) delete double_vect;
  type=UNINIT_VAR;
}

//copy from double vect
double_vect_t var_t::operator=(double_vect_t &in)
{
  reset();
  
  //cout<<"Using copy creator from double vect for var "<<name<<endl;
  
  type=DOUBLE_VECT_VAR;
  double_vect=new double_vect_t;
  
  (*double_vect)=in;
  
  return in;
}

//copy from double
double var_t::operator=(double in)
{
  reset();
  
  //cout<<"Using copy creator from double for var "<<name<<endl;
  
  type=DOUBLE_NUMB_VAR;
  double_numb=in;
  
  return in;
}

//crash if the two passed vector contains different number of entries
void check_comparable(double_vect_t &in1,double_vect_t &in2)
{
  int in1_size=in1.size();
  int in2_size=in2.size();
  if(in1_size!=in2_size)
    CRASH("asked to combine two vector having different number of entries (%d,%d)",in1_size,in2_size);
}

//addition
double_vect_t double_vect_t::operator+(double_vect_t in2)
{
  check_comparable(*this,in2);
  
  //transform
  double_vect_t out;out.resize(this->size());
  transform(this->begin(),this->end(),in2.begin(),out.begin(),plus<double>());return out;
}

//subtraction
double_vect_t double_vect_t::operator-(double_vect_t in2)
{
  check_comparable(*this,in2);
  
  //transform
  double_vect_t out;out.resize(this->size());
  transform(this->begin(),this->end(),in2.begin(),out.begin(),minus<double>());return out;
}

//multiplication
double_vect_t double_vect_t::operator*(double_vect_t in2)
{
  check_comparable(*this,in2);
  
  //transform
  double_vect_t out;out.resize(this->size());
  transform(this->begin(),this->end(),in2.begin(),out.begin(),multiplies<double>());return out;
}

//division
double_vect_t double_vect_t::operator/(double_vect_t in2)
{
  check_comparable(*this,in2);
  
  //transform
  double_vect_t out;out.resize(this->size());
  transform(this->begin(),this->end(),in2.begin(),out.begin(),divides<double>());return out;
}

//addition
double_vect_t double_vect_t::operator+(double in2)
{
  //transform
  double_vect_t out;out.resize(this->size());
  transform(this->begin(),this->end(),out.begin(),bind(plus<double>(),_1,in2));return out;
}

//subtraction
double_vect_t double_vect_t::operator-(double in2)
{
  //transform
  double_vect_t out;out.resize(this->size());
  transform(this->begin(),this->end(),out.begin(),bind(minus<double>(),_1,in2));return out;
}

//multiplication
double_vect_t double_vect_t::operator*(double in2)
{
  //transform
  double_vect_t out;out.resize(this->size());
  transform(this->begin(),this->end(),out.begin(),bind(multiplies<double>(),_1,in2));return out;
}

//division
double_vect_t double_vect_t::operator/(double in2)
{
  //transform
  double_vect_t out;out.resize(this->size());
  transform(this->begin(),this->end(),out.begin(),bind(divides<double>(),_1,in2));return out;
}

//unary minus
double_vect_t double_vect_t::operator-()
{
  double_vect_t out;out.resize(this->size());
  transform(this->begin(),this->end(),out.begin(),bind(minus<double>(),0.0,_1));return out;
}

//return the average of the vector
double double_vect_t::average()
{
  int size=this->size();
  double ave=0;
  for(int i=0;i<size;i++) ave+=(*this)[i];
  
  return ave/size;
}

//return the average of the vector
double double_vect_t::naive_error()
{
  int size=this->size();
  if(size<2) return 0;
  else
    {
      //compute ave and err
      double ave=0,err=0;
      for(int i=0;i<size;i++)
	{
	  double x=(*this)[i];
	  ave+=x;
	  err+=x*x;
	}
      
      ave/=size;
      err/=size;
      err-=ave*ave;
      err=sqrt(err/(size-1));
      
      return err;
    }
}

//addition
double_vect_t operator+(double in1,double_vect_t in2)
{
  //transform
  double_vect_t out;out.resize(in2.size());
  transform(in2.begin(),in2.end(),out.begin(),bind(plus<double>(),in1,_1));return out;
}

//subtraction
double_vect_t operator-(double in1,double_vect_t in2)
{
  //transform
  double_vect_t out;out.resize(in2.size());
  transform(in2.begin(),in2.end(),out.begin(),bind(minus<double>(),in1,_1));return out;
}

//multiplication
double_vect_t operator*(double in1,double_vect_t in2)
{
  //transform
  double_vect_t out;out.resize(in2.size());
  transform(in2.begin(),in2.end(),out.begin(),bind(multiplies<double>(),in1,_1));return out;
}

//division
double_vect_t operator/(double in1,double_vect_t in2)
{
  //transform
  double_vect_t out;out.resize(in2.size());
  transform(in2.begin(),in2.end(),out.begin(),bind(divides<double>(),in1,_1));return out;
}

//boh
double my_pow(double a,double b)
{return pow(a,b);}

//power
double_vect_t pow(double_vect_t in1,double in2)
{
  //transform
  double_vect_t out;out.resize(in1.size());
  transform(in1.begin(),in1.end(),out.begin(),bind(my_pow,_1,in2));return out;
}
