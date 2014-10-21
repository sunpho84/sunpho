#ifndef _DEBUG_HPP
#define _DEBUG_HPP

#include <chrono>
#include <cmath>
#include <ostream>

using namespace std;

void internal_crash(int line,const char *file,const char *templ,...);

class timing_t
{
  long long int tot=0;
  std::chrono::time_point<std::chrono::high_resolution_clock> start_moment;
public:
  void reset(){tot=0;}
  timing_t(){reset();}
  void start(){start_moment=std::chrono::high_resolution_clock::now();}
  void stop(){tot+=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-start_moment).count();}
  friend ostream &operator<<(ostream &of,timing_t &t);
};

inline ostream &operator<<(ostream &of,timing_t &t)
{of<<(double)t.tot*pow(10,-9)<<" s";return of;}

#endif
