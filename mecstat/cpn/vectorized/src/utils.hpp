#ifndef _UTILS_HPP
#define _UTILS_HPP

#define CRASH(...) crash(__LINE__,__FILE__,__VA_ARGS__)

void crash(int line,const char *file,const char *templ,...);

#endif
