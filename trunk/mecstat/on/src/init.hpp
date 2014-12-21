#ifndef _INIT_HPP
#define _INIT_HPP

#ifndef EXTERN_INIT
 #define EXTERN_INIT extern
#endif

void read_input(const char *path);
void init();

EXTERN_INIT int init_time;

#endif
