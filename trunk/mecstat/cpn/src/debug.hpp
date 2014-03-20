#ifndef _DEBUG_HPP
#define _DEBUG_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <mpi.h>
#include <stdio.h>

#define MASTER_PRINTF(...) simul->master_fprintf(stdout,__VA_ARGS__)
#define CRASH_HARDLY(...) internal_crash(true,__LINE__,__FILE__,__VA_ARGS__)
#define CRASH_SOFTLY(...) internal_crash(false,__LINE__,__FILE__,__VA_ARGS__)
#define SHOUT(...) internal_shout(__LINE__,__FILE__,__VA_ARGS__)

#define MAX_VERBOSITY_LV 3

//add verbosity macro
#if MAX_VERBOSITY_LV>=1
#define VERBOSITY_LV1 (simul->verbosity_lv>=1)
#else
 #define VERBOSITY_LV1 0
#endif
#if MAX_VERBOSITY_LV>=2
#define VERBOSITY_LV2 (simul->verbosity_lv>=2)
#else
 #define VERBOSITY_LV2 0
#endif
#if MAX_VERBOSITY_LV>=3
#define VERBOSITY_LV3 (simul->verbosity_lv>=3)
#else
 #define VERBOSITY_LV3 0
#endif

//wrappers for verbosity_lv?
#define verbosity_lv1_master_printf(...) do{if(VERBOSITY_LV1) master_printf(__VA_ARGS__);}while(0)
#define verbosity_lv2_master_printf(...) do{if(VERBOSITY_LV2) master_printf(__VA_ARGS__);}while(0)
#define verbosity_lv3_master_printf(...) do{if(VERBOSITY_LV3) master_printf(__VA_ARGS__);}while(0)

void signal_handler(int sig);
void internal_crash(bool strength,int line,const char *file,const char *templ,...);
void internal_shout(int line,const char *file,const char *templ,...);

#endif
