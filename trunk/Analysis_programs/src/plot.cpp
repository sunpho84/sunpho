#pragma once

#define JACK
#define VTYPE jvec
#define TYPE jack
#define N njack
#include "plot_comm.cpp"
#undef N
#undef TYPE
#undef VTYPE
#undef JACK

#define BOOT
#define VTYPE bvec
#define TYPE boot
#define N nboot
#include "plot_comm.cpp"
#undef N
#undef TYPE
#undef VTYPE
#undef BOOT
