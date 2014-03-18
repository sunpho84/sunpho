#ifndef _UTILS_HPP
#define _UTILS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

int log2N(int N);
int factorize(int *list,int N);
void take_last_characters(char *out,const char *in,int size);
void fprintf_friendly_units(FILE *fout,int quant,int orders,const char *units);
void fprintf_friendly_filesize(FILE *fout,int quant);

#endif
