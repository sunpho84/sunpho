#ifndef _STAPLES_HPP
#define _STAPLES_HPP

void topo_staple(dcomplex &staple,dcomplex *l,int s,int mu);
void site_staple(dcomplex *staple,dcomplex *z,dcomplex *l,int site);
void link_staple(dcomplex &staple,dcomplex *z,int site,int mu);

#endif
