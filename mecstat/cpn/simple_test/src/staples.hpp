#ifndef _STAPLES_HPP
#define _STAPLES_HPP

void ch_pot_site_staple(dcomplex *staple,dcomplex *z,dcomplex *l,int site);
void topo_staple(dcomplex &staple,dcomplex *l,int s,int mu);
void compute_topo_staples(dcomplex *staple,dcomplex *l);
void site_staple(dcomplex *staple,dcomplex *z,dcomplex *l,int site);
void link_staple(dcomplex &staple,dcomplex *z,int site,int mu);

#endif
