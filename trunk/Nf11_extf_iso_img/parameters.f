c     include fixed defintions
#include "definitions.f" 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                          tunable parameters                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     set the quantization condition for the magnetic field
c     if relaxed=0 the field is proportional to 1/L
c     if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c     lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
#include "dimensioni.f"
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number
#ifdef ficp
      parameter(nvol_eo=nvol)
#else
      parameter(nvol_eo=nvolh)
#endif
c     color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c     algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi)  !algorithm for inverter: i_multi or i_singol



