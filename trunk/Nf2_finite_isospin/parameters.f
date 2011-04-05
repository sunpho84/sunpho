c     include fixed defintions
      include "definitions.f" 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                          tunable parameters                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     lattice geometry
      integer nspa,nx,ny,nz,nt,nvol,nvolh
      include "lato"
      parameter(nx=nspa,ny=nspa,nz=nspa,nt=4) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number

c     color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c     debug flag: (0=no), (1=print debug info), (2=extreme verbose)
      integer debug
      parameter(debug=2)
      
c     algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi)  !algorithm for inverter: i_multi or i_singol



