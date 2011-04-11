      integer ncol,ncolm1,ncol2
      integer nperm
      integer nx,ny,nz,nt,nvol
      integer nvolh
      real pigr
      integer nmr
      integer n_quark_for_flavour
      integer unit_each

!     definitions of algorithm selection parameter
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
      
!     EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE
      parameter(EO=0,OE=1)

c      parameter(nx=8,ny=8,nz=8,nt=4) !lattice geometry
      parameter(nx=4,ny=4,nz=4,nt=4) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol) !color definition
      parameter(nperm=ncol**(ncol-3)+10000)
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(nmr=15) !maximum oredr for rational expansion
      parameter(pigr=3.141592654)

!     debug flag: (0=no), (1=print debug info), (2=extreme verbose)
      integer debug
      parameter(debug=0)
      
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan) !algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi (multishift)
