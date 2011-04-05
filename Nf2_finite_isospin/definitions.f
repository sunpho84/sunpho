cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                      internal fixed definitions                        c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $     ,u_magnet
      integer u_ferm_1,u_ferm_2,u_ferm_3,u_ferm_4

      parameter(u_ferm_1=11)
      parameter(u_ferm_2=12)
      parameter(u_ferm_3=13)
      parameter(u_ferm_4=14)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $     ,u_rhmc=34,u_meas=35,u_magnet=36)
      
C     lattice2
      character(LEN=30) lattice2
      parameter(lattice2="lattice2")

c     important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c     EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE
      parameter(EO=0,OE=1)

c     parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each
      integer nmr
      parameter(n_quark_for_flavour=2) !number of quark for flavour - keep 1
      parameter(unit_each=5)  !number of MD micro_step between reunitarization
      parameter(save_each=50000) !number of trajectory between two saved ones
      parameter(nmr=15)       !maximum oredr for rational expansion

c     definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
