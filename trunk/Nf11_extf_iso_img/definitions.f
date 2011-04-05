cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                      internal fixed definitions                        c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $     ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $     ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $     ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $     ,u_rhmc=34,u_meas=35,u_magnet=36)

c     important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c     EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c     parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5)  !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15)       !maximum oredr for rational expansion

c     definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
