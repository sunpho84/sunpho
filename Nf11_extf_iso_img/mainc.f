# 1 "main.f"
# 1 "<built-in>"
# 1 "<command line>"
# 1 "main.f"
CC=====================================================
      program su_n_1_1_extf
CC
CC Program for simulating SU(3) gauge theories
CC with 1+1 staggered flavors in presence of imaginary
CC chemical potential coupled to each quark, a
CC static electromagnetic background field, and a real
CC chemical potential coupled to isospin
CC
CC The coordinate superindex is
CC i = ix + (iy-1)*nx + (iz-1)*nx*ny + (it-1)*nx*ny*nt
CC written by Massimo D'Elia
CC version 1.0 - 2007-2008
CC
CC=====================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 18 "main.f" 2

      integer time_to_run
      integer n_rand,n_rand_save
      integer nstep_md,n_traj
      integer acc,imposed_close
      integer termalizza,scala_ogni_step
      real tiniz,tfinal
      real dt_md,dt_md_save

      complex ieps,iepsq,ieps2q,ieps3q
      common/param2/dt_md,nstep_md,ieps,iepsq,ieps2q,ieps3q


c-------------------------------------------
c initialization operations
c-------------------------------------------

      call cpu_time(tiniz)

      call init(time_to_run,n_rand,n_traj,termalizza)

c-------------------------------------------------------
c START MONTE CARLO and MEASUREMENTS
c-------------------------------------------------------

      n_rand_save=n_rand
      dt_md_save=dt_md

 300 n_traj=n_traj+1

c decides if to termalize or not
      if(n_traj.le.termalizza) then
         n_rand=1
         scala_ogni_step=1
         dt_md=dt_md_save/3
         write(*,*) "Termalization ON, I'll scale at each microstep"
      else
         n_rand=n_rand_save
         scala_ogni_step=0
         dt_md=dt_md_save
      endif

c generate a new coniguration
      call rhmc_step(acc,n_traj,scala_ogni_step)

c save one configuration each: 'save_each'
c if(int(n_traj/save_each)*save_each.eq.n_traj) then
c FINIRE
c call write_lattice_save(n_traj)
c endif

c measure
      call measure(acc,n_traj,n_rand)

c decides if to close or to continue
      call cpu_time(tfinal)
      if(imposed_close().ne.1.and.(tfinal-tiniz).le.time_to_run) then
         goto 300
      endif

c-------------------------------------------


c-------------------------------------------------------
c closing operations
c-------------------------------------------------------

      call close_all(n_traj)

      write(*,*) "Execution time: ",int(tfinal-tiniz)," seconds"

c-------------------------------------------

      stop
      end

cc=========================================================================

# 1 "init.f" 1
c=========================================================
      subroutine load_par(init_flag,termalizza,time_to_run,n_rand
     $ ,immu_quark,immu_iso,remu,file_rhmc,val_ext_f)
c written by Sanfo
c load (by standard input) simulation parameters
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 9 "init.f" 2

c arguments
      integer init_flag
      integer termalizza
      integer time_to_run
      integer n_rand
      real immu_quark,immu_iso,remu
      character(LEN=30) file_rhmc
      real val_ext_f

c common-block passed variables
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue
      real dt_md
      integer nstep_md
      complex ieps,iepsq,ieps2q,ieps3q
      common/param2/dt_md,nstep_md,ieps,iepsq,ieps2q,ieps3q

c internal variables
      real dt
      character(LEN=10) cosa

! 0 -> cold; 1 -> hot; 2|3 -> stored
      read(*,'(a,i1)',ERR=201,END=201) cosa,init_flag
      if(cosa.eq."cho") goto 102
 201 write(*,*) "Errore riga 1 (cho)"
      stop

! 0 -> no, 1 -> yes
 102 read(*,'(a,i1)',ERR=202,END=202) cosa,termalizza
      if(cosa.eq."term") goto 103
 202 write(*,*) "Errore riga 2 (term)"
      stop

! coupling in the fundamental
 103 read(*,*,ERR=203,END=203) cosa,beta
      if(cosa.eq."beta") goto 104
 203 write(*,*) "Errore riga 3 (beta)"
      stop

! fermion mass
 104 read(*,*,ERR=204,END=204) cosa,mass
      if(cosa.eq."mass") goto 105
 204 write(*,*) "Errore riga 4 (mass)"
      stop

! total running time
 105 read(*,*,ERR=205,END=205) cosa,time_to_run
      if(cosa.eq."time") goto 106
 205 write(*,*) "Errore riga 5 (time)"
      stop

! number of md steps per trajectory
 106 read(*,'(a,i10)',ERR=206,END=206) cosa,nstep_md
      if(cosa.eq."nstep") goto 107
 206 write(*,*) "Errore riga 6 (nstep)"
      stop

! trajectory length
 107 read(*,*,ERR=207,END=207) cosa,dt
      if(cosa.eq."dt") goto 108
 207 write(*,*) "Errore riga 7 (dt)"
      stop

! stopping criterion for inverter
 108 read(*,*,ERR=208,END=208) cosa,residue
      if(cosa.eq."residue") goto 109
 208 write(*,*) "Errore riga 8 (residue)"
      stop

! number of random vectors
 109 read(*,'(a,i10)',ERR=209,END=209) cosa,n_rand
      if(cosa.eq."n_rand") goto 110
 209 write(*,*) "Errore riga 9 (n_rand)"
      stop

! imaginary isospin chemical potential/(pi T)
 110 read(*,*,ERR=210,END=210) cosa,immu_iso
      if(cosa.eq."immu_iso") goto 111
 210 write(*,*) "Errore riga 10 (immu_iso)"
      stop

! imaginary quark chemical potential/(pi T)
 111 read(*,*,ERR=211,END=211) cosa,immu_quark
      if(cosa.eq."immu_quark") goto 112
      write(*,*) "q",cosa,"q"
      call flush(6)
      stop
 211 write(*,*) "Errore riga 11 (immu_quark)"
      stop

! finite real isospin chemical potential/(pi T)
 112 read(*,*,ERR=212,END=212) cosa,remu
      if(cosa.eq."remu_iso") goto 113
 212 write(*,*) "Errore riga 12 (remu_iso)"
      stop

! external field value
 113 read(*,*,ERR=213,END=213) cosa,val_ext_f
c 113 read(*,'(a,f10.10)',ERR=213,END=213) cosa,val_ext_f
      if(cosa.eq."extf") goto 114
 213 write(*,*) "Errore riga 13 (extf)"
      stop

 114 continue
!! file for rhmc expansion



      file_rhmc="rhmc4"


      mass2=mass*mass
      dt_md=dt/nstep_md
# 131 "init.f"
      write(*,*)
      write(*,*) " ----Simulation parameters----"
      if(init_flag==0) then
         write(*,*) " Cold start"
      endif
      if(init_flag==1) then
         write(*,*) " Hot start"
      endif
      if(init_flag==2) then
         write(*,*) " Start from 'lattice' file"
      endif
      if(init_flag==3) then
         write(*,*) " Start from 'lattice2' file"
      endif
      write(*,*) " Coupling constant: ",beta
      write(*,*) " Mass, squared mass: ",mass,mass2
      write(*,*) " Will run for at least: ",time_to_run," seconds"
      write(*,*) " MD micro-step for each trajectory: ",nstep_md
      write(*,*) " Trajectories length: ",dt
      write(*,*) " MD micro-step length: ",dt_md
      write(*,*) " Residual for inverter: ",residue



      write(*,*) " Chiral measurement switched off at compilation time"






      write(*,*) " Isospin potential (in PI unit): ",immu_iso
      write(*,*) " Quark potential (in PI unit): ",immu_quark



      write(*,*)
     $ " Real isospin potential switched off at compilation time"
      if(remu.ne.0) then
         write(*,*) "Error: real isospin potential different from 0!"
         stop
      endif
# 183 "init.f"
      write(*,*) " Magnetic field switched off at compilation time"
      if(val_ext_f.ne.0) then
         write(*,*) "Error: external field different from 0!"
         stop
      endif

      write(*,*) " Rational expansion parameter file: ",file_rhmc

      write(*,*)
      write(*,*) " ----Algorithms----"
      select case (alg_inv)
      case (i_singol)
         write(*,*) " Inverter: normal CJ"
      case (i_multi)
         write(*,*) " Inverter: multi_shift"
      end select

      select case (alg_din)
      case (d_leapfrog)
         write(*,*) " MD integrator: leapfrog"
      case (d_omelyan)
         write(*,*) " MD integrator: omelyan"
      end select

      call flush(6)

      return
      end



c=========================================================
      subroutine load_rhmc(file_rhmc)
c written by Sanfo
c load rhmc parameters
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 221 "init.f" 2

c arguments
      character(LEN=30) file_rhmc

c common-block passed variables
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue

      integer nterm_a
      real cost_a,pole_a,coef_a
      common/rhmc_a/nterm_a,cost_a,pole_a(nmr),coef_a(nmr)
      integer nterm_h
      real cost_h,pole_h,coef_h
      common/rhmc_h/nterm_h,cost_h,pole_h(nmr),coef_h(nmr)

      real scala,min_rhmc,max_rhmc
      real ocost_a,opole_a,ocoef_a,z_a
      real ocost_h,opole_h,ocoef_h,z_h
      common/rhmc_o/scala,min_rhmc,max_rhmc,
     $ ocost_a,opole_a(nmr),ocoef_a(nmr),z_a,
     $ ocost_h,opole_h(nmr),ocoef_h(nmr),z_h


c internal variables
      integer i


cc-----------------------------
c LOAD RHMC PARAMETER
cc-----------------------------
      open(u_rhmc,file=file_rhmc,status='old')

! load extremes of validity of rhmc approximation
      read(u_rhmc,*) min_rhmc,max_rhmc

! load action part parameter
      read(u_rhmc,*) nterm_a
      read(u_rhmc,*) z_a
      read(u_rhmc,*) ocost_a
      cost_a=ocost_a
      do i = 1,nterm_a
         read(u_rhmc,*) opole_a(i)
         read(u_rhmc,*) ocoef_a(i)
         pole_a(i)=opole_a(i)
         coef_a(i)=ocoef_a(i)
      enddo

! load heat bath part
      read(u_rhmc,*) nterm_h
      read(u_rhmc,*) z_h
      read(u_rhmc,*) ocost_h
      cost_h=ocost_h
      do i = 1,nterm_h
         read(u_rhmc,*) opole_h(i)
         read(u_rhmc,*) ocoef_h(i)
         pole_h(i)=opole_h(i)
         coef_h(i)=ocoef_h(i)
      enddo
      close(u_rhmc)





      if((nterm_a.gt.nmr).or.(nterm_h.gt.nmr)) then
         write(*,*) "Attention, the order requested in the file"
         write(*,*) "is greater than the one defined in parameters.f"
         stop
      endif

      return
      end

c=========================================================
      subroutine print_rhmc()
c written by Sanfo
c print rhmc parameters
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 301 "init.f" 2

c common-block passed variables
      integer nterm_a
      real cost_a,pole_a,coef_a
      common/rhmc_a/nterm_a,cost_a,pole_a(nmr),coef_a(nmr)
      integer nterm_h
      real cost_h,pole_h,coef_h
      common/rhmc_h/nterm_h,cost_h,pole_h(nmr),coef_h(nmr)

      real scala,min_rhmc,max_rhmc
      real ocost_a,opole_a,ocoef_a,z_a
      real ocost_h,opole_h,ocoef_h,z_h
      common/rhmc_o/scala,min_rhmc,max_rhmc,
     $ ocost_a,opole_a(nmr),ocoef_a(nmr),z_a,
     $ ocost_h,opole_h(nmr),ocoef_h(nmr),z_h

      integer i

      write(*,*)
      write(*,*) "----Calcolo dell'azione e dinamica(-nf/4)----"
      write(*,*) " Numero termini:",nterm_a
      write(*,*) " Grado sviluppo:",z_a
      write(*,*) " Termine costante:",cost_a
      do i = 1,nterm_a
         write(*,*) " n°, polo, coef: ",i,pole_a(i),coef_a(i)
      enddo
      write(*,*)
      write(*,*) "----Heat Bath(nf/8)----"
      write(*,*) " Numero termini:",nterm_h
      write(*,*) " Grado sviluppo:",z_h
      write(*,*) " Termine costante:",cost_h
      do i = 1,nterm_h
         write(*,*) " n°, polo, coef: ",i,pole_h(i),coef_h(i)
      enddo
      write(*,*)


      return
      end

c=========================================================
      subroutine set_const(immu_quark,immu_iso,remu)
c written by Sanfo
c set various parameter regarding MD and chemical potentials
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 348 "init.f" 2

c arguments
      real immu_quark,immu_iso,remu

c common-block passed variables
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue
      integer ncfact,perm,sign
      common/perm/ncfact,perm(nperm,ncol),sign(nperm)
      complex potc1,potc2
      common/param3/potc1,potc2
      real dt_md
      integer nstep_md
      complex ieps,iepsq,ieps2q,ieps3q
      common/param2/dt_md,nstep_md,ieps,iepsq,ieps2q,ieps3q

c internal variables
      integer i
      real immu1,immu2,remu_fis

cc-----------------------------
cc compute ncol!
      ncfact=1
      do i=1,ncol
         ncfact=ncfact*i
      enddo
cc-----------------------------


cc-----------------------------
cc get back to single species imaginary chemical potentials
cc multiply by pi and divide by nt
      immu1=pigr*(immu_quark+immu_iso)/float(nt)
      immu2=pigr*(immu_quark-immu_iso)/float(nt)
      remu_fis=pigr*remu/float(nt)
cc-----------------------------

cc-----------------------------
cc set dt and sub-multiples
      iepsq=cmplx(0.0,dt_md/4.0)
      ieps2q=2.0*iepsq
      ieps3q=3.0*iepsq
      ieps=4.0*iepsq

      potc1=complex(+remu_fis,immu1)
      potc2=complex(-remu_fis,immu2)

cc-----------------------------


cc-----------------------------
cc print parameters
      write(*,*)
      write(*,*) "----Chemical potential in 'physical' unit----"
      write(*,*) " Chemical potential for U quark: ", immu1
      write(*,*) " Chemical potential for D quark: ", immu2
      write(*,*) " Real isospin chemical potential: ", remu_fis
      write(*,*)

cc-----------------------------

      return
      end

c=========================================================
      subroutine init(time_to_run,n_rand,n_traj,termalizza)
c written by Sanfo
c initialize the simulation
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 419 "init.f" 2

c arguments
      integer time_to_run
      integer n_rand
      integer n_traj

c common-block passed variables
      complex u,usave
      common/field/u(ncol,ncol,4,nvol)
      common/field1/usave(ncol,ncol,4,nvol)

c internal variables
      integer init_flag
      integer termalizza
      real immu_quark,immu_iso,remu
      character(LEN=30) file_rhmc
      real val_ext_f


      write(*,*)
      write(*,*) "----System parameter----"
      write(*,*) " Number of colors: ", ncol
      write(*,*) " Lattice size: ",nx,ny,nz,nt
      write(*,*)

cc-----------------------------
c READ INPUT & RHMC
cc-----------------------------

      call load_par(init_flag,termalizza,time_to_run,n_rand,immu_quark,
     $ immu_iso,remu,file_rhmc,val_ext_f)
      call load_rhmc(file_rhmc)
      call set_const(immu_quark,immu_iso,remu)

c-------------------------------------------

c-------------------------------------------
c INITIALIZATION OPERATION
c-------------------------------------------

      call ranstart !! initialize random number generator
      call generate_permutations !! table of ncol permutations
      call geometry !! set up lattice geometry
      call initialize_lattice(init_flag,n_traj) !! initialize lattice



      call addrem_stagphase
      call copyconf(usave,u) !! save configuration at first

c-------------------------------------------
c Observables file.
c
c *Explenation of content of observable files*
c -------------------------------------------------
c
c In each file the first colomn is the trajectory number.
c
c -meas_out:
c 2 acceptance of actual MD trajectory (0=discarded, 1 accepted)
c 3 spatial plaquette
c 4 temporal plaquette
c 5 real part of polyakov loop
c 6 imaginary part of polyakov loop
c
c -magnet_out:
c 2 real part of magnetic susceptibility
c 3 imag part
c
c Each "ferm" file, is a different observable, calculated via an
c approximated noisy estimator.
c
c In each ferm file the four columns are arranged as follow:
c 2 real part mean of the estimates over current config
c 3 imag part
c 4 squared mean of real part
c 5 squared mean of imag part
c
c The fist figure in the filename is relative to the quark (1=up, 2=down)
c
c The second figure is the observable, as follow:
c 1 chiral condensate
c 2 energy density
c 3 density of quark
c 4 pressure density
c 5 6 7 electic current
c


      open(u_meas,file='meas_out',status='unknown')
# 530 "init.f"
      return
      end
# 97 "main.f" 2
# 1 "generic_subroutines.f" 1
c============================================================================
c Control if close
c written by Sanfo
c============================================================================
      function imposed_close()

      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 9 "generic_subroutines.f" 2

      integer imposed_close

      open(u_control,file='control_file',status='unknown')
      read(u_control,*,end=312) imposed_close
 312 close(u_control)

      return
      end



c============================================================================
      subroutine close_all(n_traj)
c Call ending subroutines
c written by Massimo D'Elia
c============================================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 28 "generic_subroutines.f" 2

c arguments
      integer n_traj


      write(*,*) "ENDING JOB ..."

c save lattice
c call write_lattice(n_traj) !binary version
      call write_lattice2(n_traj,"lattice2") !text version, intersystem

c save random number generator status
      call ranfinish

c create an empty file "allright"
      open(u_allright,file='allright',status='unknown')
      close(u_allright)

      write(*,*) "JOB ENDED"

      return
      end


c============================================================================
      subroutine write_lattice(n_traj)
c Save link-array in a binary file
c written by Massimo D'Elia
c============================================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 59 "generic_subroutines.f" 2

c arguments
      integer n_traj

c common-block passed variables
      complex u
      common/field/u(ncol,ncol,4,nvol)

      call addrem_stagphase

      open(u_lattice,file='lattice',status='old',form='unformatted')
      write(u_lattice) u
      write(u_lattice) n_traj
      close(u_lattice)

      call addrem_stagphase

      return
      end

c============================================================================
      subroutine write_lattice2(n_traj,file_latt)
c Save link-array in a binary file
c written by Massimo D'Elia
c============================================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 86 "generic_subroutines.f" 2

c arguments
      integer n_traj
      character(LEN=30) file_latt

c common-block passed variables
      complex u
      common/field/u(ncol,ncol,4,nvol)

c internal variables
      integer icol1,icol2,idir,ivol
      complex c1

      call addrem_stagphase

      open(u_lattice,file=file_latt,status='unknown')
      write(u_lattice,*) n_traj
      do icol1=1,ncol
         do icol2=1,ncol
            do idir=1,4
               do ivol=1,nvol
                  c1=u(icol1,icol2,idir,ivol)
                  write(u_lattice,*) real(c1),aimag(c1)
               enddo
            enddo
         enddo
      enddo
      close(u_lattice)

      call addrem_stagphase

      return
      end

c=========================================================
      subroutine geometry
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
c set lattice geometry
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 129 "generic_subroutines.f" 2
      integer perm
      common/perm/ncfact,perm(nperm,ncol),sign(nperm)
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue
      integer sind,coor,forw,back,sum,parity,par
      integer forweo,backeo,forweo2,backeo2,sindeo,sindeoh,cooreo
      common/geo2/forweo2(nvol,4),backeo2(nvol,4)
      common/ge/sind(nx,ny,nz,nt),coor(nvol,4),forw(nvol,4),back(nvol,4)
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),
     $ parity(nx,ny,nz,nt),cooreo(nvol,4),forweo(nvol,4),backeo(nvol,4)
      common/gegeo/eotolex(nvol),lextoeo(nvol)
      common/stagphase/eta(4,nvol)
      complex u,usave,staple
      common/field/u(ncol,ncol,4,nvol)
      common/field1/usave(ncol,ncol,4,nvol)
      common/staple/staple(ncol,ncol,4,nvol)
      real gen_fact
      common/cartan/gen_fact(ncol)

c to begin let us define the factors multiplying
c the diagonal cartan generators, used in gauss.f

      do icol = 1,ncol-1
         x = float(icol*(icol + 1)/2)
         gen_fact(icol) = 1.0/sqrt(x)
      enddo


      do it = 1,nt
         do iz = 1,nz
            do iy = 1,ny
               do ix = 1,nx
                  sum = ix + iy + iz + it
                  par = sum - 2*(sum/2)
                  i = ix + (iy-1)*nx + (iz-1)*nx*ny + (it-1)*nx*ny*nz
                  j = (i-1)/2 + 1 + par*nvolh
                  sind(ix,iy,iz,it) = i
                  sindeo(ix,iy,iz,it) = j
                  parity(ix,iy,iz,it) = par
                  sindeoh(ix,iy,iz,it) = (i-1)/2 + 1
                  eotolex(sindeo(ix,iy,iz,it)) = sind(ix,iy,iz,it)
                  lextoeo(sind(ix,iy,iz,it)) = sindeo(ix,iy,iz,it)
                  coor(i,1) = ix
                  coor(i,2) = iy
                  coor(i,3) = iz
                  coor(i,4) = it
                  cooreo(j,1) = ix
                  cooreo(j,2) = iy
                  cooreo(j,3) = iz
                  cooreo(j,4) = it
               enddo
            enddo
         enddo
      enddo

      do it = 1,nt
         do iz = 1,nz
            do iy = 1,ny
               do ix = 1,nx
                  ixp = ix + 1
                  iyp = iy + 1
                  izp = iz + 1
                  itp = it + 1
                  if(ix.eq.nx) ixp = 1
                  if(iy.eq.ny) iyp = 1
                  if(iz.eq.nz) izp = 1
                  if(it.eq.nt) itp = 1
                  ixm = ix - 1
                  iym = iy - 1
                  izm = iz - 1
                  itm = it - 1
                  if(ixm.eq.0) ixm = nx
                  if(iym.eq.0) iym = ny
                  if(izm.eq.0) izm = nz
                  if(itm.eq.0) itm = nt
c funzioni spostamento per i siti cartesiani
                  forw(sind(ix,iy,iz,it),1) = sind(ixp,iy,iz,it)
                  forw(sind(ix,iy,iz,it),2) = sind(ix,iyp,iz,it)
                  forw(sind(ix,iy,iz,it),3) = sind(ix,iy,izp,it)
                  forw(sind(ix,iy,iz,it),4) = sind(ix,iy,iz,itp)
                  back(sind(ix,iy,iz,it),1) = sind(ixm,iy,iz,it)
                  back(sind(ix,iy,iz,it),2) = sind(ix,iym,iz,it)
                  back(sind(ix,iy,iz,it),3) = sind(ix,iy,izm,it)
                  back(sind(ix,iy,iz,it),4) = sind(ix,iy,iz,itm)
c funzioni spostamento per i siti even-odd
                  forweo(sindeo(ix,iy,iz,it),1) = sindeoh(ixp,iy,iz,it)
                  forweo(sindeo(ix,iy,iz,it),2) = sindeoh(ix,iyp,iz,it)
                  forweo(sindeo(ix,iy,iz,it),3) = sindeoh(ix,iy,izp,it)
                  forweo(sindeo(ix,iy,iz,it),4) = sindeoh(ix,iy,iz,itp)
                  backeo(sindeo(ix,iy,iz,it),1) = sindeoh(ixm,iy,iz,it)
                  backeo(sindeo(ix,iy,iz,it),2) = sindeoh(ix,iym,iz,it)
                  backeo(sindeo(ix,iy,iz,it),3) = sindeoh(ix,iy,izm,it)
                  backeo(sindeo(ix,iy,iz,it),4) = sindeoh(ix,iy,iz,itm)

                  forweo2(sindeo(ix,iy,iz,it),1) = sindeo(ixp,iy,iz,it)
                  forweo2(sindeo(ix,iy,iz,it),2) = sindeo(ix,iyp,iz,it)
                  forweo2(sindeo(ix,iy,iz,it),3) = sindeo(ix,iy,izp,it)
                  forweo2(sindeo(ix,iy,iz,it),4) = sindeo(ix,iy,iz,itp)
                  backeo2(sindeo(ix,iy,iz,it),1) = sindeo(ixm,iy,iz,it)
                  backeo2(sindeo(ix,iy,iz,it),2) = sindeo(ix,iym,iz,it)
                  backeo2(sindeo(ix,iy,iz,it),3) = sindeo(ix,iy,izm,it)
                  backeo2(sindeo(ix,iy,iz,it),4) = sindeo(ix,iy,iz,itm)

c definition of staggered phases
                  eta(1,sindeo(ix,iy,iz,it)) = 1.0
                  sum = ix - 1
                  eta(2,sindeo(ix,iy,iz,it)) = (-1.0)**sum
                  sum = ix - 1 + iy - 1
                  eta(3,sindeo(ix,iy,iz,it)) = (-1.0)**sum
                  sum = ix - 1 + iy - 1 + iz - 1
                  eta(4,sindeo(ix,iy,iz,it)) = (-1.0)**sum
                  if(it.eq.nt) then !! antiperiodic boundary conditions
                eta(4,sindeo(ix,iy,iz,it)) = -eta(4,sindeo(ix,iy,iz,it))
                  endif
               enddo
            enddo
         enddo
      enddo

      return
      end

c=========================================================
      subroutine initialize_lattice(init_flag,n_traj)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc 0 -> cold start; 1 -> hot start;
cc 2 -> bin stored config; 3 -> text stored config
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 261 "generic_subroutines.f" 2
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue
      integer sind,coor,forw,back,parity
      integer forweo,backeo,forweo2,backeo2,sindeo,sindeoh,cooreo
      common/geo2/forweo2(nvol,4),backeo2(nvol,4)
      common/ge/sind(nx,ny,nz,nt),coor(nvol,4),forw(nvol,4),back(nvol,4)
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),
     $ parity(nx,ny,nz,nt),cooreo(nvol,4),
     $ forweo(nvol,4),backeo(nvol,4)
      common/gegeo/eotolex(nvol),lextoeo(nvol)
      common/stagphase/eta(4,nvol)
      complex u
      common/field/u(ncol,ncol,4,nvol)
      integer idir,ind,ix,iy,iz,it,n_traj

      if(init_flag.eq.0) then
         ind = 0
         do it = 1,nt
            do iz = 1,nz
               do iy = 1,ny
                  do ix = 1,nx
                     ind = ind + 1
                     do idir = 1,4
                        call one(u(1,1,idir,ind))
                     enddo
                  enddo
               enddo
            enddo
         enddo
         n_traj = 0
      endif

      if(init_flag.eq.1) then
         ind = 0
         do it = 1,nt
            do iz = 1,nz
               do iy = 1,ny
                  do ix = 1,nx
                     ind = ind + 1
                     do idir = 1,4
                        call sun_random(u(1,1,idir,ind))
                     enddo
                  enddo
               enddo
            enddo
         enddo
         n_traj = 0
      endif

      if(init_flag.eq.2) then
         open(u_lattice,file='lattice',status='old',form='unformatted')
         read(u_lattice) u
         read(u_lattice) n_traj
         close(u_lattice)
      endif

      if(init_flag.eq.3) then
         open(u_lattice,file='lattice2',status='old')
         read(u_lattice,*) n_traj
         do icol1=1,ncol
            do icol2=1,ncol
               do idir=1,4
                  do ivol=1,nvol
                     read(u_lattice,*) ur,ui
                     u(icol1,icol2,idir,ivol)=cmplx(ur,ui)
                  enddo
               enddo
            enddo
         enddo
         close(u_lattice)
      endif

      if(init_flag.gt.3) then
         write(*,*) "BAD INIT_FLAG! use 0,1,2,3 (cold,hot,stored)"
      endif

      return
      end
# 508 "generic_subroutines.f"
c=========================================================
      subroutine addrem_stagphase()
c written by Massimo D'Elia
c version 1.0 - 2007/2008
c add or remove staggered phases
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 516 "generic_subroutines.f" 2

c common-block passed variables
      integer sindeo,sindeoh,parity,cooreo,forweo,backeo
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),
     $ parity(nx,ny,nz,nt),cooreo(nvol,4),
     $ forweo(nvol,4),backeo(nvol,4)
      real eta
      common/stagphase/eta(4,nvol)
      complex u
      common/field/u(ncol,ncol,4,nvol)

c internal variables
      integer idir,ix,iy,iz,it,icol1,icol2,eos

      do it=1,nt
         do iz=1,nz
            do iy=1,ny
               do ix=1,nx
                  eos=sindeo(ix,iy,iz,it)
                  do idir=1,4
                     do icol2=1,ncol
                        do icol1=1,ncol
                           u(icol1,icol2,idir,eos)=
     $ eta(idir,eos)*u(icol1,icol2,idir,eos)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
      end


c=========================================================
      subroutine normalize_lattice()
c written by Massimo D'Elia
c version 1.0 - 2007/2008
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 560 "generic_subroutines.f" 2
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue
      integer sind,coor,forw,back,parity
      integer forweo,backeo,forweo2,backeo2,sindeo,sindeoh,cooreo
      common/geo2/forweo2(nvol,4),backeo2(nvol,4)
      common/ge/sind(nx,ny,nz,nt),coor(nvol,4),forw(nvol,4),back(nvol,4)
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),
     $ parity(nx,ny,nz,nt),cooreo(nvol,4),
     $ forweo(nvol,4),backeo(nvol,4)
      common/gegeo/eotolex(nvol),lextoeo(nvol)
      common/stagphase/eta(4,nvol)
      complex u
      common/field/u(ncol,ncol,4,nvol)
      integer idir

         do ivol = 1,nvol
            do idir = 1,4
               call unitarize(u(1,1,idir,ivol))
            enddo
         enddo

      return
      end


c=========================================================
      subroutine normalize1_lattice()
c written by Massimo D'Elia
c version 1.0 - 2007/2008
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 593 "generic_subroutines.f" 2
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue
      integer sind,coor,forw,back,parity
      integer forweo,backeo,forweo2,backeo2,sindeo,sindeoh,cooreo
      common/geo2/forweo2(nvol,4),backeo2(nvol,4)
      common/ge/sind(nx,ny,nz,nt),coor(nvol,4),forw(nvol,4),back(nvol,4)
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),
     $ parity(nx,ny,nz,nt),cooreo(nvol,4),
     $ forweo(nvol,4),backeo(nvol,4)
      common/gegeo/eotolex(nvol),lextoeo(nvol)
      common/stagphase/eta(4,nvol)
      complex u
      common/field/u(ncol,ncol,4,nvol)
      integer idir

         do ivol = 1,nvol
            do idir = 1,4
               call unitarize1(u(1,1,idir,ivol))
            enddo
         enddo

      return
      end



c=========================================================
      subroutine generate_permutations
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc generate all permutations of the first ncol numbers and
cc their signs using an awful stupid algorithm
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 629 "generic_subroutines.f" 2
      integer perm
      common/perm/ncfact,perm(nperm,ncol),sign(nperm)

      integer perm1(ncol),perm11(ncol,ncol)
      integer lfact(ncol),lst(ncol),jref(ncol)
      real sign1(ncol)

      do i1 = 1,ncol
         perm1(i1) = i1
         do i2 = 1,ncol
            perm11(i1,i2) = i2
         enddo
         lst(i1) = 1
         sign1(i1) = 1.
         jref(i1) = 1
      enddo

      k2 = ncfact
      do j = 1,ncol
         lfact(j) = k2
c write(*,*) j,k2
         k2 = k2/(ncol - j + 1)
      enddo

      iperm = 1

 500 continue

      signaux = 1.0
      do k = 1,ncol
         perm(iperm,k) = perm1(k)
         signaux = signaux*sign1(k)
c write(*,*) iperm, signaux, perm1(k)
      enddo
      sign(iperm) = signaux


      do j = 2,ncol

         k2 = (iperm/lfact(j))*lfact(j) - iperm
         if(k2.eq.0) then
            do k = 1,ncol
               perm11(j,k) = perm11(jref(j),k)
            enddo
            sign1(j) = -1.0
            do k = j+1,ncol
               jref(k) = j
               sign1(k) = 1.0
            enddo
            k3 = perm11(j,j-1+lst(j))
            perm11(j,j-1+lst(j)) = perm11(j,j-1)
            perm11(j,j-1) = k3
            lst(j) = lst(j) + 1
            do k1 = j+1,ncol
               lst(k1) = 1
            enddo
            do k = 1,ncol
               perm1(k) = perm11(j,k)
            enddo
            goto 510
         endif


      enddo



 510 iperm = iperm + 1
c write(*,*) "----------------"
      if (iperm.lt.ncfact) goto 500

      signaux = 1.0
      do k = 1,ncol
         perm(iperm,k) = perm1(k)
         signaux = signaux * sign1(k)
c write(*,*) iperm, signaux, perm1(k)
      enddo
      sign(iperm) = signaux


      return
      end
c=========================================================

c============================================================================
c RANDOM NUMBER GENERATOR: standard ran2 from numerical recipes
c============================================================================
      function ran2()
      implicit real (a-h,o-z)
      implicit integer (i-n)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real ran2,am,eps,rnmx
      parameter(im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     & ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,
     & ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,
     & rnmx=1.-eps)
      integer idum2,j,k,iv,iy
      common /dasav/ idum,idum2,iv(ntab),iy
c save iv,iy,idum2
c data idum2/123456789/, iv/NTAB*0/, iy/0/

      if(idum.le.0) then
         idum=max0(-idum,1)
         idum2=idum
         do j=ntab+8,1,-1
            k=idum/iq1
            idum=ia1*(idum-k*iq1)-k*ir1
            if(idum.lt.0) idum=idum+im1
            if(j.le.ntab) iv(j)=idum
         enddo
         iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if(idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1) iy=iy+imm1
      ran2=min(am*iy,rnmx)
c write(*,*) "CHIAMATA",ran2
      return
      end

c=============================================================================
      subroutine ranstart
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 761 "generic_subroutines.f" 2

      common /dasav/ idum,idum2,iv(32),iy

      open(u_random, file='randomseed', status='unknown')
      read(u_random,*) idum
      read(u_random,*,end=117) idum2
      do i=1,32
         read(u_random,*) iv(i)
      enddo
      read(u_random,*) iy
      close(u_random)
      goto 118 !!takes account of the first start
 117 if(idum.ge.0) idum = -idum -1 !!
      close(u_random)
 118 continue !!

      return
      end

c=============================================================================
      subroutine ranfinish
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 785 "generic_subroutines.f" 2

      common /dasav/ idum,idum2,iv(32),iy

      open(u_random, file='randomseed', status='unknown')
      write(u_random,*) idum
      write(u_random,*) idum2
      do i=1,32
         write(u_random,*) iv(i)
      enddo
      write(u_random,*) iy
      close(u_random)

      return
      end
c=============================================================================
# 98 "main.f" 2
# 1 "sun_subroutines.f" 1
c=========================================================
      subroutine mmult(u,u1,u2)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc multiplication of two su(n) matrices: u = u1 * u2
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 10 "sun_subroutines.f" 2
      complex u(ncol,ncol),u1(ncol,ncol),u2(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
            u(i,j) = (0.,0.)
            do l = 1,ncol
               u(i,j) = u(i,j) + u1(i,l)*u2(l,j)
            enddo
         enddo
      enddo

      return
      end
c=========================================================
c=========================================================
      subroutine mmult_add(u,u1,u2)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc multiplication of two su(n) matrices: u = u1 * u2
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 33 "sun_subroutines.f" 2
      complex u(ncol,ncol),u1(ncol,ncol),u2(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
c u(i,j) = (0.,0.)
            do l = 1,ncol
               u(i,j) = u(i,j) + u1(i,l)*u2(l,j)
            enddo
         enddo
      enddo

      return
      end
c=========================================================
c=========================================================
      subroutine hmmult(u,u1,u2)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc multiplication of two su(n) matrices: u = u1~ * u2
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 56 "sun_subroutines.f" 2
      complex u(ncol,ncol),u1(ncol,ncol),u2(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
            u(i,j) = (0.,0.)
            do l = 1,ncol
               u(i,j) = u(i,j) + CONJG(u1(l,i))*u2(l,j)
            enddo
         enddo
      enddo

      return
      end
c=========================================================
c=========================================================
      subroutine hmmult_add(u,u1,u2)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc addition of multiplication of two su(n) matrices: u = u1~ * u2
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 79 "sun_subroutines.f" 2
      complex u(ncol,ncol),u1(ncol,ncol),u2(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
            do l = 1,ncol
               u(i,j) = u(i,j) + CONJG(u1(l,i))*u2(l,j)
            enddo
         enddo
      enddo

      return
      end
c=========================================================
c=========================================================
      subroutine mhmult(u,u1,u2)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc multiplication of two su(n) matrices: u = u1 * u2~
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 101 "sun_subroutines.f" 2
      complex u(ncol,ncol),u1(ncol,ncol),u2(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
            u(i,j) = (0.,0.)
            do l = 1,ncol
               u(i,j) = u(i,j) + u1(i,l)*CONJG(u2(j,l))
            enddo
         enddo
      enddo

      return
      end
c=========================================================

c=========================================================
      subroutine mhmult_add(u,u1,u2)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc addition of multiplication of two su(n) matrices: u = u + u1 * u2~
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 125 "sun_subroutines.f" 2
      complex u(ncol,ncol),u1(ncol,ncol),u2(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
            do l = 1,ncol
               u(i,j) = u(i,j) + u1(i,l)*CONJG(u2(j,l))
            enddo
         enddo
      enddo

      return
      end
c=========================================================

c=========================================================
      subroutine madd(u,u1,u2)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc addition of two su(n) matrices: u = u1 + u2
cc can be called with identical arguments
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 149 "sun_subroutines.f" 2
      complex u(ncol,ncol),u1(ncol,ncol),u2(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
               u(i,j) = u1(i,j) + u2(i,j)
         enddo
      enddo

      return
      end
c=========================================================


c=========================================================
      subroutine lincomb(u,a,u1,b,u2)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc addition of two su(n) matrices: u = a*u1 + b*u2
cc can be called with identical arguments
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 172 "sun_subroutines.f" 2
      complex u(ncol,ncol),u1(ncol,ncol),u2(ncol,ncol)
      real a,b

      do i = 1,ncol
         do j = 1,ncol
               u(i,j) = a*u1(i,j) + b*u2(i,j)
         enddo
      enddo

      return
      end
c=========================================================
c=========================================================
      subroutine rmult(u,a,u1)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc multiplication of su(n) matrix by real : u = a*u1
cc can be called with identical arguments
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 194 "sun_subroutines.f" 2
      complex u(ncol,ncol),u1(ncol,ncol)
      real a

      do i = 1,ncol
         do j = 1,ncol
               u(i,j) = a*u1(i,j)
         enddo
      enddo

      return
      end
c=========================================================

c=========================================================
      subroutine cmult(u,a,u1)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc multiplication of su(n) matrix by real : u = a*u1
cc can be called with identical arguments
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 217 "sun_subroutines.f" 2
      complex u(ncol,ncol),u1(ncol,ncol)
      complex a

      do i = 1,ncol
         do j = 1,ncol
               u(i,j) = a*u1(i,j)
         enddo
      enddo

      return
      end
c=========================================================

c=========================================================
      subroutine unitarize(u)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc project a (ncol,ncol) matrix onto SU(n)
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 239 "sun_subroutines.f" 2
      complex d,c(ncol),u(ncol,ncol)
      real rnorm

CC ORTONORMALIZATION ROW BY ROW

      do i = 1,ncol

c compute the scalar product with previous rows
         do j = 1,i-1
            c(j) = (0.,0.)
            do k = 1,ncol
               c(j) = c(j) + CONJG(u(j,k))*u(i,k)
            enddo
         enddo

c ortogonalize with respect to previous rows
         do j = 1,i-1
            do k = 1,ncol
               u(i,k) = u(i,k) - c(j)*u(j,k)
            enddo
         enddo

c normalize
         rnorm = 0.0
         do k = 1,ncol
            rnorm = rnorm + Real(u(i,k))**2 + aimag(u(i,k))**2
         enddo
         rnorm = 1./sqrt(rnorm)
         do k = 1,ncol
            u(i,k) = rnorm*u(i,k)
         enddo


      enddo !! close the loop on rows

c compute the determinant
      call det(d,u)
c correct last row to have determinant = 1
      do k = 1,ncol
         u(ncol,k) = conjg(d)*u(ncol,k)
      enddo

      return
      end
c=========================================================

c=========================================================
      subroutine unitarize1(u)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc project a (ncol,ncol) matrix onto SU(n)
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 294 "sun_subroutines.f" 2
      complex c(ncol),u(ncol,ncol)
      real rnorm

CC ORTONORMALIZATION ROW BY ROW

      do i = 1,ncol

c compute the scalar product with previous rows
         do j = 1,i-1
            c(j) = (0.,0.)
            do k = 1,ncol
               c(j) = c(j) + CONJG(u(j,k))*u(i,k)
            enddo
         enddo

c ortogonalize with respect to previous rows
         do j = 1,i-1
            do k = 1,ncol
               u(i,k) = u(i,k) - c(j)*u(j,k)
            enddo
         enddo

c normalize
         rnorm = 0.0
         do k = 1,ncol
            rnorm = rnorm + Real(u(i,k))**2 + aimag(u(i,k))**2
         enddo
         rnorm = 1./sqrt(rnorm)
         do k = 1,ncol
            u(i,k) = rnorm*u(i,k)
         enddo


      enddo !! close the loop on rows

c compute the determinant
c call det(d,u)
c correct last row to have determinant = 1
c do k = 1,ncol
c u(ncol,k) = conjg(d)*u(ncol,k)
c enddo

      return
      end
c=========================================================

c=========================================================
      subroutine det(d,u)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc compute the determinant of an (ncol,ncol) matrix
cc using the sum over permutations
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 350 "sun_subroutines.f" 2
      integer perm
      common/perm/ncfact,perm(nperm,ncol),sign(nperm)
      complex u(ncol,ncol),d,c

      d = cmplx(0.,0.)

      do i = 1,ncfact
        c = sign(i)
        do j = 1,ncol
           c = c * u(j,perm(i,j))
        enddo
        d = d + c
      enddo

      return
      end
c=========================================================

c=========================================================
      subroutine one(u)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc define the unit matrix
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 377 "sun_subroutines.f" 2
      complex u(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
            u(i,j) = cmplx(0.0,0.0)
         enddo
         u(i,i) = (1.0,0.0)
      enddo

      return
      end

c=========================================================
      subroutine zero(u)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc define the zero matrix
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 398 "sun_subroutines.f" 2
      complex u(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
            u(i,j) = cmplx(0.0,0.0)
         enddo
      enddo

      return
      end
c=========================================================
      subroutine equal(ua,ub)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc ua = ub for su(N) matrices
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 417 "sun_subroutines.f" 2
      complex ua(ncol,ncol),ub(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
            ua(i,j) = ub(i,j)
         enddo
      enddo

      return
      end

c=========================================================
      subroutine equalh(ua,ub)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc ua = ub~ for su(N) matrices
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 437 "sun_subroutines.f" 2
      complex ua(ncol,ncol),ub(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
            ua(i,j) = CONJG(ub(j,i))
         enddo
      enddo

      return
      end



c=========================================================
      subroutine ctrace(trace,ua)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc trace = complex trace of ua
c=========================================================

      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 460 "sun_subroutines.f" 2
      complex ua(ncol,ncol)
      complex trace

      trace = (0.,0.)
      do i = 1, ncol
         trace = trace + ua(i,i)
      enddo

      return
      end
c=========================================================
      subroutine rtrace(trace,ua)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc trace = real trace of ua
c=========================================================

      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 480 "sun_subroutines.f" 2
      complex ua(ncol,ncol)
      real trace

      trace = 0.
      do i = 1, ncol
         trace = trace + Real(ua(i,i))
      enddo

      return
      end

c=========================================================
      subroutine ctrace_uu(ctrace,ua,ub)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc ctrace = complex trace of ua * ub
c=========================================================

      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 501 "sun_subroutines.f" 2
      complex ua(ncol,ncol),ub(ncol,ncol)
      complex ctrace
      integer i,j

      ctrace = (0.,0.)
      do i = 1, ncol
         do j = 1,ncol
            ctrace = ctrace + ua(i,j)*ub(j,i)
         enddo
      enddo

      return
      end
c=========================================================
      subroutine ctrace_uuh(ctrace,ua,ub)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc ctrace = complex trace of ua * ub~
c=========================================================

      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 524 "sun_subroutines.f" 2
      complex ua(ncol,ncol),ub(ncol,ncol)
      complex ctrace
      integer i,j

      ctrace = (0.,0.)
      do i = 1, ncol
         do j = 1,ncol
            ctrace = ctrace + ua(i,j)*CONJG(ub(i,j))
         enddo
      enddo

      return
      end
c=========================================================
      subroutine rtrace_uu(rtrace,ua,ub)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc ctrace = real trace of ua * ub~
c=========================================================

      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 547 "sun_subroutines.f" 2
      complex ua(ncol,ncol),ub(ncol,ncol)
      complex ctrace
      real rtrace

      ctrace = (0.,0.)
      do i = 1, ncol
         do j = 1,ncol
            ctrace = ctrace + ua(i,j)*ub(j,i)
         enddo
      enddo
      rtrace = real(ctrace)
      return
      end

c=========================================================
      subroutine rtrace_uuh(rtrace,ua,ub)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc ctrace = real trace of ua * ub~
c=========================================================

      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 571 "sun_subroutines.f" 2
      complex ua(ncol,ncol),ub(ncol,ncol)
      complex ctrace
      real rtrace

      ctrace = (0.,0.)
      do i = 1, ncol
         do j = 1,ncol
            ctrace = ctrace + ua(i,j)*CONJG(ub(i,j))
         enddo
      enddo
      rtrace = real(ctrace)
      return
      end

c=========================================================
      subroutine TA(ua,ub)
c written by Massimo D'Elia
c version 1.0 - 01/2008
cc ua is the traceless antihermitean part of ub
c=========================================================

      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 595 "sun_subroutines.f" 2

c argomenti
      complex ua(ncol,ncol),ub(ncol,ncol)

c variabili interne
      complex ctrace

c parametri
      parameter(xinvncol=1.0/ncol)


      ctrace=cmplx(0,0)
      do i=1,ncol
         do j=1,ncol
            ua(i,j)=0.5*(ub(i,j)-conjg(ub(j,i)))
         enddo
         ctrace=ctrace+ua(i,i)
      enddo
      do i=1,ncol
         ua(i,i)=ua(i,i)-xinvncol*ctrace
      enddo

      return
      end
# 99 "main.f" 2
# 1 "sun_ext_subroutines.f" 1
c=========================================================
      subroutine expmat(u,h)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc exponential of su(n) subroutine up to sixth order
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 10 "sun_ext_subroutines.f" 2
      complex u(ncol,ncol),h(ncol,ncol),hp1(ncol,ncol),hp2(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
            u(i,j) = h(i,j)
         enddo
         u(i,i) = u(i,i) + (1.0,0.0)
      enddo

c (1 + x) done

      call mmult(hp1(1,1),h(1,1),h(1,1))

      do i = 1,ncol
         do j = 1,ncol
               u(i,j) = u(i,j) + 0.50000000000*hp1(i,j)
         enddo
      enddo

c second order done

      call mmult(hp2(1,1),hp1(1,1),h(1,1))

      do i = 1,ncol
         do j = 1,ncol
               u(i,j) = u(i,j) + 0.1666666667*hp2(i,j)
         enddo
      enddo

c third order done

      call mmult(hp1(1,1),hp2(1,1),h(1,1))

      do i = 1,ncol
         do j = 1,ncol
               u(i,j) = u(i,j) + 0.0416666667*hp1(i,j)
         enddo
      enddo

c fourth order done

      call mmult(hp2(1,1),hp1(1,1),h(1,1))

      do i = 1,ncol
         do j = 1,ncol
               u(i,j) = u(i,j) + 0.0083333333*hp2(i,j)
         enddo
      enddo

c fifth order done

      call mmult(hp1(1,1),hp2(1,1),h(1,1))

      do i = 1,ncol
         do j = 1,ncol
               u(i,j) = u(i,j) + 0.0013888889*hp1(i,j)
         enddo
      enddo

c sixth order done

      call unitarize1(u(1,1))

      return
      end
c=========================================================

c=========================================================
      subroutine sun_random(u_ran)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc create a random su(n) matrix by multiplying
cc n(n-1)/2 su(2) random matrices
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 86 "sun_ext_subroutines.f" 2
      complex u_l(ncol,ncol),u_ran(ncol,ncol),u_aux(ncol,ncol)

      integer i1,i2
      real u0,alpha,phi,ran2,sintheta,costheta,u1,u2,u3

CC generate i1,i2 indices of SU(2) subgroup

      call one(u_ran)
      do i1 = 1,ncol
         do i2 = i1+1,ncol
            call equal(u_aux,u_ran)

CC generate u0,u1,u2,u3 random on the four dim. sphere

            u0 = 1.0 - 2.0*ran2()
            alpha = sqrt(1 - u0**2)
            phi = 2.0*pigr*ran2()
            costheta = 1.0 - 2.0*ran2()
            sintheta = sqrt(1.0 - costheta**2)
            u3 = alpha*costheta
            u1 = alpha*sintheta*cos(phi)
            u2 = alpha*sintheta*sin(phi)

cc define u_l as unit matrix ...
            call one(u_l)

cc ... and then modify the elements in the chosen su(2) subgroup
            u_l(i1,i1) = CMPLX(u0,u3)
            u_l(i1,i2) = CMPLX(u2,u1)
            u_l(i2,i1) = CMPLX(-u2,u1)
            u_l(i2,i2) = CMPLX(u0,-u3)

            call mmult(u_ran,u_l,u_aux)

         enddo
      enddo

      return
      end
# 100 "main.f" 2
# 1 "action.f" 1
c=========================================================
      subroutine uaction(action,beta)
c written by Massimo D'Elia
c version 1.0 - 01/08
c local plaquette for molecular dynamics updating
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 9 "action.f" 2

c arguments
      real action(nvol),beta

c internal variables
      integer ivol
      real locplaq(nvol)


      call local_plaq(locplaq)

      do ivol=1,nvol
         action(ivol)=action(ivol)+locplaq(ivol)*beta
      enddo


      return
      end


c=========================================================
      subroutine qaction_flav(action,phi,potc,ud)
c written by Massimo D'Elia
c version 1.0 - 01/08
c local action for ud flavour
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 37 "action.f" 2

c arguments
      real action(nvol)
      complex phi(ncol,nvol_eo)
      complex potc
      integer ud

c common-block passed variables
      integer nterm
      real cost,pole,coef
      common/rhmc_a/nterm,cost,pole(nmr),coef(nmr)

c internal variables
      integer ivol,icol
      complex chi(ncol,nvol_eo)
      real loc_action

      call multi_shift_summed_inverter
     $ (chi,phi,
     $ potc,
     $ cost,pole,coef,
     $ nterm,ud)

      do ivol=1,nvol_eo
         loc_action=0
         do icol=1,ncol
            loc_action=
     $ loc_action+real(chi(icol,ivol)*conjg(phi(icol,ivol)))
         enddo
         action(ivol)=action(ivol)+loc_action
      enddo

      return
      end


c=========================================================
      subroutine haction(action,h)
c written by Massimo D'Elia
c version 1.0 - 01/08
c link momentum - part system action
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 81 "action.f" 2

c arguments
      real action(nvol)
      complex h(ncol,ncol,4,nvol)

c internal variables
      integer ivol,idir
      real rtr
      real loc_action


      do ivol=1,nvol
         loc_action=0
         do idir=1,4
            call rtrace_uu(rtr,h(1,1,idir,ivol),h(1,1,idir,ivol))
            loc_action=loc_action+rtr
         enddo
         action(ivol)=action(ivol)+loc_action*0.5
      enddo


      return
      end




c=========================================================
      subroutine action(az)
c written by Massimo D'Elia
c version 1.0 - 01/08
c calculate locally the total action of the system
c IMPORTANT: the action have to be calculated locally
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 117 "action.f" 2

c arguments
      real az(nvol)

c common-block passed variables
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue
      complex p_w,ipdot
      common/momenta/p_w(ncol,ncol,4,nvol),ipdot(ncol,ncol,4,nvol)
      complex phi1,phi2,chi1,chi2
      common/pseudof1/phi1(ncol,nvol_eo),chi1(ncol,nvol_eo)
      common/pseudof2/phi2(ncol,nvol_eo),chi2(ncol,nvol_eo)
      complex potc1,potc2
      common/param3/potc1,potc2

      integer ivol

      do ivol=1,nvol
         az(ivol)=0
      enddo

      call uaction(az,beta)
      call haction(az,p_w)
      call qaction_flav(az,phi1,potc1,1)
      call qaction_flav(az,phi2,potc2,2)


      return
      end



c=========================================================
      subroutine dif_action(dif,az_2,az_1)
c written by Massimo D'Elia
c version 1.0 - 01/08
c calculate the global differenxe between two locally-calculated action
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 157 "action.f" 2

c arguments
      real dif
      real az_2(nvol),az_1(nvol)

c internal variables
      integer ivol


      dif=0
      do ivol=1,nvol
         dif=dif+az_2(ivol)-az_1(ivol)
      enddo


      return
      end
# 101 "main.f" 2
# 1 "dirac_matrix.f" 1
c===================================================
      subroutine m2d(y,x,potc,mass)
c written by Sanfo
c calculate y=(m^dag m)*x
c===================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 8 "dirac_matrix.f" 2

c arguments
      complex y(ncol,nvol_eo),x(ncol,nvol_eo)
      complex potc
      real mass

c internal variables
      integer ivol,icol
      complex h(ncol,nvolh)
      real mass2

      mass2=mass**2

c even part
      call D(OE, h,x,potc)
      call D(OE_DAG,y,h,potc)




      do ivol=1,nvolh
         do icol=1,ncol
            y(icol,ivol)=mass2*x(icol,ivol)+y(icol,ivol)



         enddo
      enddo
# 52 "dirac_matrix.f"
      return
      end

c=========================================================
      subroutine D(eooe,w,v,potc)
c written by Massimo D'Elia, unified by Zumbo
c apply the even/odd part or the odd/even part of the Dirac
c Matrix, depending on the value of the passed variable eooe
c OE and EO variable defined in "parameters.f"
c version 1.0 -
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 65 "dirac_matrix.f" 2

c arguments
      integer eooe
      complex w(ncol,nvolh),v(ncol,nvolh)
      complex potc

c common blocks
      integer sindeo,sindeoh,parity,cooreo,forweo,backeo
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),parity(nx,ny
     $ ,nz,nt),cooreo(nvol,4), forweo(nvol,4),backeo(nvol,4)
      integer forweo2,backeo2
      common/geo2/forweo2(nvol,4),backeo2(nvol,4)
      complex u
      common/field/u(ncol,ncol,4,nvol)

c internal variables
      integer base_ivol,ivol,icol,idir
      complex w1(ncol)

      complex w2(ncol),w3(ncol)

      complex A,B,C

      A=0
      B=0
      C=0

      select case (eooe)

      case(EO,OE)
         A=0.5
         B=+0.5*exp(+potc)
         C=-0.5*exp(-potc)
      case(OE_DAG,EO_DAG)
         A=-0.5
         B=-0.5*exp(-conjg(potc))
         C=+0.5*exp(+conjg(potc))
      case(EO_P_OE_DAG,OE_P_EO_DAG)
         A=0
         B=+0.5*(exp(potc)-exp(-conjg(potc)))
         C=+0.5*(exp(conjg(potc))-exp(-potc))
      end select

      do base_ivol=1,nvolh

c select which part of D to apply
         if(eooe.eq.EO.or.eooe.eq.OE_DAG.or.eooe.eq.EO_P_OE_DAG) then
            ivol=base_ivol
         else
            ivol=base_ivol+nvolh
         endif

         if(eooe.ne.EO_P_OE_DAG.and.eooe.ne.OE_P_EO_DAG) then

c reset temporary variable
            do icol=1,ncol
               w1(icol)=cmplx(0,0)
            enddo




            do idir=1,3

               call vmult_add(w1,u(1,1,idir,ivol),v(1,forweo(ivol
     $ ,idir)))
               call vhmult_sub(w1,u(1,1,idir,backeo2(ivol,idir)), v(1
     $ ,backeo(ivol,idir)))
            enddo

         endif

         call vmult(w2,u(1,1,4,ivol),v(1,forweo(ivol,4)))
         call vhmult(w3,u(1,1,4,backeo2(ivol,4)),v(1,backeo(ivol,4)))


c write final value
         do icol=1,ncol
            w(icol,base_ivol)=A*w1(icol)

     $ +B*w2(icol)+C*w3(icol)

         enddo

      enddo

      return
      end
# 102 "main.f" 2
# 1 "force.f" 1
c=========================================================
      subroutine compute_utpdt(iepsatt)
c written by Massimo D'Elia
c version 1.0 - 2007/2008
c evolve all the links for a "iepsatt" MD time
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 9 "force.f" 2

c arguments
      complex iepsatt

c common blocks
      complex u,p_w,ipdot
      common/field/u(ncol,ncol,4,nvol)
      common/momenta/p_w(ncol,ncol,4,nvol),ipdot(ncol,ncol,4,nvol)

c internal variables
      complex h(ncol,ncol),u1(ncol,ncol),u2(ncol,ncol)
      integer ivol,idir

      do ivol=1,nvol
         do idir=1,4
            call cmult(h,iepsatt,p_w(1,1,idir,ivol))
            call expmat(u1,h)
            call mmult(u2,u1,u(1,1,idir,ivol))
            call equal(u(1,1,idir,ivol),u2)
         enddo
      enddo !! end of starting step

      return
      end


c=========================================================
      subroutine compute_ptpdt(iepsatt,scale_expansion)
c calculate p_w=p_w-iepsatt*ipdot to evolve link momenta
c i.e calculate v(t+dt)=v(t)+a*dt
c written by Massimo D'Elia
c version 1.0 - 2007/2008
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 44 "force.f" 2

c arguments
      complex iepsatt
      integer scale_expansion

c common blocks
      complex p_w,ipdot
      common/momenta/p_w(ncol,ncol,4,nvol),ipdot(ncol,ncol,4,nvol)

c internal variables
      integer ivol,idir,icol1,icol2


      call compute_ipdot(scale_expansion)

      do ivol=1,nvol
         do idir=1,4
            do icol1=1,ncol
               do icol2=1,ncol
                  p_w(icol1,icol2,idir,ivol)=p_w(icol1,icol2,idir,ivol)-
     $ iepsatt*ipdot(icol1,icol2,idir,ivol)
               enddo
            enddo
         enddo
      enddo

      return
      end

c=========================================================
      subroutine compute_ipdot(scale_expansion)
c written by Massimo D'Elia
c compute the force a(t)
c version 1.0 - 2007/2008
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 81 "force.f" 2

c arguments
      integer scale_expansion

c common blocks
      complex p_w,ipdot
      common/momenta/p_w(ncol,ncol,4,nvol),ipdot(ncol,ncol,4,nvol)
      complex u
      common/field/u(ncol,ncol,4,nvol)

c internal variables
      complex ipdot_g(ncol,ncol,4,nvol),ipdot_q(ncol,ncol,4,nvol)
      integer ivol,idir,icol1,icol2
      complex u1(ncol,ncol),u2(ncol,ncol)


      if(scale_expansion.eq.1) then
         call scale_rhmc
      endif





      call compute_gluonic_force(ipdot_g)
      call compute_quark_force(ipdot_q)

c sum gluonic and quark force and calculate Trace Anti-hermitian part
      do ivol=1,nvol
         do idir=1,4
            do icol1=1,ncol
               do icol2=1,ncol
                  u1(icol1,icol2)
     $ =ipdot_g(icol1,icol2,idir,ivol)
     $ +ipdot_q(icol1,icol2,idir,ivol)
               enddo
            enddo
            call mmult(u2(1,1),u(1,1,idir,ivol),u1(1,1))
            call TA(ipdot(1,1,idir,ivol),u2(1,1))
         enddo
      enddo
# 131 "force.f"
      return
      end


c=========================================================
      subroutine compute_gluonic_force(ipdot)
c written by Sanfo
c calculte gluonic part of the force
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 142 "force.f" 2

c arguments
      complex ipdot(ncol,ncol,4,nvol)

c common blocks
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue
      complex staple
      common/staple/staple(ncol,ncol,4,nvol)

c internal variables
      integer ivol,idir
      real r_1

      r_1=beta/float(ncol)

      call compute_staples
      do ivol=1,nvol
         do idir=1,4
            call rmult(ipdot(1,1,idir,ivol),r_1,staple(1,1,idir,ivol))
         enddo
      enddo

      return
      end

c=========================================================
      subroutine compute_quark_force(ipdot)
c written by Sanfo
c calculate total quark force
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 175 "force.f" 2

c arguments
      complex ipdot(ncol,ncol,4,nvol)

c common blocks
      integer forweo,backeo,sindeo,sindeoh,cooreo,parity
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),
     $ parity(nx,ny,nz,nt),cooreo(nvol,4),
     $ forweo(nvol,4),backeo(nvol,4)
      complex potc1,potc2
      common/param3/potc1,potc2
      complex phi1,phi2,fuf1,fuf2
      common/pseudof1/phi1(ncol,nvol_eo),fuf1(ncol,nvol_eo)
      common/pseudof2/phi2(ncol,nvol_eo),fuf2(ncol,nvol_eo)
      real pole,coef,cost
      integer nterm
      common/rhmc_a/nterm,cost,pole(nmr),coef(nmr)






c internal variables
      integer ie,ivol,icol,idir,iterm,icol1,icol2
      complex chi1_t(ncol,nvol_eo,nmr),chi2_t(ncol,nvol_eo,nmr)
      complex w11(ncol),w12(ncol)
      complex w21(ncol),w22(ncol)
      complex v1_o(ncol,nvolh,nmr),v2_o(ncol,nvolh,nmr)
      complex c1,c2
      complex exp_1,exp_2
# 218 "force.f"
      call multi_shift_inverter
     $ (chi1_t,phi1,
     $ potc1,
     $ pole,
     $ nterm,1)

      call multi_shift_inverter
     $ (chi2_t,phi2,
     $ potc2,
     $ pole,
     $ nterm,2)






      exp_1=exp(+potc1)
      exp_2=exp(+potc2)





      do iterm=1,nterm


! calculation of force, without or with EO improvement
         call D(OE,v1_o(1,1,iterm),chi1_t(1,1,iterm),potc1)




      enddo







      do iterm=1,nterm


! calculation of force, without or with EO improvement
         call D(OE,v2_o(1,1,iterm),chi2_t(1,1,iterm),potc2)





      enddo





      do ivol=1,nvol
         ie=ivol-nvolh
         do idir=1,4






            do icol1=1,ncol
               do icol2=1,ncol
                  ipdot(icol1,icol2,idir,ivol)=0
               enddo
            enddo

            do iterm=1,nterm
               do icol=1,ncol
                  if(ivol.le.nvolh) then
                     w11(icol)=v1_o(icol,forweo(ivol,idir),iterm)
                     w21(icol)=v2_o(icol,forweo(ivol,idir),iterm)
                     w12(icol)=CONJG(chi1_t(icol,ivol,iterm))
                     w22(icol)=CONJG(chi2_t(icol,ivol,iterm))
# 307 "force.f"
                  else
                     w11(icol)=-chi1_t(icol,forweo(ivol,idir),iterm)!minus sign
                     w21(icol)=-chi2_t(icol,forweo(ivol,idir),iterm)
                     w12(icol)=CONJG(v1_o(icol,ie,iterm)) !! implemented here
                     w22(icol)=CONJG(v2_o(icol,ie,iterm))
# 320 "force.f"
                  endif


                  if(idir.eq.4) then
# 333 "force.f"
                        w11(icol)=exp_1*w11(icol)
                        w21(icol)=exp_2*w21(icol)





                  endif

               enddo

               do icol1=1,ncol
# 356 "force.f"
                  c1=w11(icol1)
                  c2=w21(icol1)







                  do icol2=1,ncol
                     ipdot(icol1,icol2,idir,ivol)=ipdot(icol1,icol2,idir
     $ ,ivol)+coef(iterm)*(c1*w12(icol2)+c2*w22(icol2)




     $ )
                  enddo
               enddo
            enddo

         enddo
      enddo

      return
      end
# 103 "main.f" 2
# 1 "updating_subroutines.f" 1
c=========================================================
      subroutine rhmc_step(acc,n_traj,scale_each_step)
c written by Sanfo
c version 1.0 - 04/08
c exec singolo rhmc step
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 9 "updating_subroutines.f" 2

c arguments
      integer acc
      integer n_traj
      integer scale_each_step

c common-block passed variables
      complex u,usave
      common/field/u(ncol,ncol,4,nvol)
      common/field1/usave(ncol,ncol,4,nvol)
      complex ieps,iepsq,ieps2q,ieps3q
      real dt_md
      integer nstep_md
      common/param2/dt_md,nstep_md,ieps,iepsq,ieps2q,ieps3q
c RHMC ACTION PART
      real az_old,az_new
      common/azione_loc/az_old(nvol),az_new(nvol)

c internal variables
      real difact,edifact
      real ran2,xrand
      integer temp







      open(u_random, file='n_steps_tune', status='unknown')
      read(u_random,*,end=111,err=111) temp
      nstep_md=temp
      go to 112
 111 write(u_random,*) nstep_md
      go to 113

 112 xrand=0
      do temp=1,tune_nstep_each
         read(u_random,*,end=113,err=113) edifact
         xrand=xrand+edifact
      enddo
      xrand=xrand/tune_nstep_each
      close(u_random)
      write(*,*) "CONTRO"
      if(xrand.gt.0.80) nstep_md=nstep_md-1
      if(xrand.lt.0.70) nstep_md=nstep_md+1
      open(u_random, file='n_steps_tune', status='unknown')
      write(u_random,*) nstep_md

 113 call create_momenta
      call create_phi

      call action(az_old)
      call dinamica(scale_each_step)
      call action(az_new)

      call dif_action(difact,az_new,az_old)


! metropolis test
      acc=1
      if(difact.gt.0.and.scale_each_step.eq.0) then
         edifact=exp(-difact)
         xrand=ran2()
      else
         edifact=1
         xrand=0
      endif

      if(edifact.gt.xrand) then
         call copyconf(usave,u) !! accept configuration
         acc=1
      else
         call copyconf(u,usave) !! reject configuration
         acc=0
      endif

      write(*,*)n_traj,difact,edifact,xrand,acc,nstep_md


      if(scale_each_step.eq.0) write(u_random,*) edifact
      close(u_random)

c call write_lattice2(n_traj)

      return
      end


c=========================================================
      subroutine create_momenta()
c written by Massimo D'Elia
c version 1.0 -
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 105 "updating_subroutines.f" 2

c common-block passed variables
      complex p_w,ipdot
      common/momenta/p_w(ncol,ncol,4,nvol),ipdot(ncol,ncol,4,nvol)

c internal variables
      integer ivol,idir

c internal parameters
      real sigma
      parameter(sigma=1)

      do ivol=1,nvol
         do idir=1,4
            call gauss_matrix(p_w(1,1,idir,ivol),sigma)
         enddo
      enddo

      return
      end

c=========================================================
      subroutine create_phi_flav(phi,potc,ud)
c written by Massimo D'Elia
c version 1.0 -
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 133 "updating_subroutines.f" 2

c arguments
      complex phi(ncol,nvol_eo)
      complex potc
      integer ud

c common-block passed variables
      integer nterm
      real cost,pole,coef
      common/rhmc_h/nterm,cost,pole(nmr),coef(nmr)

c internal variables
      complex rnd(ncol,nvol_eo)
      integer ivol

c internal parameters
      real sigma
      parameter(sigma=1)

      do ivol=1,nvol_eo
         call gauss_vector(rnd(1,ivol),sigma)
      enddo

      call multi_shift_summed_inverter
     $ (phi,rnd,
     $ potc,
     $ cost,pole,coef,
     $ nterm,ud)

      return
      end

c=========================================================
      subroutine create_phi
c written by Massimo D'Elia
c version 1.0 -
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 172 "updating_subroutines.f" 2

c common-block passed variables
      complex potc1,potc2
      common/param3/potc1,potc2
      complex phi1,phi2,chi1,chi2
      common/pseudof1/phi1(ncol,nvol_eo),chi1(ncol,nvol_eo)
      common/pseudof2/phi2(ncol,nvol_eo),chi2(ncol,nvol_eo)

      call scale_rhmc

      call create_phi_flav(phi1,potc1,1)
      call create_phi_flav(phi2,potc2,2)

      return
      end



c=========================================================
      subroutine leapfrog(scale_each_step)
c written by Massimo D'Elia
c version 1.0 - 2007/2008
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 197 "updating_subroutines.f" 2

c arguments
      integer scale_each_step

c common-block passed variables
      complex ieps,iepsq,ieps2q,ieps3q
      real dt_md
      integer nstep_md
      common/param2/dt_md,nstep_md,ieps,iepsq,ieps2q,ieps3q

c internal variables
      integer istep





c----------------------------------------------------------------------------
c PRIMO MEZZO STEP: TROVA U(-dt/2)
c----------------------------------------------------------------------------

c CALCOLA U(-dt/2) ovvero p(-dt/2)=p(0)-v(0)*dt/2
c NB: v(0)=p_w è stato estratto con create_momenta
      call compute_utpdt(ieps2q)

c normalize after first half step
      call addrem_stagphase
      call normalize1_lattice
      call addrem_stagphase

c----------------------------------------------------------------------------
c CICLO PRINCIPALE
c----------------------------------------------------------------------------

      do istep=1,nstep_md-1 !! starting main M.D. loop






c CALCOLA P(t+dt) overo v(t+dt)=v(t)+a(t)*dt
         call compute_ptpdt(ieps,scale_each_step)

c CALCOLA U(t+dt/2) ovvero p(t+dt/2)=p(t-dt/2)+v(t)*dt
         call compute_utpdt(ieps)

c NORMALIZZA OGNI unit_each STEP
         if((istep-unit_each*(istep/unit_each)).eq.0) then
            call addrem_stagphase
            call normalize1_lattice
            call addrem_stagphase
         endif

      enddo !! enddo over M.D. steps minus one

c----------------------------------------------------------------------------
c ULTIMO MEZZO STEP: RITROVA U(dt)
c----------------------------------------------------------------------------

c CALCOLA P(t+dt)
      call compute_ptpdt(ieps,scale_each_step)

c CALCOLA U(t'=t+dt)
      call compute_utpdt(ieps2q)

      return
      end


c=========================================================
      subroutine omelyan(scale_each_step)
c written by Sanfo
c version 1.0 - 2007/2008
c questo è una modifica del leapfrog, per info vedi
c cond-mat/0110438v1
c per ogni singolo step dovrei calcolare
c
c v1 = v(t) + a[r(t)]*lambda*dt
c r1 = r(t) + v1*dt/2
c v2 = v1 + a[r1]*(1 -2*lambda)*dt
c r(t + dt) = r1 + v2*dt/2
c v(t + h) = v2 + a[r(t + dt)]*lambda*dt
c
c ma invece ottimizzo un poco mettendo insieme la prima e l'ultima riga
c per fare questo devo distinguere il primo e l'ultimo step
c
c solo 1°step:
c v1 = v(t) + a[r(t)]*lambda*dt
c
c r1 = r(t) + v1*dt/2
c v2 = v1 + a[r1]*(1 -2*lambda)*dt
c r(t + dt) = r1 + v2*dt/2
c
c non l'ultimo step: 
c v1 = v2 + a[r(t + dt)]*2*lambda*dt
c
c ultimo step:
c v(t + h) = v2 + a[r(t + dt)]*lambda*dt
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 299 "updating_subroutines.f" 2

c arguments
      integer scale_each_step

c common-block passed variables
      complex ieps,iepsq,ieps2q,ieps3q
      real dt_md
      integer nstep_md
      common/param2/dt_md,nstep_md,ieps,iepsq,ieps2q,ieps3q

c internal variables
      complex ldt,dldt,udldt
      integer istep

c parameters
      real lambda
      parameter(lambda=0.1931833)


      ldt=ieps*lambda
      dldt=2*ldt
      udldt=ieps-dldt





c----------------------------------------------------------------------------
c PRIMO MEZZO STEP: TROVA P(t+lambda*dt)
c----------------------------------------------------------------------------

c CALCOLA P(t+lambda*dt) ovvero v1=v(t)+a[r(t)]*lambda*dt
      call compute_ptpdt(ldt,scale_each_step)

c normalize after first half step
      call addrem_stagphase
      call normalize1_lattice
      call addrem_stagphase

c----------------------------------------------------------------------------
c CICLO PRINCIPALE
c----------------------------------------------------------------------------

      do istep=1,nstep_md !! starting main M.D. loop





c CALCOLA U(t+dt/2) ovvero r1=r(t)+v1*dt/2
         call compute_utpdt(ieps2q)

c CALCOLA P(t+(1-2*lambda)*dt) overo v2=v1+a[r1]*(1-2*lambda)*dt
         call compute_ptpdt(udldt,scale_each_step)

c CALCOLA U(t+dt/2) ovvero r(t+dt)=r1+v2*dt/2
         call compute_utpdt(ieps2q)

         if(istep.eq.nstep_md) then
c CALCOLA P(t+dt) ovvero v(t+dt)=v2+a[r(t+dt)]*lambda*dt
            call compute_ptpdt(ldt,scale_each_step)
         else
c CALCOLA P(t+dt) ovvero v1'=v2+a[r(t+dt)]*2*lambda*dt
            call compute_ptpdt(dldt,scale_each_step)
         endif

c NORMALIZZA OGNI unit_each STEP
         if((istep-unit_each*(istep/unit_each)).eq.0) then
            call addrem_stagphase
            call normalize1_lattice
            call addrem_stagphase
         endif

      enddo

      return
      end


c=========================================================
      subroutine dinamica(scale_each_step)
c written by Sanfo
c chiama l'opportuna funzione di update
c=========================================================
      implicit none

# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 386 "updating_subroutines.f" 2

c arguments
      integer scale_each_step

c common blocks
      complex ieps,iepsq,ieps2q,ieps3q
      real dt_md
      integer nstep_md
      common/param2/dt_md,nstep_md,ieps,iepsq,ieps2q,ieps3q

      dt_md=1.0/nstep_md
      iepsq=cmplx(0.0,dt_md/4.0)
      ieps2q=2.0*iepsq
      ieps3q=3.0*iepsq
      ieps=4.0*iepsq

      if(scale_each_step.eq.1) then
         call leapfrog(scale_each_step)
      else
         select case (alg_din)
         case (d_leapfrog)
            call leapfrog(scale_each_step)
         case (d_omelyan)
            call omelyan(scale_each_step)
         CASE DEFAULT
            write(*,*) "Integratore sconosciuto"
            stop
         end select
      end if

c----------------------------------------------------------------------------
c NORMALIZZAZIONE FINALE
c----------------------------------------------------------------------------

! normalize before ending trajectory
      call addrem_stagphase
      call normalize_lattice
      call addrem_stagphase

      call scale_rhmc

      return
      end



c=========================================================
      subroutine copyconf(u_out,u_in)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 439 "updating_subroutines.f" 2

c arguments
      complex u_in(ncol,ncol,4,nvol),u_out(ncol,ncol,4,nvol)

c internal variables
      integer icol1,icol2,idir,ivol

      do icol1=1,ncol
        do icol2=1,ncol
          do idir=1,4
            do ivol=1,nvol
              u_out(icol1,icol2,idir,ivol)=u_in(icol1,icol2,idir,ivol)
            enddo
          enddo
        enddo
      enddo

      return
      end
# 104 "main.f" 2
# 1 "measure_subroutines.f" 1
c=========================================================
      subroutine measure(acc,n_traj,n_rand)
c written by Sanfo
c version 1.0 - 04/2008
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 8 "measure_subroutines.f" 2

c arguments
      integer acc
      integer n_traj
      integer n_rand

c common blocks
      complex potc1,potc2
      common/param3/potc1,potc2

c internal variables
      real plaq_s,plaq_t,plaq2_s,plaq2_t
      complex ploop
# 38 "measure_subroutines.f"
      n_rand=n_rand

c measure section
      call plaquette(plaq_s,plaq_t,plaq2_s,plaq2_t)
      call polyakov(ploop)
# 55 "measure_subroutines.f"
c the file writing is collapsed in order to minimize possible failure
      write(u_meas,*) n_traj,acc,plaq_s,plaq_t,real(ploop),aimag(ploop)
# 87 "measure_subroutines.f"
      return
      end


c=========================================================
      subroutine local_plaq(locplaq)
c written by Massimo D'Elia
c version 1.0 - 01/08

c calculate the plaquette site by site. Used to avoid numerical
c rounding during the calculation of the diference between new and
c old gauge action (so do not change even if it seems unoptimized)
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 102 "measure_subroutines.f" 2

c arguments
      real locplaq(nvol)

c common blocks
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue
      integer forweo2,backeo2
      common/geo2/forweo2(nvol,4),backeo2(nvol,4)
      complex u
      common/field/u(ncol,ncol,4,nvol)

c internal variables
      integer i,mu,nu
      complex v1(ncol,ncol),v2(ncol,ncol)
      complex ctr
      real ctr_tot

c parameters
      real xinvncol
      parameter(xinvncol=1.0/ncol)


      do i=1,nvol

         ctr_tot=0

         do mu=1,3

            do nu=mu+1,4 !loop on orthogonal directions

               call mmult(v1,u(1,1,nu,i),u(1,1,mu,forweo2(i,nu)))
               call mmult(v2,u(1,1,mu,i),u(1,1,nu,forweo2(i,mu)))
               call ctrace_uuh(ctr,v1,v2)

               ctr_tot=ctr_tot+real(ctr)

            enddo !nu

            locplaq(i)=6+xinvncol*ctr_tot !include segno fasi stag

         enddo !! mu

      enddo !! i

      return
      end



c=========================================================
      subroutine compute_staples
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
cc compute staples to be used in the updating
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 160 "measure_subroutines.f" 2

c common blocks
      integer sind,coor,forw,back
      common/ge/sind(nx,ny,nz,nt),coor(nvol,4),forw(nvol,4),back(nvol,4)
      integer forweo2,backeo2
      common/geo2/forweo2(nvol,4),backeo2(nvol,4)
      complex u
      common/field/u(ncol,ncol,4,nvol)
      complex staple
      common/staple/staple(ncol,ncol,4,nvol)

c internal variables
      integer ind,idir,jdir
      complex v1(ncol,ncol)

      do ind=1,nvol

         do idir=1,4

            call zero(staple(1,1,idir,ind))
c loop on orthogonal directions
            do jdir=1,4

               if(jdir.ne.idir) then

c staple in the forward direction
                  call mmult(v1,u(1,1,jdir,ind),u(1,1,idir,forweo2(ind
     $ ,jdir)))
                  call mhmult_add(staple(1,1,idir,ind),u(1,1,jdir
     $ ,forweo2(ind,idir)),v1)
c staple in the backward direction
                  call mmult(v1,u(1,1,idir,backeo2(ind,jdir)),u(1,1
     $ ,jdir,forweo2(backeo2(ind,jdir),idir)))
                  call hmmult_add(staple(1,1,idir,ind),v1,u(1,1,jdir
     $ ,backeo2(ind,jdir)))

               endif

            enddo

         enddo !! idir
      enddo !! ivol

      return
      end

c=========================================================
      subroutine plaquette(plaq_s,plaq_t,plaq2_s,plaq2_t)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
c measures fundamental and adjoint plaquette
c=========================================================

      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 215 "measure_subroutines.f" 2

c arguments
      real plaq_s,plaq_t,plaq2_s,plaq2_t

c common blocks
      integer forweo2,backeo2
      common/geo2/forweo2(nvol,4),backeo2(nvol,4)
      complex u
      common/field/u(ncol,ncol,4,nvol)

c local variables
      integer ind,idir,jdir
      complex v1(ncol,ncol),v2(ncol,ncol)
      complex ctr


      ind=0
      plaq_s=0.0 !! spatial and temporal plaquette
      plaq_t=0.0
      plaq2_s=0.0 !! spatial and temporal plaquette squared
      plaq2_t=0.0

      do ind =1,nvol

         do idir=1,3

c loop on orthogonal directions

            do jdir=idir+1,4

               call mmult(v1,u(1,1,jdir,ind),u(1,1,idir,forweo2(ind
     $ ,jdir)))
               call mmult(v2,u(1,1,idir,ind),u(1,1,jdir,forweo2(ind
     $ ,idir)))
               call ctrace_uuh(ctr,v1,v2)

               if (jdir.eq.4) then
                  plaq_t= plaq_t+Real(ctr)
                  plaq2_t= plaq2_t+Real(ctr)**2+aimag(ctr)**2
               else
                  plaq_s=plaq_s+Real(ctr)
                  plaq2_s=plaq2_s+Real(ctr)**2+aimag(ctr)**2
               endif

            enddo !! jdir: closing do over orthogonal directions

         enddo !! idir: closing do over directions

      enddo !! ind

      plaq_s=-plaq_s/(ncol*3.0*nvol) !! minus sign due to
      plaq_t=-plaq_t/(ncol*3.0*nvol) !! staggered phases
      plaq2_s=plaq2_s/(ncol*ncol*3.0*nvol)
      plaq2_t=plaq2_t/(ncol*ncol*3.0*nvol)

      return
      end

c=========================================================
      subroutine polyakov(poly_c)
c written by Massimo D'Elia
c version 1.0 - 11/07/2004
c=========================================================

      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 281 "measure_subroutines.f" 2

c arguments
      complex poly_c

c common blocks
      integer sindeo,sindeoh,parity,cooreo,forweo,backeo
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),
     $ parity(nx,ny,nz,nt),cooreo(nvol,4),forweo(nvol,4),backeo(nvol,4)
      complex u
      common/field/u(ncol,ncol,4,nvol)

c internal variables
      integer ix,iy,iz,it,idir
      complex v1(ncol,ncol),poly_mat(ncol,ncol)
      complex ctr


      poly_c=cmplx(0.0,0.0)

      do iz=1,nz
         do iy=1,ny
            do ix=1,nx

               call one(poly_mat)
               idir=4

               do it=1,nt
                  call mmult(v1,poly_mat,u(1,1,idir,sindeo(ix,iy,iz
     $ ,it)))
                  call equal(poly_mat,v1)
               enddo

               call ctrace(ctr,poly_mat)
               poly_c=poly_c+ctr

            enddo !! ix
         enddo !! iy
      enddo !! iz

c minus due to staggered phase
      poly_c=-poly_c/float(ncol*nx*ny*nz)

      return
      end
# 105 "main.f" 2
# 1 "vecn_subroutines.f" 1
c=========================================================
      subroutine vmult(v,u,w)
c written by Massimo D'Elia
c version 1.0 -
cc multiplication of an su(n) matrices by a su(n) vector: v = u * w
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 10 "vecn_subroutines.f" 2
      complex u(ncol,ncol),v(ncol),w(ncol)

      do i = 1,ncol
         v(i) = (0.,0.)
         do l = 1,ncol
            v(i) = v(i) + u(i,l)*w(l)
         enddo
      enddo

      return
      end
c=========================================================
c=========================================================
      subroutine vhmult(v,u,w)
c written by Massimo D'Elia
cc multiplication of an su(n) matrices by a su(n) vector: v = u~ * w
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 30 "vecn_subroutines.f" 2
      complex u(ncol,ncol),v(ncol),w(ncol)

      do i = 1,ncol
         v(i) = (0.,0.)
         do l = 1,ncol
            v(i) = v(i) + CONJG(u(l,i))*w(l)
         enddo
      enddo

      return
      end
c=========================================================
c=========================================================
      subroutine vmult_add(v,u,w)
c written by Massimo D'Elia
c version 1.0 -
cc add the multiplication of an su(n) matrices by a su(n) vector: v = u * w
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 51 "vecn_subroutines.f" 2
      complex u(ncol,ncol),v(ncol),w(ncol)

      do i = 1,ncol
         do l = 1,ncol
            v(i) = v(i) + u(i,l)*w(l)
         enddo
      enddo

      return
      end
c=========================================================
c=========================================================
      subroutine vhmult_add(v,u,w)
c written by Massimo D'Elia
cc add the multiplication of an su(n) matrices by a su(n) vector: v = u~ * w
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 70 "vecn_subroutines.f" 2
      complex u(ncol,ncol),v(ncol),w(ncol)

      do i = 1,ncol
         do l = 1,ncol
            v(i) = v(i) + CONJG(u(l,i))*w(l)
         enddo
      enddo

      return
      end
c=========================================================
c=========================================================
      subroutine vmult_sub(v,u,w)
c written by Massimo D'Elia
c version 1.0 -
cc add the multiplication of an su(n) matrices by a su(n) vector: v = u * w
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 90 "vecn_subroutines.f" 2
      complex u(ncol,ncol),v(ncol),w(ncol)

      do i = 1,ncol
         do l = 1,ncol
            v(i) = v(i) - u(i,l)*w(l)
         enddo
      enddo

      return
      end
c=========================================================
c=========================================================
      subroutine vhmult_sub(v,u,w)
c written by Massimo D'Elia
cc add the multiplication of an su(n) matrices by a su(n) vector: v = u~ * w
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 109 "vecn_subroutines.f" 2
      complex u(ncol,ncol),v(ncol),w(ncol)

      do i = 1,ncol
         do l = 1,ncol
            v(i) = v(i) - CONJG(u(l,i))*w(l)
         enddo
      enddo

      return
      end
c=========================================================
c=========================================================
      subroutine tensn(u,v,w)
c written by Massimo D'Elia
cc tensor product of two sun vectors
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 128 "vecn_subroutines.f" 2
      complex u(ncol,ncol),v(ncol),w(ncol)

      do i = 1,ncol
         do j = 1,ncol
            u(i,j) = v(i) * w(j)
         enddo
      enddo

      return
      end
c=========================================================
c=========================================================
      subroutine tensn_add(u,v,w)
c written by Massimo D'Elia
cc tensor product of two sun vectors
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 147 "vecn_subroutines.f" 2
      complex u(ncol,ncol),v(ncol),w(ncol)

      do i = 1,ncol
         do j = 1,ncol
            u(i,j) = u(i,j) + v(i) * w(j)
         enddo
      enddo

      return
      end
c=========================================================
# 106 "main.f" 2
# 1 "inverter.f" 1
c=========================================================
      subroutine singol_inverter(chi,phi,potc,pole)
c written by Massimo D'Elia
c version 1.0 - 2007/2008
c chi = (MM~)^-1 phi
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 10 "inverter.f" 2

c argomenti
      complex chi(ncol,nvol_eo),phi(ncol,nvol_eo)
      complex potc
      real pole

c variabili passate da common block
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue

c variabili interne
      complex s(ncol,nvol_eo),r(ncol,nvol_eo),p(ncol,nvol_eo)
      complex c1
      real delta,omega,lambda,gammag
      integer iter,riter

c parametri
      parameter(sigma=1)
      parameter(niter=2000,rniter=5)
      parameter(one_q=0.25)

      do ivol=1,nvol_eo
         do icol=1,ncol
            chi(icol,ivol)=cmplx(0,0)
         enddo
      enddo

      riter=0

 1200 continue !! loop over congrad to verify truerest

      call m2d(s,chi,potc,pole)

      delta=0
      do ivol=1,nvol_eo
         do icol=1,ncol
            c1=phi(icol,ivol)-s(icol,ivol)
            p(icol,ivol)=c1
            r(icol,ivol)=c1
            delta=delta+real(c1)**2+aimag(c1)**2
         enddo
      enddo

!!-------------------------------------------------------------------------

      iter=0

 1201 continue

      iter=iter+1

      call m2d(s,p,potc,pole)

      alpha=0
      do ivol=1,nvol_eo
         do icol=1,ncol
            alpha=alpha+real(s(icol,ivol)*conjg(p(icol,ivol)))
         enddo
      enddo

      omega=delta/alpha
      lambda=0

      do ivol=1,nvol_eo
         do icol=1,ncol
            chi(icol,ivol)=chi(icol,ivol)+omega*p(icol,ivol)
            c1=r(icol,ivol)-omega*s(icol,ivol)
            r(icol,ivol)=c1
            lambda=lambda+real(c1)**2+aimag(c1)**2
         enddo
      enddo

      gammag=lambda/delta
      delta=lambda
      do ivol=1,nvol_eo
         do icol=1,ncol
            p(icol,ivol)=r(icol,ivol)+gammag*p(icol,ivol)
         enddo
      enddo

      if(lambda.gt.residue.and.iter.lt.niter) goto 1201

      call m2d(s,chi,potc,pole)

      lambda=0
      do ivol=1,nvol_eo
         do icol=1,ncol
            c1=phi(icol,ivol)-s(icol,ivol)
            lambda=lambda+real(c1)**2+aimag(c1)**2
         enddo
      enddo

      riter=riter+1

      if(lambda.gt.residue.and.riter.lt.rniter) goto 1200

      return
      end



c=========================================================
      subroutine multi_shift_inverter(x,b,potc,pole,nterm,ud)
c written by Sanfo
c chiama il corretto inverter
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 118 "inverter.f" 2

c argomenti
      complex x(ncol,nvol_eo,nmr),b(ncol,nvol_eo)
      complex potc
      real pole(nmr)
      integer nterm
      integer ud


      if(alg_inv.eq.i_singol) then
         call multi_shift_fuffa_inverter(x,b,potc,pole,nterm,ud)
      endif

      if(alg_inv.eq.i_multi) then
         call multi_shift_true_inverter(x,b,potc,pole,nterm,ud)
      endif

      return
      end


c=========================================================
      subroutine multi_shift_true_inverter(x,b,potc,pole,nterm,ud)
c written by Sanfo
c x_s = (m+X)^-1 b
c lo calcola usando il multi-shift conjugate gradient
c in cui X=-0.25*DOE*DEO
c la variabile cost contiene la costante additiva
c il vettore pole contiene i poli
c il vettore coef contiene i coefficienti
c nterm sono i termini effettivi (max passato dal file parameters.f)
c NB: la massa fisica è già inclusa nei poli
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 153 "inverter.f" 2

c argomenti
      complex x(ncol,nvol_eo,nmr),b(ncol,nvol_eo)
      complex potc
      real pole(nmr)
      integer nterm
      integer ud

c variabili passate da common block
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue

c variabili interne
      complex s(ncol,nvol_eo),r(ncol,nvol_eo),p(ncol,nvol_eo)
      complex ps(ncol,nvol_eo,nmr)
      real zps(nmr),zas(nmr),zfs(nmr),betas(nmr),alphas(nmr)
      complex c1
      real rr,rfrf,pap,alpha
      real betap,betaa
      integer iterm
      integer icol,ivol
      integer iter
      integer flag(nmr)

c parametri
      real one_q
      integer niter
      parameter(one_q=0.25)
      parameter(niter=2000)




      ud=ud !to avoid warning






! x=0
      do iterm=1,nterm
         do ivol=1,nvol_eo
            do icol=1,ncol
               x(icol,ivol,iterm)=0
            enddo
          enddo
          flag(iterm)=1
      enddo

! -p=b
! -r=b
! -calcolo Rr=(r,r)
      rr=0
      do ivol=1,nvol_eo
         do icol=1,ncol
            c1=b(icol,ivol)
            p(icol,ivol)=c1
            r(icol,ivol)=c1
            rr=rr+real(c1)**2+aimag(c1)**2
         enddo
      enddo

! -betaa=1
      betaa=1

! -ps=b
! -zps=zas=1
! -alphas=0
      do iterm=1,nterm
         do ivol=1,nvol_eo
            do icol=1,ncol
               ps(icol,ivol,iterm)=b(icol,ivol)
            enddo
         enddo
         zps(iterm)=1
         zas(iterm)=1
         alphas(iterm)=0
      enddo

! -alpha=0
      alpha=0

!----------------------------------------------------------------------

      do iter=1,niter

! -s=Ap
         call m2d(s,p,potc,mass)

! -pap=(p,s)=(p,Ap)
         pap=0
         do ivol=1,nvol_eo
            do icol=1,ncol
               pap=pap+real(s(icol,ivol)*conjg(p(icol,ivol)))
            enddo
         enddo

! calcola betaa=rr/pap=(r,r)/(p,Ap)
         betap=betaa
         betaa=-rr/pap

! calcola
! -zfs
! -betas
! -x
         do iterm=1,nterm
         if(flag(iterm).eq.1) then
            zfs(iterm)=
     $ zas(iterm)*zps(iterm)*betap/
     $ (betaa*alpha*(zps(iterm)-zas(iterm))+
     $ zps(iterm)*betap*(1-pole(iterm)*betaa))
            betas(iterm)=betaa*zfs(iterm)/zas(iterm)
            do ivol=1,nvol_eo
               do icol=1,ncol
                  x(icol,ivol,iterm)=x(icol,ivol,iterm)-
     $ betas(iterm)*ps(icol,ivol,iterm)
               enddo
            enddo
         endif
         enddo

! calcolo
! -r'=r+betaa*s=r+beta*Ap
! -rfrf=(r',r')
         rfrf=0
         do ivol=1,nvol_eo
            do icol=1,ncol
               c1=r(icol,ivol)+betaa*s(icol,ivol)
               r(icol,ivol)=c1
               rfrf=rfrf+real(c1)**2+aimag(c1)**2
            enddo
         enddo

! calcola alpha=rfrf/rr=(r',r')/(r,r)
         alpha=rfrf/rr

! calcola p'=r'+alpha*p
         do ivol=1,nvol_eo
            do icol=1,ncol
               p(icol,ivol)=r(icol,ivol)+alpha*p(icol,ivol)
            enddo
         enddo

! calcola alphas=alpha*zfs*betas/zas*beta
         do iterm=1,nterm
         if(flag(iterm).eq.1) then
            alphas(iterm)=alpha*zfs(iterm)*betas(iterm)/
     $ (zas(iterm)*betaa)
! calcola ps'=r'+alpha*p
            do icol=1,ncol
               do ivol=1,nvol_eo
                  ps(icol,ivol,iterm)=zfs(iterm)*r(icol,ivol)+
     $ alphas(iterm)*ps(icol,ivol,iterm)
               enddo
            enddo

! check sui residui
            if(rr*zfs(iterm).lt.residue) then
               flag(iterm)=0
            endif

! passa tutti gli f in quelli attuali
            zps(iterm)=zas(iterm)
            zas(iterm)=zfs(iterm)
         endif
         enddo
         rr=rfrf

! check sul residuo
         if(rfrf<residue) exit
      enddo
# 334 "inverter.f"
      return
      end


c=========================================================
      subroutine multi_shift_fuffa_inverter(x,b,potc,pole,nterm,ud)
c written by Sanfo
c version 1.0 - 2007/2008
c esegue lo shift termine a termine e lo schiaffain x:
c x_i = (pole_i+MM~)^-1 b
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 347 "inverter.f" 2

c argomenti
      complex x(ncol,nvol_eo,nmr)
      complex b(ncol,nvol_eo)
      complex potc
      real pole(nmr)
      integer nterm
      integer ud

c variabili interne
      integer iterm
# 366 "inverter.f"
      ud=ud !to avoid warning

      do iterm=1,nterm
         call singol_inverter(x(1,1,iterm),b,potc,pole(iterm))
      enddo




      return
      end

c=========================================================
      subroutine multi_shift_summed_inverter
     $ (x,b,potc,cost,pole,coef,nterm,ud)
c written by Sanfo
c version 1.0 - 2008
c somma i vari termini calcolati sopra con i loro coefficienti
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 387 "inverter.f" 2

c argomenti
      complex x(ncol,nvol_eo),b(ncol,nvol_eo)
      complex potc
      real cost
      real pole(nmr),coef(nmr)
      integer nterm
      integer ud

c variabili interne
      complex s(ncol,nvol_eo,nmr)
      integer ivol,icol,iterm





      call multi_shift_inverter(s,b,potc,pole,nterm,ud)

      do ivol=1,nvol_eo
         do icol=1,ncol
            x(icol,ivol)=cost*b(icol,ivol)
            do iterm=1,nterm
               x(icol,ivol)=x(icol,ivol)+coef(iterm)*s(icol,ivol,iterm)
            enddo
         enddo
      enddo

      return
      end
# 107 "main.f" 2
# 1 "gauss_su3.f" 1
c =========================================================================
      subroutine gauss_matrix(h,sigma)
c =========================================================================
c IT'S GOOD ONLY FOR SU(3) !!!
!!===========================================================================!!
!! A gauss_matrix is a linear combination of eight LAMBDA matrices. The !!
!! coefficients are real numbers with gaussian distribution. We have !!
!! f(z) = 1/(sigma pi) * exp[ - z z* / sigma] , z complex, !!
!! so that the integral f(z)dz dz* is one. !!
!! Further we have < z z* > = sigma . !!
!! !!
!! We first create 8 gaussian numbers !!
!! with sigma=1 and then we use the eight real components of this vector !!
!! to get a linear combination of the eight LAMDA matrices. Note, that the !!
!! expectation value of each coefficient is 1/2 because we created complex !!
!! numbers with expectation value 1. !!
!! !!
!! Our choice of the LAMBDA matrices is so, that the trace of the product !!
!! of two LAMBDAS A and B is 2 if A = B and 0 else: !!
!! !!
!! | 0 1 0 | | 0 -i 0 | | 1 0 0 | !!
!! lambda1 = | 1 0 0 | , lambda2 = | i 0 0 | , lambda3 = | 0 -1 0 | !!
!! | 0 0 0 | | 0 0 0 | | 0 0 0 | !!
!! !!
!! | 0 0 1 | | 0 0 -i | | 0 0 0 | !!
!! lambda4 = | 0 0 0 | , lambda5 = | 0 0 0 | , lambda6 = | 0 0 1 | !!
!! | 1 0 0 | | i 0 0 | | 0 1 0 | !!
!! !!
!! | 0 0 0 | | 1 0 0 | !!
!! lambda7 = | 0 0 -i | , lambda8 = 1/sqrt(3) | 0 1 0 | . !!
!! | 0 i 0 | | 0 0 -2 | !!
!! !!
!!===========================================================================!!

      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 38 "gauss_su3.f" 2
      parameter(twopi = 6.28318530718)
      parameter(one_by_sqrt3 = 0.57735026919)
      parameter(two_by_sqrt3 = 1.154700538379)
      complex h(ncol,ncol)
      real temp(ncol2),phi,radius,raux
      real sigma

      do igen = 1,ncol2/2
         phi = twopi * ran2()
         raux = - sigma * log(ran2())
         radius = sqrt(raux)
         temp(2*igen-1) = radius * cos(phi)
         temp(2*igen) = radius * sin(phi)
      enddo

      x1 = temp(3) + one_by_sqrt3 * temp(8)
      x2 = 0.0
      h(1,1) = cmplx(x1,x2)

      h(1,2) = cmplx(temp(1),-temp(2))

      h(1,3) = cmplx(temp(4),-temp(5))

      h(2,1) = cmplx(temp(1),temp(2))

      x1 = - temp(3) + one_by_sqrt3 * temp(8)
      x2 = 0.0
      h(2,2) = cmplx(x1,x2)

      h(2,3) = cmplx(temp(6),-temp(7))

      h(3,1) = cmplx(temp(4),temp(5))

      h(3,2) = cmplx(temp(6),temp(7))

      x1 = - two_by_sqrt3 * temp(8)
      x2 = 0.0
      h(3,3) = cmplx(x1,x2)

      return
      end

c =========================================================================
      subroutine gauss_vector(v,sigma)
c =========================================================================
!!===========================================================================!!
!! A gauss vector has complex components z which are gaussian distributed !!
!! with <z~ z> = sigma !!
!!===========================================================================!!


      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 92 "gauss_su3.f" 2
      parameter(twopi = 6.28318530718)
c parameter(sigma = 1.0000000)
      complex v(ncol)
      real temp1,temp2,phi,radius,raux
      real sigma

      do icol = 1,ncol
         phi = twopi * ran2()
         raux = - sigma * log(ran2())
         radius = sqrt(raux)
         temp1 = radius * cos(phi)
         temp2 = radius * sin(phi)
         v(icol) = cmplx(temp1,temp2)
      enddo

      return
      end
# 108 "main.f" 2
# 1 "eigenvalue.f" 1
c=========================================================
      subroutine eigenval(autmin,autmax,potc,ud)
c written by Sanfo
c find minimum and maximum eigenvalue of even\even part of M^+M
c by using an approximated (and probably not stable) algorithm
c This is used to find the appropriate interval where to scale the
c rational expansion
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 11 "eigenvalue.f" 2

c arguments
      real autmin,autmax
      complex potc
      integer ud

c common blocks
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue

c internal variables
      integer aut
      complex y(ncol,nvol_eo),x(ncol,nvol_eo)
      real norm
      integer icol,ivol,iter
      integer niter

c parameters
      real sigma
      parameter(sigma=1)





      ud=ud !to avoid warning

      do aut=1,2
         do ivol=1,nvol_eo
            call gauss_vector(x(1,ivol),sigma)
         enddo

         if(aut.eq.1) then
            niter=50
         else
            niter=5
         endif

         do iter=1,niter

            call m2d(y,x,potc,mass)

            if(aut.eq.2) then
               do ivol=1,nvol_eo
                  do icol=1,ncol
                     y(icol,ivol)=autmax*x(icol,ivol)-y(icol,ivol)
                  enddo
               enddo
            endif

            norm=0
            do ivol=1,nvol_eo
               do icol=1,ncol
                  norm=norm+real(y(icol,ivol))**2+aimag(y(icol,ivol))**2
               enddo
            enddo
            norm=(norm/nvol_eo/ncol)**0.5
            do ivol=1,nvol_eo
               do icol=1,ncol
                  x(icol,ivol)=y(icol,ivol)/norm
               enddo
            enddo





         enddo

         if(aut.eq.1) then
            autmax=norm
         else
            autmin=autmax-norm
         endif

      enddo






      return
      end


c=========================================================
      subroutine scale_rhmc
c written by Sanfo
c version 1.0 - 05/08
c riscala i parametri dell'espansione adattandoli
c all'intervallo opportuno.
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 106 "eigenvalue.f" 2

c common blocks
      complex potc1,potc2
      common/param3/potc1,potc2

c AZIONE
c approssimazione di (m+x)^(-nf/4)
c deve essere usato per il calcolo dell'azione, che avviene qui,
c per il calcolo di Polyakov e quant'altro, che avviene in CHIRAL
c e per le equazioni della dinamica, che avvengono in compute_ipdot
      integer nterm_a
      real cost_a,pole_a,coef_a
      common/rhmc_a/nterm_a,cost_a,pole_a(nmr),coef_a(nmr)

c HEAT-BATH
c approssimazione di (m+x)^(nf/8)
c viene usato per la generazione dei campi phi tramite heat-bath
c nella funzione: create_phi
      integer nterm_h
      real cost_h,pole_h,coef_h
      common/rhmc_h/nterm_h,cost_h,pole_h(nmr),coef_h(nmr)

c ESPANSIONI ORIGINALI DELL'RHMC
c sono quelle che poi vengono scalate
      real scala,min_rhmc,max_rhmc
      real ocost_a,opole_a,ocoef_a,z_a
      real ocost_h,opole_h,ocoef_h,z_h
      common/rhmc_o/scala,min_rhmc,max_rhmc,
     $ ocost_a,opole_a(nmr),ocoef_a(nmr),z_a,
     $ ocost_h,opole_h(nmr),ocoef_h(nmr),z_h

c variabili interne
      real autmin1,autmax1
      real autmin2,autmax2
      real autmin,autmax
      real scala_az,scala_hz
      integer iterm

c calcola gli autovalori massimo e minimo per le due matrici separatamente
      call eigenval(autmin1,autmax1,potc1,1)
      call eigenval(autmin2,autmax2,potc2,2)

c trova il minimo ed il massimo "assoluti"
      if(autmax1.gt.autmax2) then
         autmax=autmax1
      else
         autmax=autmax2
      endif
      if(autmin1.lt.autmin2) then
         autmin=autmin1
      else
         autmin=autmin2
      endif




c setta i fattori di scala
      scala=autmax*1.1
      scala_hz=scala**z_h
      scala_az=scala**z_a

c check: stan nell'intervallo?
      if(((autmin/scala).lt.min_rhmc).or.((autmax/scala).gt.max_rhmc))
     $ then
         write(*,*) "ATTENTION! Approximation not valid!!!"
         write(*,*) "Interval of validity of approximation:",
     $ min_rhmc,max_rhmc
         write(*,*) "Interval of eigenvalue:",
     $ autmin,autmax
      endif

c scala il tutto
      cost_h=ocost_h*scala_hz
      do iterm=1,nterm_h
         pole_h(iterm)=opole_h(iterm)*scala
         coef_h(iterm)=ocoef_h(iterm)*scala*scala_hz
      enddo
      cost_a=ocost_a*scala_az
      do iterm=1,nterm_a
         pole_a(iterm)=opole_a(iterm)*scala
         coef_a(iterm)=ocoef_a(iterm)*scala*scala_az
      enddo

      return
      end
# 109 "main.f" 2
# 1 "test.f" 1
c=========================================================
      subroutine test_reversibility(scale_each_step)
c written by Sanfo
c testa la reversibilit delle equazioni usate per la dinamica
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 8 "test.f" 2

c argomenti
      integer scale_each_step

c variabili passate da blocchi common
      complex u
      common/field/u(ncol,ncol,4,nvol)
      complex p_w,ipdot
      common/momenta/p_w(ncol,ncol,4,nvol),ipdot(ncol,ncol,4,nvol)

c variabili interne
      real dif_act01,dif_act02
      real normp_w,normu
      integer nudof,npdof,nfdof,nhdof
      parameter(nudof=4*8*nvol,npdof=nudof,nfdof=3*nvol_eo,
     $ nhdof=nudof+npdof+nfdof)
      complex ur(ncol,ncol,4,nvol)
      complex p_wr(ncol,ncol,4,nvol)

c VARIABILI CHE CONTENGONO L'AZIONE
      real az_0(nvol),az_1(nvol),az_2(nvol)

c crea la configurazione iniziale
      call create_momenta
      call create_phi
      call action(az_0)

c copia la configurazione iniziale
      call copyconf(p_wr,p_w)
      call copyconf(ur,u)

c va avanti,calcola l'azione,torna indietro e ricalcola l'azione
      call dinamica(scale_each_step)
      call action(az_1)

      call revert_momenta

      call dinamica(scale_each_step)
      call revert_momenta
      call action(az_2)

c calcola la differenza tra le norme dei campi u e p_w e delle azioni
      call diff_norm(normp_w,p_wr,p_w)
      call diff_norm(normu,ur,u)
      call dif_action(dif_act01,az_0,az_1)
      call dif_action(dif_act02,az_0,az_2)

      write(*,*) "Azione1-azione0:",dif_act01
      write(*,*) "Azione2-azione0:",dif_act02

      write(*,*) "|P''-P|/npdof^0.5:",normp_w/real(npdof)**0.5
      write(*,*) "|U''-U|/nudof^0.5:",normu/real(nudof)**0.5
      write(*,*) "|H''-H|/nhdof^0.5:",dif_act02/real(nhdof)**0.5
      write(*,*) "|H''-H|/|H'-H|:",dif_act02/dif_act01

      return
      end



c=========================================================
      subroutine revert_momenta()
c written by Sanfo
c inverte il segno agli H in modo da tornare indietro
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 76 "test.f" 2

c variabili passate con blocchi common
      complex p_w,ipdot
      common/momenta/p_w(ncol,ncol,4,nvol),ipdot(ncol,ncol,4,nvol)

      do ivol=1,nvol
         do idir=1,4
            do icol1=1,ncol
               do icol2=1,ncol
                  p_w(icol1,icol2,idir,ivol)=-p_w(icol1,icol2,idir,ivol)
               enddo
            enddo
         enddo
      enddo

      return
      end


c=========================================================
      subroutine diff_norm(norm,mr,m)
c written by Sanfo
c calcola una fuffa norma della differenza tra il campo ottenuto
c ritornando indietro e quello originale per la coppia di matrici
c passate da argomenti
c=========================================================

      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 106 "test.f" 2

c argomenti
      real norm
      complex mr(ncol,ncol,4,nvol)
      complex m(ncol,ncol,4,nvol)

c variabili interne
      complex c1


      norm=0

      do ivol=1,nvol
         do idir=1,4
            do icol1=1,ncol
               do icol2=1,ncol
                  c1=mr(icol1,icol2,idir,ivol)-m(icol1,icol2,idir,ivol)
                  norm=norm+real(c1*conjg(c1))
               enddo
            enddo
         enddo
      enddo

      norm=norm**0.5

      return
      end


c=========================================================
      subroutine back_campo(mr,m)
c written by Sanfo
c copia il campo m in mr
c=========================================================

      implicit real (a-h,o-z)
      implicit integer (i-n)
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 144 "test.f" 2

c argomenti
      complex mr(ncol,ncol,4,nvol)
      complex m(ncol,ncol,4,nvol)


      do ivol=1,nvol
         do idir=1,4
            do icol1=1,ncol
               do icol2=1,ncol
                  mr(icol1,icol2,idir,ivol)=m(icol1,icol2,idir,ivol)
               enddo
            enddo
         enddo
      enddo

      return
      end



c=========================================================
      subroutine test_invertibility()
c written by Sanfo
c testa la reversibilit delle equazioni usate per la dinamica
c=========================================================
      implicit none
# 1 "parameters.f" 1
c include fixed defintions
# 1 "definitions.f" 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c internal fixed definitions c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer u_control,u_allright,u_lattice,u_random,u_rhmc,u_meas
     $ ,u_magnet
      integer u_ferm_11,u_ferm_12,u_ferm_13,u_ferm_14 ,u_ferm_15
     $ ,u_ferm_16,u_ferm_17,u_ferm_21,u_ferm_22,u_ferm_23,u_ferm_24
     $ ,u_ferm_25,u_ferm_26,u_ferm_27

      parameter(u_ferm_11=11,u_ferm_21=21)
      parameter(u_ferm_12=12,u_ferm_22=22)
      parameter(u_ferm_13=13,u_ferm_23=23)
      parameter(u_ferm_14=14,u_ferm_24=24)
      parameter(u_ferm_15=15,u_ferm_25=25)
      parameter(u_ferm_16=16,u_ferm_26=26)
      parameter(u_ferm_17=17,u_ferm_27=27)
      parameter(u_control=30,u_allright=31,u_lattice=32,u_random=33
     $ ,u_rhmc=34,u_meas=35,u_magnet=36)

c important numbers
      real pigr,my_phone_in_la_sapienza
      parameter(pigr=3.141592654,my_phone_in_la_sapienza=0649913487)

c EVEN-ODD or ODD-EVEN (used to apply Dirac Matrix on a vector)
      integer EO,OE,EO_DAG,OE_DAG,EO_P_OE_DAG,OE_P_EO_DAG
      parameter(EO=0,OE=1,EO_DAG=2,OE_DAG=3,EO_P_OE_DAG=4,OE_P_EO_DAG=5)

c parameter for various part of the code
      integer n_quark_for_flavour
      integer unit_each,save_each,tune_nstep_each
      integer nmr
      parameter(tune_nstep_each=25) !number of traj for tuning nstep
      parameter(n_quark_for_flavour=1) !number of quark for flavour - keep 1
      parameter(unit_each=5) !number of MD micro_step between reunitarization
      parameter(save_each=5000) !number of trajectory between two saved ones
      parameter(nmr=15) !maximum oredr for rational expansion

c definitions of algorithm selection parameter (select from parameters.f)
      integer d_leapfrog,d_omelyan
      integer i_singol,i_multi
      parameter(d_leapfrog=0,d_omelyan=1)
      parameter(i_singol=0,i_multi=1)
# 3 "parameters.f" 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c tunable parameters c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set the quantization condition for the magnetic field
c if relaxed=0 the field is proportional to 1/L
c if relaxed=1 the field is proportional to 1/L^2
      integer relaxed
      parameter(relaxed=1)

c lattice geometry
      integer nx,ny,nz,nvol,nvolh,nvol_eo
# 1 "dimensioni.f" 1
       integer lato,nt
       parameter(lato=2,nt=4)
# 17 "parameters.f" 2
      parameter(nx=lato,ny=lato,nz=lato) !lattice geometry
      parameter(nvol=nx*ny*nz*nt,nvolh=nvol/2) !lattice site and e.site number



      parameter(nvol_eo=nvolh)

c color definition
      integer ncol,ncolm1,ncol2,nperm
      parameter(ncol=3,ncolm1=ncol-1,ncol2 =ncol*ncol)
      parameter(nperm=ncol**(ncol-3)+10000)

c algorithm choice
      integer alg_din
      integer alg_inv
      parameter(alg_din=d_omelyan)!algorithm for MD: d_omelyan or d_leapfrog
      parameter(alg_inv=i_multi) !algorithm for inverter: i_multi or i_singol
# 172 "test.f" 2

c variabili passate da blocchi common
      complex u
      common/field/u(ncol,ncol,4,nvol)
      complex p_w,ipdot
      common/momenta/p_w(ncol,ncol,4,nvol),ipdot(ncol,ncol,4,nvol)
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue
      complex potc1,potc2
      common/param3/potc1,potc2
      integer nterm
      real cost,pole,coef
      common/rhmc_a/nterm,cost,pole(nmr),coef(nmr)

c variabili interne
      complex temp(ncol,nvol_eo)
      complex origi(ncol,nvol_eo)
      complex final(ncol,nvol_eo)
      integer i,icol,ivol
      complex c1
      real norm

c parametri
      real sigma
      parameter(sigma=1)

      call scale_rhmc

c estrae un vettore a caso
      do ivol=1,nvol_eo
         call gauss_vector(temp(1,ivol),sigma)
      enddo

c copia il vettore estratto in quello originale
      do icol=1,ncol
         do ivol=1,nvol_eo
            origi(icol,ivol)=temp(icol,ivol)
         enddo
      enddo

c applica l'inverter 4 volte, così ottiene temp=m^-1orig
      do i=1,4
         call multi_shift_summed_inverter
     $ (final,temp,
     $ potc1,cost,pole,coef,
     $ nterm,1)

         do icol=1,ncol
            do ivol=1,nvol_eo
               temp(icol,ivol)=final(icol,ivol)
            enddo
         enddo
      enddo

c riapplica M così in teoria avrei quello originale



      call m2d(final,temp,potc1,mass)



c stampa il vettore originale, quello finale e quello intermedio
c intanto calcola la norma della differenza tra l'originale ed il finale
      norm=0
      do icol=1,ncol
         do ivol=1,nvol_eo
c write(*,*) real(origi(icol,ivol)),real(final(icol,ivol)),
c $ real(temp(icol,ivol))
            c1=origi(icol,ivol)-final(icol,ivol)
            norm=norm+real(c1*conjg(c1))
         enddo
      enddo
      norm=norm**0.5
      norm=norm/nvol_eo/ncol
      write(*,*) "norma della differenza:",norm

      return
      end
# 110 "main.f" 2
