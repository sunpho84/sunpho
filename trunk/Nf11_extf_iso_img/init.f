c=========================================================
      subroutine load_par(init_flag,termalizza,time_to_run,n_rand
     $     ,immu_quark,immu_iso,remu,file_rhmc,val_ext_f)
c  written by Sanfo
c  load (by standard input) simulation parameters
c=========================================================
      implicit none
#include "parameters.f"

c     arguments
      integer init_flag
      integer termalizza
      integer time_to_run
      integer n_rand
      real immu_quark,immu_iso,remu
      character(LEN=30) file_rhmc
      real val_ext_f

c     common-block passed variables
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue
      real dt_md
      integer nstep_md
      complex ieps,iepsq,ieps2q,ieps3q
      common/param2/dt_md,nstep_md,ieps,iepsq,ieps2q,ieps3q

c     internal variables
      real dt
      character(LEN=10) cosa

!     0 -> cold; 1 -> hot; 2|3 -> stored
      read(*,'(a3,i10)',ERR=201,END=201) cosa,init_flag 
      if(cosa.eq."cho") goto 102
 201  write(*,*) "Errore riga 1 (cho)"
      stop

!     0 -> no, 1 -> yes
 102  read(*,'(a4,i10)',ERR=202,END=202) cosa,termalizza 
      if(cosa.eq."term") goto 103
 202  write(*,*) "Errore riga 2 (term)"
      stop

!     coupling in the fundamental
 103  read(*,'(a4,f10.20)',ERR=203,END=203) cosa,beta
      if(cosa.eq."beta") goto 104
 203  write(*,*) "Errore riga 3 (beta)"
      stop

!     fermion mass
 104  read(*,'(a4,f10.20)',ERR=204,END=204) cosa,mass
      if(cosa.eq."mass") goto 105
 204  write(*,*) "Errore riga 4 (mass)"
      stop

!     total running time
 105  read(*,*,ERR=205,END=205) cosa,time_to_run
      if(cosa.eq."time") goto 106
 205  write(*,*) "Errore riga 5 (time)"
      stop

!     number of md steps per trajectory
 106  read(*,'(a5,i10)',ERR=206,END=206) cosa,nstep_md
      if(cosa.eq."nstep") goto 107
 206  write(*,*) "Errore riga 6 (nstep)"
      stop

!     trajectory length
 107  read(*,'(a2,f10.20)',ERR=207,END=207) cosa,dt
      if(cosa.eq."dt") goto 108
 207  write(*,*) "Errore riga 7 (dt)"
      stop

!     stopping criterion for inverter
 108  read(*,'(a7,f10.20)',ERR=208,END=208) cosa,residue
      if(cosa.eq."residue") goto 109
 208  write(*,*) "Errore riga 8 (residue)"
      stop

!     number of random vectors 
 109  read(*,'(a6,i10)',ERR=209,END=209) cosa,n_rand
      if(cosa.eq."n_rand") goto 110
 209  write(*,*) "Errore riga 9 (n_rand)"
      stop

!     imaginary isospin chemical potential/(pi T) 
 110  read(*,'(a8,f10.20)',ERR=210,END=210) cosa,immu_iso
      if(cosa.eq."immu_iso") goto 111
 210  write(*,*) "Errore riga 10 (immu_iso)"
      stop

!     imaginary quark chemical potential/(pi T) 
 111  read(*,'(a10,f10.20)',ERR=211,END=211) cosa,immu_quark
      if(cosa.eq."immu_quark") goto 112
      stop
 211  write(*,*) "Errore riga 11 (immu_quark)"
      stop

!     finite real isospin chemical potential/(pi T)
 112  read(*,'(a8,f10.20)',ERR=212,END=212) cosa,remu
      write(*,*) cosa
      if(cosa.eq."remu_iso") goto 113
 212  write(*,*) "Errore riga 12 (remu_iso)"
      stop

!     external field value
 113  read(*,*,ERR=213,END=213) cosa,val_ext_f
c 113  read(*,'(a,f10.10)',ERR=213,END=213) cosa,val_ext_f
      if(cosa.eq."extf") goto 114
 213  write(*,*) "Errore riga 13 (extf)"
      stop

 114  continue
!! file for rhmc expansion
#ifdef ficp
      file_rhmc="rhmc8"
#else
      file_rhmc="rhmc4"
#endif

      mass2=mass*mass
      dt_md=dt/nstep_md

#ifdef chiral_meas
      if(n_rand.eq.0) then
         write(*,*) "Can't set random vector number to 0"
         n_rand=1
      endif
#endif

      write(*,*)
      write(*,*) " ----Simulation parameters----"
      if(init_flag==0) then
         write(*,*) " Cold start",init_flag
      endif
      if(init_flag==1) then
         write(*,*) " Hot start",init_flag
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
#ifdef chiral_meas
      write(*,*) " Random vector number: ",n_rand
#else
      write(*,*) " Chiral measurement switched off at compilation time"
#endif

#ifdef isotropic
      write(*,*) " Isotropic simulation on"
#endif

      write(*,*) " Isospin potential (in PI unit): ",immu_iso
      write(*,*) " Quark potential (in PI unit): ",immu_quark
#ifdef ficp
      write(*,*) " Real isospin potential (in PI unit): ",remu
#else
      write(*,*)
     $     " Real isospin potential switched off at compilation time"
      if(remu.ne.0) then
         write(*,*) "Error: real isospin potential different from 0!"
         stop
      endif
#endif
#ifdef mf
      if(relaxed.eq.0) then
         write(*,*) " External field: ",val_ext_f," (1/L unit)"
      elseif(relaxed.eq.1) then
         write(*,*) " External field: ",val_ext_f," (1/L^2 unit)"
      else
         write(*,*) " Unkwnown quantization format!"
      endif
#else
      write(*,*) " Magnetic field switched off at compilation time"
      if(val_ext_f.ne.0) then
         write(*,*) "Error: external field different from 0!"
         stop
      endif
#endif
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
c  written by Sanfo
c  load rhmc parameters
c=========================================================
      implicit none
#include "parameters.f"

c     arguments
      character(LEN=30) file_rhmc

c     common-block passed variables
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
     $     ocost_a,opole_a(nmr),ocoef_a(nmr),z_a,
     $     ocost_h,opole_h(nmr),ocoef_h(nmr),z_h


c     internal variables
      integer i


cc-----------------------------
c     LOAD RHMC PARAMETER 
cc-----------------------------
      open(u_rhmc,file=file_rhmc,status='old')

!     load extremes of validity of rhmc approximation
      read(u_rhmc,*) min_rhmc,max_rhmc

!     load action part parameter
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
      
!     load heat bath part
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

#ifdef debug1
      call print_rhmc
#endif
      
      if((nterm_a.gt.nmr).or.(nterm_h.gt.nmr)) then
         write(*,*) "Attention, the order requested in the file"
         write(*,*) "is greater than the one defined in parameters.f"
         stop
      endif

      return
      end

c=========================================================
      subroutine print_rhmc()
c  written by Sanfo
c  print rhmc parameters
c=========================================================
      implicit none
#include "parameters.f"

c     common-block passed variables
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
     $     ocost_a,opole_a(nmr),ocoef_a(nmr),z_a,
     $     ocost_h,opole_h(nmr),ocoef_h(nmr),z_h

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
c  written by Sanfo
c  set various parameter regarding MD and chemical potentials
c=========================================================
      implicit none
#include "parameters.f"

c     arguments
      real immu_quark,immu_iso,remu

c     common-block passed variables
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

c     internal variables
      integer i
      real immu1,immu2,remu_fis

cc-----------------------------
cc compute  ncol!
      ncfact=1
      do i=1,ncol
         ncfact=ncfact*i
      enddo
cc-----------------------------


cc-----------------------------
cc get back to single species imaginary chemical potentials
cc  multiply by pi and divide by nt      
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
c  written by Sanfo
c  initialize the simulation
c=========================================================
      implicit none
#include "parameters.f"

c     arguments
      integer time_to_run
      integer n_rand
      integer n_traj

c     common-block passed variables
      complex u,usave
      common/field/u(ncol,ncol,4,nvol)
      common/field1/usave(ncol,ncol,4,nvol)

c     internal variables
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
c     READ INPUT & RHMC
cc-----------------------------

      call load_par(init_flag,termalizza,time_to_run,n_rand,immu_quark,
     $     immu_iso,remu,file_rhmc,val_ext_f)
      call load_rhmc(file_rhmc)
      call set_const(immu_quark,immu_iso,remu)

c-------------------------------------------

c-------------------------------------------
c     INITIALIZATION OPERATION
c-------------------------------------------     

      call ranstart                 !! initialize random number generator
      call generate_permutations    !! table of ncol permutations
      call geometry                 !! set up lattice geometry
      call initialize_lattice(init_flag,n_traj)  !! initialize lattice
#ifdef mf
      call initialize_extf(val_ext_f) !! initialize external field
#endif
      call addrem_stagphase
      call copyconf(usave,u) !! save configuration at first

c-------------------------------------------
c     Observables file. 
c
c     *Explenation of content of observable files*
c     -------------------------------------------------
c     
c     In each file the first colomn is the trajectory number.
c
c     -meas_out:
c      2 acceptance of actual MD trajectory (0=discarded, 1 accepted)
c      3 spatial plaquette
c      4 temporal plaquette
c      5 real part of polyakov loop
c      6 imaginary part of polyakov loop
c
c     -magnet_out:
c      2 real part of magnetic susceptibility
c      3 imag part
c
c     Each "ferm" file, is a different observable, calculated via an
c     approximated noisy estimator.
c
c     In each ferm file the four columns are arranged as follow:
c      2 real part mean of the estimates over current config
c      3 imag part
c      4 squared mean of real part
c      5 squared mean of imag part
c
c     The fist figure in the filename is relative to the quark (1=up, 2=down)
c
c     The second figure is the observable, as follow:
c      1 chiral condensate
c      2 energy density
c      3 density of quark
c      4 pressure density
c      5 6 7 electic current
c      


      open(u_meas,file='meas_out',status='unknown')

#ifdef chiral_meas
      open(u_ferm_11,file='ferm11_out',status='unknown')
      open(u_ferm_21,file='ferm21_out',status='unknown')
      open(u_ferm_22,file='ferm22_out',status='unknown')
      open(u_ferm_12,file='ferm12_out',status='unknown')
      open(u_ferm_13,file='ferm13_out',status='unknown')
      open(u_ferm_23,file='ferm23_out',status='unknown')
      open(u_ferm_14,file='ferm14_out',status='unknown')
      open(u_ferm_24,file='ferm24_out',status='unknown')
#endif

#ifdef mf
      open(u_ferm_15,file='ferm15_out',status='unknown')
      open(u_ferm_25,file='ferm25_out',status='unknown')
      open(u_ferm_16,file='ferm16_out',status='unknown')
      open(u_ferm_26,file='ferm26_out',status='unknown')
      open(u_ferm_17,file='ferm17_out',status='unknown')
      open(u_ferm_27,file='ferm27_out',status='unknown')
      open(u_magnet,file='magnet_out',status='unknown')
#endif
      return
      end
