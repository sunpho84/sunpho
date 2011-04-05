c=========================================================
      subroutine load_par(init_flag,termalizza,time_to_run,n_rand
     $     ,remu_iso,file_rhmc)
c  written by Sanfo
c  load (by standard input) simulation parameters
c=========================================================
      implicit none
      include "parameters.f"

c     arguments
      integer init_flag
      integer termalizza
      integer time_to_run
      integer n_rand
      real remu_iso
      character(LEN=30) file_rhmc

c     common blocks
      integer nran,iad_flag
      real beta_f,beta_a,rho,pi
      common/param/nran,beta_f,beta_a,rho,pi,iad_flag
      real mass,mass2,residue
      common/param2/mass,mass2,residue
      real dt_md
      integer nstep_md
      complex ieps,iepsq,ieps2q,ieps3q
      common/param3/dt_md,nstep_md,ieps,iepsq,ieps2q,ieps3q

c     internal variables
      real dt


      beta_a=0
      
      read(*,*) init_flag       !! 0 -> cold; 1 -> hot; 2|3 -> stored
      read(*,*) termalizza      !! 0 -> no, 1 -> yes
      read(*,*) beta_f          !! coupling in the fundamental
      read(*,*) mass            !! fermion mass
      read(*,*) time_to_run     !! total running time
      read(*,*) nstep_md        !! number of md steps per trajectory
      read(*,*) dt              !! trajectory length
      read(*,*) residue         !! stopping criterion for inverter
      read(*,*) n_rand          !! number of random vectors 
      remu_iso=0        !! imaginary isospin chemical potential/(pi T) 
      read(*,*) file_rhmc       !! file for rhmc expansion

      mass2=mass*mass
      dt_md=dt/nstep_md
      if(n_rand.eq.0) then
         write(*,*) "Can't set random vector number to 0"
         n_rand=1
      endif

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
      write(*,*) " Coupling constant: ",beta_f
      write(*,*) " Mass, squared mass: ",mass,mass2
      write(*,*) " Will run for at least: ",time_to_run," seconds"
      write(*,*) " MD micro-step for each trajectory: ",nstep_md
      write(*,*) " Trajectories length: ",dt
      write(*,*) " MD micro-step length: ",dt_md
      write(*,*) " Residual for inverter: ",residue
      write(*,*) " Random vector number: ",n_rand
      write(*,*) " Isospin potential (in PI unit): ",remu_iso
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

      return
      end



c=========================================================
      subroutine load_rhmc(file_rhmc)
c  written by Sanfo
c  load rhmc parameters
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"

c     arguments
      character(LEN=30) file_rhmc

c     common blocks
      real mass,mass2,residue
      common/param2/mass,mass2,residue

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

      if(debug.ge.1) then
         call print_rhmc
      endif
      
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
      include "parameters.f"

c     common blocks
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


      write(*,*)
      write(*,*) "----Calculation of action and dynamic(-nf/2)----"
      write(*,*) " Number of terms:",nterm_a
      write(*,*) " Degree of expansion:",z_a
      write(*,*) " Constant term:",cost_a
      do i = 1,nterm_a
         write(*,*) " n°, pole, coef: ",i,pole_a(i),coef_a(i)
      enddo
      write(*,*)
      write(*,*) "----Heat Bath(nf/4)----"
      write(*,*) " Number of terms:",nterm_h
      write(*,*) " Degree of expansion:",z_h
      write(*,*) " Constant term:",cost_h
      do i = 1,nterm_h
         write(*,*) " n°, pole, coef: ",i,pole_h(i),coef_h(i)
      enddo
      write(*,*)


      return
      end

c=========================================================
      subroutine set_const(remu_iso)
c  written by Sanfo
c  set various parameter regarding MD and chemical potentials
c=========================================================
      implicit none
      include "parameters.f"

c     arguments
      real remu_iso
      
c     common-block passed variables
      integer nran,iad_flag
      real beta_f,beta_a,rho,pi
      common/param/nran,beta_f,beta_a,rho,pi,iad_flag
      real mass,mass2,residue
      common/param2/mass,mass2,residue
      integer ncfact,perm,sign
      common/perm/ncfact,perm(nperm,ncol),sign(nperm)
      real remu
      real emu,ememu
      common/param4/remu,emu,ememu
      real dt_md
      integer nstep_md
      complex ieps,iepsq,ieps2q,ieps3q
      common/param3/dt_md,nstep_md,ieps,iepsq,ieps2q,ieps3q

c     internal variables
      integer i


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
      remu=pigr*remu_iso/float(nt)
cc-----------------------------

cc-----------------------------
cc set dt and sub-multiples
      iepsq=cmplx(0.0,dt_md/4.0)
      ieps2q=2.0*iepsq
      ieps3q=3.0*iepsq
      ieps=4.0*iepsq
      
      emu  = exp(remu)
      ememu = exp(-remu)
      
cc-----------------------------
      
      
cc-----------------------------
cc print parameters
      write(*,*)
      write(*,*) "----Chemical potential in 'physical' unit----"
      write(*,*) " Chemical potential for U quark: ", emu
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
      include "parameters.f"

c     arguments
      integer time_to_run
      integer n_rand
      integer n_traj

c     common blocks
      complex u,usave
      common/field/u(ncol,ncol,4,nvol)
      common/field1/usave(ncol,ncol,4,nvol)

c     internal variables
      integer init_flag
      integer termalizza
      real remu_iso
      character(LEN=30) file_rhmc
      integer f
      integer indicefile(5)
      character(LEN=15) nomefile(5)
      data indicefile /u_meas,u_ferm_1,u_ferm_2,u_ferm_3,u_ferm_4/
      data nomefile /"meas_out","ferm1_out","ferm2_out","ferm3_out"
     $     ,"ferm4_out"/

      write(*,*)
      write(*,*) "----System parameter----"
      write(*,*) " Number of colors: ", ncol
      write(*,*) " Lattice size: ",nx,ny,nz,nt
      write(*,*)
      
cc-----------------------------
c     READ INPUT & RHMC
cc-----------------------------

      call load_par(init_flag,termalizza,time_to_run,n_rand,remu_iso
     $     ,file_rhmc)
      call load_rhmc(file_rhmc)
      call set_const(remu_iso)

c-------------------------------------------

c-------------------------------------------
c     INITIALIZATION OPERATION
c-------------------------------------------     

      call ranstart                 !! initialize random number generator
      call generate_permutations    !! table of ncol permutations
      call generate_su2rhotables    !! table of random su2 for metro and hb
      call geometry                 !! set up lattice geometry
      call initialize_lattice(init_flag,n_traj)  !! initialize lattice
      call addrem_stagphase
      call copyconf(usave,u) !! save configuration at first

c     -------------------------------------------
c     Observables file. 
c
c     *Explenation of content of observable files*
c     --------------------------------------------
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
c     The fist figure in the filename is relative to the quark (1=up)
c
c     The second figure is the observable, as follow:
c      1 chiral condensate
c      2 energy density
c      3 density of quark
c      4 pressure density
c      

      do f=1,5
         open(indicefile(f),file=nomefile(f),status='unknown')
      enddo


      return
      end
