CC=====================================================
      program su_n_quark_chem_pot
CC
CC  Program for simulating SU(N_c) gauge theories
CC  with 2 staggered flavors in presence of equal 
CC  imaginary  chemical potential coupled to each
CC  quark
CC
CC  includes fundamental [and adjoint] couplings
CC  The number of permutations of N_c limits the code
CC  to N_c = 10.
CC
CC  The coordinate superindex is
CC  i = ix + (iy-1)*nx + (iz-1)*nx*ny + (it-1)*nx*ny*nt
CC  written by Massimo D'Elia
CC  version 1.0 - 2007-2008
CC
CC===================================================== 

      implicit none
      include "parameters.f"

c     common block definitions
      integer ncfact,perm,sign
      common/perm/ncfact,perm(nperm,ncol),sign(nperm)

      integer nran,iad_flag
      real beta_f,beta_a,rho,pi
      common/param/nran,beta_f,beta_a,rho,pi,iad_flag

      real mass,mass2,residue
      common/param2/mass,mass2,residue

      real dt_md
      integer nstep_md
      complex ieps,iepsq,ieps2q,ieps3q
      common/param3/dt_md,nstep_md,ieps,iepsq,ieps2q,ieps3q

      real remu
      real emu,ememu
      common/param4/remu,emu,ememu

      integer forweo2,backeo2
      common/geo2/forweo2(nvol,4),backeo2(nvol,4)

      integer sind,coor,forw,back
      common/ge/sind(nx,ny,nz,nt),coor(nvol,4),forw(nvol,4),back(nvol,4)

      integer sindeo,sindeoh,parity,cooreo,forweo,backeo
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),parity(nx,ny
     $     ,nz,nt),cooreo(nvol,4),forweo(nvol,4),backeo(nvol,4)

      integer eotolex,lextoeo
      common/gegeo/eotolex(nvol),lextoeo(nvol)

      complex u
      common/field/u(ncol,ncol,4,nvol)

      complex usave
      common/field1/usave(ncol,ncol,4,nvol)

      integer eta
      common/stagphase/eta(4,nvol)

      complex staple
      common/staple/staple(ncol,ncol,4,nvol)

      complex p_w,ipdot
      common/momenta/p_w(ncol,ncol,4,nvol),ipdot(ncol,ncol,4,nvol)

      complex phi_e,phi_o,chi_e,chi_o
      common/pseudof/phi_e(ncol,nvolh),phi_o(ncol,nvolh),chi_e(ncol
     $     ,nvolh),chi_o(ncol,nvolh)

      complex ue
      common/fielde/ue(4,nvol)

      real az_old_u,az_old_h,az_old_q1,az_old_q2
      common/azione_old/az_old_u(nvol),az_old_h(nvol),
     $     az_old_q1(nvol),az_old_q2(nvol)

      real az_new_u,az_new_h,az_new_q1,az_new_q2
      common/azione_new/az_new_u(nvol),az_new_h(nvol),
     $     az_new_q1(nvol),az_new_q2(nvol)

c     ACTION
c     approximation of (x)^(-1/2)
      integer nterm_a
      real cost_a,pole_a,coef_a
      common/rhmc_a/nterm_a,cost_a,pole_a(nmr),coef_a(nmr)

c     HEAT-BATH
c     approximation of (x)^(1/4)
      integer nterm_h
      real cost_h,pole_h,coef_h
      common/rhmc_h/nterm_h,cost_h,pole_h(nmr),coef_h(nmr)

c     ORGINAL EXPANSION
c     these are original expansion, dinamically scaled in the program
      real scala,min_rhmc,max_rhmc
      real ocost_a,opole_a,ocoef_a,z_a
      real ocost_h,opole_h,ocoef_h,z_h
      common/rhmc_o/scala,min_rhmc,max_rhmc,
     $     ocost_a,opole_a(nmr),ocoef_a(nmr),z_a,
     $     ocost_h,opole_h(nmr),ocoef_h(nmr),z_h


c     internal variables
      character(LEN=30) file_latt

      integer time_to_run
      real tiniz,tinistep,tfinal
      integer imposed_close

      integer termalizza,scala_ogni_step
      real dt_md_save

      integer n_rand,n_rand_save ! number of random vector to use

      integer n_traj
      integer acc

      integer metro_test
      common/com_metro/metro_test

c-------------------------------------------
c initialization operations
c-------------------------------------------     

      call cpu_time(tiniz)
      tinistep=0

      call init(time_to_run,n_rand,n_traj,termalizza)

      n_rand_save=n_rand
      dt_md_save=dt_md

c-------------------------------------------------------
c START MONTE CARLO and MEASUREMENTS
c-------------------------------------------------------
      
 300  n_traj=n_traj+1
      
c     decides if to termalize or not                                                                                                                        
      if(n_traj.le.termalizza) then
         n_rand=1
         scala_ogni_step=1
         dt_md=dt_md_save/2
         metro_test=0
      else
         n_rand=n_rand_save
         scala_ogni_step=0
         dt_md=dt_md_save
         metro_test=1
      endif

      write(*,*) "trajec:",1
      write(*,*) "n_step:",nstep_md
      write(*,*) "n_vect:",n_rand
      write(*,*) "metropolis:",metro_test

c     generate a new coniguration
      call rhmc_step(acc,n_traj,scala_ogni_step)

c     measure 
      call measure(acc,n_traj,n_rand)

c     save the generated field
c      call write_lattice2(n_traj,lattice2)

c     save one separate configuration each: 'save_each'
      if(int(n_traj/save_each)*save_each.eq.n_traj) then
         write(unit=file_latt,fmt=*) n_traj
         file_latt='lattice_save/'//file_latt
         call write_lattice2(n_traj,file_latt)
      endif

c     decides if to close or to continue
      call cpu_time(tfinal)
      if(debug.ge.1) then
         write(*,*) "Total lasted time for this config: ",int(tfinal
     $        -tinistep)
         write(*,*)
         write(*,*)
         write(*,*)
      endif
      tinistep=tfinal
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

      include "init.f"
      include "generic_subroutines.f"
      include "sun_subroutines.f"
      include "sun_ext_subroutines.f"
      include "action.f"
      include "dirac_matrix.f"
      include "force.f"
      include "updating_subroutines.f"
      include "measure_subroutines.f"
      include "vecn_subroutines.f"
      include "inverter.f"
      include "gauss_su3.f"
      include "eigenvalue.f"
c      include "test.f"
      
