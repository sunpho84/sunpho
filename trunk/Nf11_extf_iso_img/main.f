CC=====================================================
      program su_n_1_1_extf
CC
CC  Program for simulating SU(3) gauge theories
CC  with 1+1 staggered flavors in presence of imaginary 
CC  chemical potential coupled to each quark, a 
CC  static electromagnetic background field, and a real 
CC  chemical potential coupled to isospin
CC
CC  The coordinate superindex is
CC  i = ix + (iy-1)*nx + (iz-1)*nx*ny + (it-1)*nx*ny*nt
CC  written by Massimo D'Elia
CC  version 1.0 - 2007-2008
CC
CC===================================================== 
      implicit none
#include "parameters.f"

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

 300  n_traj=n_traj+1
      
c     decides if to termalize or not                                                                                                                        
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

c     generate a new coniguration
      call rhmc_step(acc,n_traj,scala_ogni_step)

c     save one configuration each: 'save_each'
c      if(int(n_traj/save_each)*save_each.eq.n_traj) then
c     FINIRE
c         call write_lattice_save(n_traj)
c      endif

c     measure 
      call measure(acc,n_traj,n_rand)

c     decides if to close or to continue
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

#include "init.f"
#include "generic_subroutines.f"
#include "sun_subroutines.f"
#include "sun_ext_subroutines.f"
#include "action.f"
#include "dirac_matrix.f"
#include "force.f"
#include "updating_subroutines.f"
#include "measure_subroutines.f"
#include "vecn_subroutines.f"
#include "inverter.f"
#include "gauss_su3.f"
#include "eigenvalue.f"
#include "test.f"
