c=========================================================
      subroutine rhmc_step(acc,n_traj,scale_each_step)
c  written by Sanfo
c  version 1.0 - 04/08
c  exec singolo rhmc step
c=========================================================
      implicit none
#include "parameters.f"

c     arguments
      integer acc
      integer n_traj
      integer scale_each_step

c     common-block passed variables
      complex u,usave
      common/field/u(ncol,ncol,4,nvol)
      common/field1/usave(ncol,ncol,4,nvol)
      complex ieps,iepsq,ieps2q,ieps3q
      real dt_md
      integer nstep_md
      common/param2/dt_md,nstep_md,ieps,iepsq,ieps2q,ieps3q
c     RHMC ACTION PART
      real az_old,az_new
      common/azione_loc/az_old(nvol),az_new(nvol)

c     internal variables
      real difact,edifact
      real ran2,xrand
      integer temp

#ifdef debug1
      write(*,*)
      write(*,*)
      write(*,*) "Init trajectory ",n_traj
#endif

      open(u_random, file='n_steps_tune', status='unknown')
      read(u_random,*,end=111,err=111) temp
      nstep_md=temp
      go to 112
 111  write(u_random,*) nstep_md
      go to 113

 112  xrand=0
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

 113  call create_momenta
      call create_phi
      
      call action(az_old)
      call dinamica(scale_each_step)
      call action(az_new)

      call dif_action(difact,az_new,az_old)
      

!     metropolis test
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

c      call write_lattice2(n_traj)

      return
      end


c=========================================================
      subroutine create_momenta()
c  written by Massimo D'Elia
c  version 1.0 - 
c=========================================================
      implicit none
#include "parameters.f"

c     common-block passed variables
      complex p_w,ipdot
      common/momenta/p_w(ncol,ncol,4,nvol),ipdot(ncol,ncol,4,nvol)
      
c     internal variables
      integer ivol,idir
      
c     internal parameters
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
c  written by Massimo D'Elia
c  version 1.0 - 
c=========================================================
      implicit none
#include "parameters.f"

c     arguments
      complex phi(ncol,nvol_eo)
      complex potc
      integer ud

c     common-block passed variables
      integer nterm
      real cost,pole,coef
      common/rhmc_h/nterm,cost,pole(nmr),coef(nmr)
     
c     internal variables
      complex rnd(ncol,nvol_eo) 
      integer ivol

c     internal parameters
      real sigma
      parameter(sigma=1)

      do ivol=1,nvol_eo
         call gauss_vector(rnd(1,ivol),sigma)
      enddo

      call multi_shift_summed_inverter
     $     (phi,rnd,
     $     potc,
     $     cost,pole,coef,
     $     nterm,ud)

      return
      end

c=========================================================
      subroutine create_phi
c  written by Massimo D'Elia
c  version 1.0 - 
c=========================================================
      implicit none
#include "parameters.f"

c     common-block passed variables
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
c  written by Massimo D'Elia
c  version 1.0 - 2007/2008
c=========================================================
      implicit none
#include "parameters.f"

c     arguments
      integer scale_each_step

c     common-block passed variables
      complex ieps,iepsq,ieps2q,ieps3q
      real dt_md
      integer nstep_md
      common/param2/dt_md,nstep_md,ieps,iepsq,ieps2q,ieps3q

c     internal variables
      integer istep

#ifdef debug1
      write(*,*) "Integrator: leapfrog"
#endif
      
c----------------------------------------------------------------------------
c                  PRIMO MEZZO STEP: TROVA U(-dt/2)
c----------------------------------------------------------------------------

c     CALCOLA U(-dt/2) ovvero p(-dt/2)=p(0)-v(0)*dt/2
c     NB: v(0)=p_w è stato estratto con create_momenta
      call compute_utpdt(ieps2q)

c     normalize after first half step
      call addrem_stagphase
      call normalize1_lattice
      call addrem_stagphase

c----------------------------------------------------------------------------
c                             CICLO PRINCIPALE
c----------------------------------------------------------------------------

      do istep=1,nstep_md-1     !! starting main M.D. loop

#ifdef debug1
         write(*,*)
         write(*,*) "Inizio step",istep
#endif
         
c     CALCOLA P(t+dt) overo v(t+dt)=v(t)+a(t)*dt
         call compute_ptpdt(ieps,scale_each_step)
         
c     CALCOLA U(t+dt/2) ovvero p(t+dt/2)=p(t-dt/2)+v(t)*dt
         call compute_utpdt(ieps)
         
c     NORMALIZZA OGNI unit_each STEP
         if((istep-unit_each*(istep/unit_each)).eq.0) then
            call addrem_stagphase
            call normalize1_lattice
            call addrem_stagphase
         endif
         
      enddo                     !! enddo over M.D. steps minus one
      
c----------------------------------------------------------------------------
c     ULTIMO MEZZO STEP: RITROVA U(dt)
c----------------------------------------------------------------------------
      
c     CALCOLA P(t+dt)
      call compute_ptpdt(ieps,scale_each_step)

c     CALCOLA U(t'=t+dt)
      call compute_utpdt(ieps2q)
      
      return
      end


c=========================================================
      subroutine omelyan(scale_each_step)
c  written by Sanfo
c  version 1.0 - 2007/2008
c     questo è una modifica del leapfrog, per info vedi 
c     cond-mat/0110438v1 
c     per ogni singolo step dovrei calcolare
c
c     v1 = v(t) + a[r(t)]*lambda*dt
c     r1 = r(t) + v1*dt/2
c     v2 = v1 + a[r1]*(1 -2*lambda)*dt
c     r(t + dt) = r1 + v2*dt/2
c     v(t + h) = v2 + a[r(t + dt)]*lambda*dt
c
c     ma invece ottimizzo un poco mettendo insieme la prima e l'ultima riga
c     per fare questo devo distinguere il primo e l'ultimo step
c
c     solo 1°step: 
c     v1 = v(t) + a[r(t)]*lambda*dt
c
c     r1 = r(t) + v1*dt/2
c     v2 = v1 + a[r1]*(1 -2*lambda)*dt
c     r(t + dt) = r1 + v2*dt/2
c
c     non l'ultimo step: 
c     v1 = v2 + a[r(t + dt)]*2*lambda*dt
c     
c     ultimo step: 
c     v(t + h) = v2 + a[r(t + dt)]*lambda*dt
c=========================================================
      implicit none
#include "parameters.f"

c     arguments
      integer scale_each_step

c     common-block passed variables
      complex ieps,iepsq,ieps2q,ieps3q
      real dt_md
      integer nstep_md
      common/param2/dt_md,nstep_md,ieps,iepsq,ieps2q,ieps3q

c     internal variables
      complex ldt,dldt,udldt
      integer istep

c     parameters
      real lambda
      parameter(lambda=0.1931833)


      ldt=ieps*lambda
      dldt=2*ldt
      udldt=ieps-dldt

#ifdef debug1
      write(*,*) "Integratore: omelyan"
#endif

c----------------------------------------------------------------------------
c                  PRIMO MEZZO STEP: TROVA P(t+lambda*dt)
c----------------------------------------------------------------------------

c     CALCOLA P(t+lambda*dt) ovvero v1=v(t)+a[r(t)]*lambda*dt
      call compute_ptpdt(ldt,scale_each_step)

c     normalize after first half step
      call addrem_stagphase
      call normalize1_lattice
      call addrem_stagphase

c----------------------------------------------------------------------------
c                             CICLO PRINCIPALE
c----------------------------------------------------------------------------

      do istep=1,nstep_md   !! starting main M.D. loop

#ifdef debug1
         write(*,*) "Inizio step",istep
#endif
         
c     CALCOLA U(t+dt/2) ovvero r1=r(t)+v1*dt/2
         call compute_utpdt(ieps2q)

c     CALCOLA P(t+(1-2*lambda)*dt) overo v2=v1+a[r1]*(1-2*lambda)*dt
         call compute_ptpdt(udldt,scale_each_step)

c     CALCOLA U(t+dt/2) ovvero r(t+dt)=r1+v2*dt/2
         call compute_utpdt(ieps2q)

         if(istep.eq.nstep_md) then
c     CALCOLA P(t+dt) ovvero v(t+dt)=v2+a[r(t+dt)]*lambda*dt
            call compute_ptpdt(ldt,scale_each_step)
         else
c     CALCOLA P(t+dt) ovvero v1'=v2+a[r(t+dt)]*2*lambda*dt
            call compute_ptpdt(dldt,scale_each_step)
         endif
         
c     NORMALIZZA OGNI unit_each STEP
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
c  written by Sanfo
c  chiama l'opportuna funzione di update
c=========================================================
      implicit none

#include "parameters.f"

c     arguments
      integer scale_each_step

c     common blocks
      complex ieps,iepsq,ieps2q,ieps3q
      real dt_md
      integer nstep_md
      common/param2/dt_md,nstep_md,ieps,iepsq,ieps2q,ieps3q

      dt_md=1.0/nstep_md
      iepsq=cmplx(0.0,dt_md/4.0)
      ieps2q=2.0*iepsq
      ieps3q=3.0*iepsq
      ieps=4.0*iepsq

      write(*,*) "Checking", nstep_md,dt_md

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
c                       NORMALIZZAZIONE FINALE
c----------------------------------------------------------------------------

!     normalize before ending trajectory
      call addrem_stagphase
      call normalize_lattice
      call addrem_stagphase

      call scale_rhmc

      return
      end



c=========================================================
      subroutine copyconf(u_out,u_in)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
c=========================================================
      implicit none
#include "parameters.f"

c     arguments
      complex u_in(ncol,ncol,4,nvol),u_out(ncol,ncol,4,nvol)

c     internal variables
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
      
