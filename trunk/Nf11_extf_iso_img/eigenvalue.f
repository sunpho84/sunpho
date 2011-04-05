c=========================================================
      subroutine eigenval(autmin,autmax,potc,ud)
c  written by Sanfo
c  find minimum and maximum eigenvalue of even\even part of M^+M
c  by using an approximated (and probably not stable) algorithm
c  This is used to find the appropriate interval where to scale the
c  rational expansion
c=========================================================
      implicit none
#include "parameters.f"

c     arguments
      real autmin,autmax
      complex potc
      integer ud

c     common blocks
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue

c     internal variables
      integer aut
      complex y(ncol,nvol_eo),x(ncol,nvol_eo)
      real norm
      integer icol,ivol,iter
      integer niter

c     parameters
      integer met,sane,marc
      real sigma
      parameter(sigma=1,sane=1,marc=2,met=marc)


#ifdef mf
      call add_extf(ud)
#else
      ud=ud !to avoid warning
#endif
      do aut=1,2
         do ivol=1,nvol_eo
            call gauss_vector(x(1,ivol),sigma)
         enddo

         if(met.eq.marc.and.aut.eq.2) then
            niter=5
         else
            if(aut.eq.1) then
               niter=50
            else
               niter=50
            endif
         endif
         
         do iter=1,niter
            
            if(aut.eq.1) then
               call m2d(y,x,potc,mass)
            else
               if(met.eq.sane) then
                  call singol_inverter(y,x,potc,mass)
               else
                  call m2d(y,x,potc,mass)
                  do ivol=1,nvol_eo
                     do icol=1,ncol
                        y(icol,ivol)=autmax*x(icol,ivol)-y(icol,ivol)
                     enddo
                  enddo
               endif
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
            
#ifdef debug2
            if(aut.eq.1) then
               write(*,*) "EIGEN_STEP:",iter,"NORM:",norm
            else
               write(*,*) "EIGEN_STEP:",iter,"NORM:",1/norm
            endif

            call flush(6)
            call flush(7)
#endif
            
         enddo
         
         if(aut.eq.1) then 
            autmax=norm
         else
            if(met.eq.sane) then
               autmin=1/norm
            else
               autmin=autmax-norm
            endif
         endif

      enddo
      

#ifdef mf
      call rem_extf(ud)
#endif
      
      return
      end


c=========================================================
      subroutine scale_rhmc
c  written by Sanfo
c  version 1.0 - 05/08
c  riscala i parametri dell'espansione adattandoli
c  all'intervallo opportuno.
c=========================================================
      implicit none
#include "parameters.f"

c     common blocks
      complex potc1,potc2
      common/param3/potc1,potc2

c     AZIONE
c     approssimazione di (m+x)^(-nf/4)
c     deve essere usato per il calcolo dell'azione, che avviene qui,
c     per il calcolo di Polyakov e quant'altro, che avviene in CHIRAL
c     e per le equazioni della dinamica, che avvengono in compute_ipdot
      integer nterm_a
      real cost_a,pole_a,coef_a
      common/rhmc_a/nterm_a,cost_a,pole_a(nmr),coef_a(nmr)

c     HEAT-BATH
c     approssimazione di (m+x)^(nf/8)
c     viene usato per la generazione dei campi phi tramite heat-bath
c     nella funzione: create_phi
      integer nterm_h
      real cost_h,pole_h,coef_h
      common/rhmc_h/nterm_h,cost_h,pole_h(nmr),coef_h(nmr)

c     ESPANSIONI ORIGINALI DELL'RHMC
c     sono quelle che poi vengono scalate
      real scala,min_rhmc,max_rhmc
      real ocost_a,opole_a,ocoef_a,z_a
      real ocost_h,opole_h,ocoef_h,z_h
      common/rhmc_o/scala,min_rhmc,max_rhmc,
     $     ocost_a,opole_a(nmr),ocoef_a(nmr),z_a,
     $     ocost_h,opole_h(nmr),ocoef_h(nmr),z_h

c     variabili interne
      real autmin1,autmax1
      real autmin2,autmax2
      real autmin,autmax
      real scala_az,scala_hz
      integer iterm

c     calcola gli autovalori massimo e minimo per le due matrici separatamente
      call eigenval(autmin1,autmax1,potc1,1)
      call eigenval(autmin2,autmax2,potc2,2)

c     trova il minimo ed il massimo "assoluti"
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

#ifdef debug1
         write(*,*) "Eigenvalue:",autmin,autmax,autmax/autmin,z_a,z_h
#endif
c     setta i fattori di scala
      scala=autmax*1.1
      scala_hz=scala**z_h
      scala_az=scala**z_a

c     check: stan nell'intervallo?
      if(((autmin/scala).lt.min_rhmc).or.((autmax/scala).gt.max_rhmc)) 
     $     then
         write(*,*) "ATTENTION! Approximation not valid!!!"
         write(*,*) "Interval of validity of approximation:",
     $        min_rhmc,max_rhmc
         write(*,*) "Interval of eigenvalue:",
     $        autmin,autmax
      endif

c     scala il tutto
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
