c=========================================================
      subroutine eigenval(autmin,autmax,emu,ememu)
c  written by Sanfo
c  find minimum and maximum eigenvalue of even\even part of M^+M
c  by using an approximated (and probably not stable) algorithm
c  This is used to find the appropriate interval where to scale the
c  rational expansion
c=========================================================
      implicit none
      include "parameters.f"

c     arguments
      real autmin,autmax
      real emu,ememu

c     common blocks
      real mass,mass2,residue
      common/param2/mass,mass2,residue

c     internal variables
      integer aut
      complex y_e(ncol,nvolh),x_e(ncol,nvolh)
      complex y_o(ncol,nvolh),x_o(ncol,nvolh)
      real norm
      integer icol,ivol,iter
      integer niter

c     parameters
      real sigma
      parameter(sigma=1)


      do aut=1,2

         do ivol=1,nvolh
            call gauss_vector(x_e(1,ivol),sigma)
            call gauss_vector(x_o(1,ivol),sigma)
         enddo
         
         if(aut.eq.1) then
            niter=50
         else
            niter=5
         endif

         do iter=1,niter
            
            if(aut.eq.1) then
               call m2d(y_e,y_o,x_e,x_o,emu,ememu,mass2)
            else
c               call singol_inverter(y_e,y_o,x_e,x_o,emu,ememu,mass2)
               call m2d(y_e,y_o,x_e,x_o,emu,ememu,mass2)
               do ivol=1,nvolh
                  do icol=1,ncol
                     y_e(icol,ivol)=autmax*x_e(icol,ivol)-y_e(icol,ivol)
                     y_o(icol,ivol)=autmax*x_o(icol,ivol)-y_o(icol,ivol)
                  enddo
               enddo
            endif
            
            norm=0
            do ivol=1,nvolh
               do icol=1,ncol
                  norm=norm+real(y_e(icol,ivol))**2+aimag(y_e(icol
     $                 ,ivol))**2
                  norm=norm+real(y_o(icol,ivol))**2+aimag(y_o(icol
     $                 ,ivol))**2
               enddo
            enddo
            norm=(norm/nvol/ncol)**0.5
            
            do ivol=1,nvolh
               do icol=1,ncol
                  x_e(icol,ivol)=y_e(icol,ivol)/norm
                  x_o(icol,ivol)=y_o(icol,ivol)/norm
               enddo
            enddo

            if(debug.ge.2) then
               write(*,*) "EIGEN_STEP:",iter,"NORM:",norm
            endif
            
         enddo
         
         if(aut.eq.1) then 
            autmax=norm
         else
c            autmin=1/norm
            autmin=autmax-norm
         endif
         
      enddo

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
      include "parameters.f"

c     common blocks
      real remu
      real emu,ememu
      common/param4/remu,emu,ememu

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
      real autmin,autmax
      real scala_az,scala_hz
      integer iterm

c     calcola gli autovalori massimo e minimo per le due matrici separatamente
      call eigenval(autmin,autmax,emu,ememu)

      if(debug.ge.1) then
         write(*,*) "Eigenvalue:",autmin,autmax,autmax/autmin,z_a,z_h
      endif

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
