c=========================================================
      subroutine singol_inverter(chi_e,chi_o,phi_e,phi_o,emu,ememu,pole)
c  written by Massimo D'Elia
c  version 1.0 - 2007/2008
c  chi = (MM~)^-1 phi
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"

c     argomenti
      complex chi_e(ncol,nvolh),phi_e(ncol,nvolh)
      complex chi_o(ncol,nvolh),phi_o(ncol,nvolh)
      real emu,ememu
      real pole
      
c     variabili passate da common block
      real mass,mass2,residue
      common/param2/mass,mass2,residue

c     variabili interne
      complex s_e(ncol,nvolh),s_o(ncol,nvolh)
      complex p_e(ncol,nvolh),p_o(ncol,nvolh)
      complex r_e(ncol,nvolh),r_o(ncol,nvolh)
      complex c1
      real delta,omega,lambda,gammag
      integer iter,riter

c     parametri
      parameter(sigma=1) 
      parameter(niter=2000,rniter=5)
      parameter(one_q=0.25)

      do ivol=1,nvolh
         do icol=1,ncol
            chi_e(icol,ivol)=cmplx(0,0)
            chi_o(icol,ivol)=cmplx(0,0)
         enddo
      enddo   

      riter=0

 1200 continue    !! loop over congrad to verify truerest
      
      call m2d(s_e,s_o,chi_e,chi_o,emu,ememu,pole**0.5)
      
      delta=0
      do ivol=1,nvolh
         do icol=1,ncol
            c1=phi_e(icol,ivol)-s_e(icol,ivol)
            p_e(icol,ivol)=c1
            r_e(icol,ivol)=c1
            delta=delta+real(c1)**2+aimag(c1)**2

            c1=phi_o(icol,ivol)-s_o(icol,ivol)
            p_o(icol,ivol)=c1
            r_o(icol,ivol)=c1
            delta=delta+real(c1)**2+aimag(c1)**2
         enddo
      enddo   

!!-------------------------------------------------------------------------

      iter=0

 1201 continue

      iter=iter+1
      
      call m2d(s_e,s_o,p_e,p_o,emu,ememu,pole)
      
      alpha=0
      do ivol=1,nvolh
         do icol=1,ncol
            alpha=alpha+real(s_e(icol,ivol)*conjg(p_e(icol,ivol)))
            alpha=alpha+real(s_o(icol,ivol)*conjg(p_o(icol,ivol)))
         enddo
      enddo   
      
      omega=delta/alpha
      lambda=0

      do ivol=1,nvolh
         do icol=1,ncol
            chi_e(icol,ivol)=chi_e(icol,ivol)+omega*p_e(icol,ivol)
            c1=r_e(icol,ivol)-omega*s_e(icol,ivol)
            r_e(icol,ivol)=c1
            lambda=lambda+real(c1)**2+aimag(c1)**2

            chi_o(icol,ivol)=chi_o(icol,ivol)+omega*p_o(icol,ivol)
            c1=r_o(icol,ivol)-omega*s_o(icol,ivol)
            r_o(icol,ivol)=c1
            lambda=lambda+real(c1)**2+aimag(c1)**2

         enddo
      enddo 
      
      gammag=lambda/delta
      delta=lambda
      do ivol=1,nvolh
         do icol=1,ncol
            p_e(icol,ivol)=r_e(icol,ivol)+gammag*p_e(icol,ivol)
            p_o(icol,ivol)=r_o(icol,ivol)+gammag*p_o(icol,ivol)
         enddo
      enddo 
      
c      write(*,*) iter,lambda
      if(lambda.gt.residue.and.iter.lt.niter) goto 1201
      
      call m2d(s_e,s_o,chi_e,chi_o,emu,ememu,pole)
      
      lambda=0
      do ivol=1,nvolh
         do icol=1,ncol
            c1=phi_e(icol,ivol)-s_e(icol,ivol)
            lambda=lambda+real(c1)**2+aimag(c1)**2

            c1=phi_o(icol,ivol)-s_o(icol,ivol)
            lambda=lambda+real(c1)**2+aimag(c1)**2
         enddo
      enddo   
  
      riter=riter+1    
      
      if(lambda.gt.residue.and.riter.lt.rniter) goto 1200

      return
      end


c=========================================================
      subroutine multi_shift_inverter(x_e,x_o,b_e,b_o,emu,ememu,pole
     $     ,nterm)
c  written by Sanfo
c  chiama il corretto inverter
c=========================================================
      implicit none
      include "parameters.f"

c     argomenti
      complex x_e(ncol,nvolh,nmr),b_e(ncol,nvolh)
      complex x_o(ncol,nvolh,nmr),b_o(ncol,nvolh)
      real emu,ememu
      real pole(nmr)
      integer nterm
      
      if(alg_inv.eq.i_singol) then
         call multi_shift_fuffa_inverter(x_e,x_o,b_e,b_o,emu,ememu,pole
     $        ,nterm)
      endif

      if(alg_inv.eq.i_multi) then
         call multi_shift_true_inverter(x_e,x_o,b_e,b_o,emu,ememu,pole
     $        ,nterm)
      endif
      
      return
      end


c=========================================================
      subroutine multi_shift_true_inverter(x_e,x_o,b_e,b_o,emu,ememu
     $     ,pole,nterm)
c  written by Sanfo
c  x_s = (m+X)^-1 b
c  lo calcola usando il multi-shift conjugate gradient
c  in cui X=-0.25*DOE*DEO
c  la variabile cost contiene la costante additiva
c  il vettore pole contiene i poli
c  il vettore coef contiene i coefficienti
c  nterm sono i termini effettivi (max passato dal file parameters.f)
c  NB: la massa fisica è già inclusa nei poli
c=========================================================
      implicit none
      include "parameters.f"

c     argomenti
      complex x_e(ncol,nvolh,nmr),b_e(ncol,nvolh)
      complex x_o(ncol,nvolh,nmr),b_o(ncol,nvolh)
      real emu,ememu
      real pole(nmr)
      integer nterm

c     variabili passate da common block
      real mass,mass2,residue
      common/param2/mass,mass2,residue

c     variabili interne
      complex s_e(ncol,nvolh),s_o(ncol,nvolh)
      complex r_e(ncol,nvolh),r_o(ncol,nvolh)
      complex p_e(ncol,nvolh),p_o(ncol,nvolh)
      complex ps_e(ncol,nvolh,nmr),ps_o(ncol,nvolh,nmr)
      real zps(nmr),zas(nmr),zfs(nmr),betas(nmr),alphas(nmr)
      complex c1
      real rr,rfrf,pap,alpha
      real betap,betaa
      integer iterm
      integer icol,ivol
      integer iter
      integer flag(nmr)

c     parametri
      real one_q
      integer niter
      parameter(one_q=0.25)
      parameter(niter=2000)

      if(debug.ge.1) then 
         write(*,*) " Use the MULTI-SHIFT inverter"
      endif

!     x=0
      do iterm=1,nterm
         do ivol=1,nvolh
            do icol=1,ncol
               x_e(icol,ivol,iterm)=0
               x_o(icol,ivol,iterm)=0
            enddo
          enddo
          flag(iterm)=1
      enddo

!     -p=b
!     -r=b
!     -calcolo Rr=(r,r)
      rr=0
      do ivol=1,nvolh
         do icol=1,ncol
            c1=b_e(icol,ivol)
            p_e(icol,ivol)=c1
            r_e(icol,ivol)=c1
            rr=rr+real(c1)**2+aimag(c1)**2

            c1=b_o(icol,ivol)
            p_o(icol,ivol)=c1
            r_o(icol,ivol)=c1
            rr=rr+real(c1)**2+aimag(c1)**2
         enddo
      enddo

c      if(debug.ge.2) then 
c         write(*,*) "RR:",rr
c      endif

!     -betaa=1
      betaa=1

!     -ps=b
!     -zps=zas=1
!     -alphas=0
      do iterm=1,nterm
         do ivol=1,nvolh
            do icol=1,ncol
               ps_e(icol,ivol,iterm)=b_e(icol,ivol)
               ps_o(icol,ivol,iterm)=b_o(icol,ivol)
            enddo
         enddo
         zps(iterm)=1
         zas(iterm)=1
         alphas(iterm)=0
      enddo

!     -alpha=0
      alpha=0
      
!----------------------------------------------------------------------
      
      do iter=1,niter

!     -s=Ap
         call m2d(s_e,s_o,p_e,p_o,emu,ememu,mass)
      
!     -pap=(p,s)=(p,Ap)
         pap=0
         do ivol=1,nvolh
            do icol=1,ncol
               pap=pap+real(s_e(icol,ivol)*conjg(p_e(icol,ivol)))
               pap=pap+real(s_o(icol,ivol)*conjg(p_o(icol,ivol)))
            enddo
         enddo   
         
!     calcola betaa=rr/pap=(r,r)/(p,Ap)
         betap=betaa
         betaa=-rr/pap

!     calcola 
!     -zfs
!     -betas
!     -x
         do iterm=1,nterm
         if(flag(iterm).eq.1) then
            zfs(iterm)=
     $           zas(iterm)*zps(iterm)*betap/
     $           (betaa*alpha*(zps(iterm)-zas(iterm))+
     $           zps(iterm)*betap*(1-pole(iterm)*betaa))
            betas(iterm)=betaa*zfs(iterm)/zas(iterm)
c            if(debug.ge.2) then
c               write(*,*) "Z",iterm,zfs(iterm),zas(iterm)
c               write(*,*) "BETAS",iterm,betas(iterm)
c            endif
            do ivol=1,nvolh
               do icol=1,ncol
                  x_e(icol,ivol,iterm)=x_e(icol,ivol,iterm)-
     $                 betas(iterm)*ps_e(icol,ivol,iterm)
                  x_o(icol,ivol,iterm)=x_o(icol,ivol,iterm)-
     $                 betas(iterm)*ps_o(icol,ivol,iterm)
               enddo
            enddo
         endif
         enddo
         
!     calcolo
!     -r'=r+betaa*s=r+beta*Ap
!     -rfrf=(r',r')
         rfrf=0
         do ivol=1,nvolh
            do icol=1,ncol
               c1=r_e(icol,ivol)+betaa*s_e(icol,ivol)
               r_e(icol,ivol)=c1
               rfrf=rfrf+real(c1)**2+aimag(c1)**2

               c1=r_o(icol,ivol)+betaa*s_o(icol,ivol)
               r_o(icol,ivol)=c1
               rfrf=rfrf+real(c1)**2+aimag(c1)**2
            enddo
         enddo

!     calcola alpha=rfrf/rr=(r',r')/(r,r)
         alpha=rfrf/rr
c         if(debug.ge.2) then 
c            write(*,*) "ALPHA",alpha
c         endif
         
!     calcola p'=r'+alpha*p
         do ivol=1,nvolh
            do icol=1,ncol
               p_e(icol,ivol)=r_e(icol,ivol)+alpha*p_e(icol,ivol)
               p_o(icol,ivol)=r_o(icol,ivol)+alpha*p_o(icol,ivol)
            enddo 
         enddo

!     calcola alphas=alpha*zfs*betas/zas*beta
         do iterm=1,nterm
         if(flag(iterm).eq.1) then
            alphas(iterm)=alpha*zfs(iterm)*betas(iterm)/
     $           (zas(iterm)*betaa)
c            if(debug.ge.2) then 
c               write(*,*) "ALPHAS",iterm,alphas(iterm)
c            endif
!     calcola ps'=r'+alpha*p
            do icol=1,ncol
               do ivol=1,nvolh
                  ps_e(icol,ivol,iterm)=zfs(iterm)*r_e(icol,ivol)+
     $                 alphas(iterm)*ps_e(icol,ivol,iterm)
                  ps_o(icol,ivol,iterm)=zfs(iterm)*r_o(icol,ivol)+
     $                 alphas(iterm)*ps_o(icol,ivol,iterm)
               enddo
            enddo

!     check sui residui
            if(rr*zfs(iterm).lt.residue) then
               flag(iterm)=0
            endif

!     passa tutti gli f in quelli attuali
            zps(iterm)=zas(iterm)
            zas(iterm)=zfs(iterm)
         endif
         enddo         
         rr=rfrf
c         write(*,*) rfrf
!     check sul residuo
         if(rfrf<residue) exit
      enddo

      if(debug.ge.1) then 
         write(*,*) "Iteration to invert:",iter
      endif

      return
      end


c=========================================================
      subroutine multi_shift_summed_inverter
     $     (x_e,x_o,b_e,b_o,emu,ememu,cost,pole,coef,nterm)
c  written by Sanfo
c  version 1.0 - 2008
c  somma i vari termini calcolati sopra con i loro coefficienti
c=========================================================
      implicit none
      include "parameters.f"

c     argomenti
      complex x_e(ncol,nvolh),x_o(ncol,nvolh)
      complex b_e(ncol,nvolh),b_o(ncol,nvolh)
      real emu,ememu
      real cost
      real pole(nmr),coef(nmr)
      integer nterm

c     variabili interne
      complex s_e(ncol,nvolh,nmr),s_o(ncol,nvolh,nmr)
      integer ivol,icol,iterm

      if(debug.ge.1) then 
         write(*,*) " Use the MULTI_SHIFT_SUMMED inverter"
      endif
      
      call multi_shift_inverter(s_e,s_o,b_e,b_o,emu,ememu,pole,nterm)

      do ivol=1,nvolh
         do icol=1,ncol
            x_e(icol,ivol)=cost*b_e(icol,ivol)
            x_o(icol,ivol)=cost*b_o(icol,ivol)
            do iterm=1,nterm
               x_e(icol,ivol)=x_e(icol,ivol)+coef(iterm)*s_e(icol,ivol
     $              ,iterm)
               x_o(icol,ivol)=x_o(icol,ivol)+coef(iterm)*s_o(icol,ivol
     $              ,iterm)
            enddo
         enddo
      enddo            

      return
      end


c=========================================================
      subroutine multi_shift_fuffa_inverter(x_e,x_o,b_e,b_o,emu,ememu
     $     ,pole,nterm)
c  written by Sanfo
c     version 1.0 - 2007/2008
c  esegue lo shift termine a termine e lo schiaffain x:
c  x_i = (pole_i+MM~)^-1 b
c=========================================================
      implicit none
      include "parameters.f"

c     argomenti
      complex x_e(ncol,nvolh,nmr),x_o(ncol,nvolh,nmr)
      complex b_e(ncol,nvolh,nmr),b_o(ncol,nvolh,nmr)
      real emu,ememu
      real cost
      real pole(nmr),coef(nmr)
      integer nterm

c     variabili interne
      integer iterm
      
      do iterm=1,nterm
         call singol_inverter(x_e(1,1,iterm),x_o(1,1,iterm),b_e,b_o,emu
     $        ,ememu,pole(iterm))
      enddo
      
      return
      end

