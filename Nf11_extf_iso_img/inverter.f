c=========================================================
      subroutine singol_inverter(chi,phi,potc,pole)
c  written by Massimo D'Elia
c  version 1.0 - 2007/2008
c  chi = (MM~)^-1 phi
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
#include "parameters.f"

c     argomenti
      complex chi(ncol,nvol_eo),phi(ncol,nvol_eo)
      complex potc
      real pole
      
c     variabili passate da common block
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue

c     variabili interne
      complex s(ncol,nvol_eo),r(ncol,nvol_eo),p(ncol,nvol_eo)
      complex c1
      real delta,omega,lambda,gammag
      integer iter,riter

c     parametri
      parameter(sigma=1) 
      parameter(niter=2000,rniter=5)

      do ivol=1,nvol_eo
         do icol=1,ncol
            chi(icol,ivol)=cmplx(0,0)
         enddo
      enddo   

      riter=0

 1200 continue    !! loop over congrad to verify truerest

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
c  written by Sanfo
c  chiama il corretto inverter
c=========================================================
      implicit none
#include "parameters.f"

c     argomenti
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
#include "parameters.f"

c     argomenti
      complex x(ncol,nvol_eo,nmr),b(ncol,nvol_eo)
      complex potc
      real pole(nmr)
      integer nterm
      integer ud

c     variabili passate da common block
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue

c     variabili interne
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

c     parametri
      real one_q
      integer niter
      parameter(one_q=0.25)
      parameter(niter=2000)

#ifdef mf
      call add_extf(ud)
#else
      ud=ud !to avoid warning
#endif

#ifdef debug1 
      write(*,*) " Uso il MULTI-SHIFT inverter"
#endif

!     x=0
      do iterm=1,nterm
         do ivol=1,nvol_eo
            do icol=1,ncol
               x(icol,ivol,iterm)=0
            enddo
          enddo
          flag(iterm)=1
      enddo

!     -p=b
!     -r=b
!     -calcolo Rr=(r,r)
      rr=0
      do ivol=1,nvol_eo
         do icol=1,ncol
            c1=b(icol,ivol)
            p(icol,ivol)=c1
            r(icol,ivol)=c1
            rr=rr+real(c1)**2+aimag(c1)**2
         enddo
      enddo

!     -betaa=1
      betaa=1

!     -ps=b
!     -zps=zas=1
!     -alphas=0
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

!     -alpha=0
      alpha=0
      
!----------------------------------------------------------------------
      
      do iter=1,niter

!     -s=Ap
         call m2d(s,p,potc,mass)
      
!     -pap=(p,s)=(p,Ap)
         pap=0
         do ivol=1,nvol_eo
            do icol=1,ncol
               pap=pap+real(s(icol,ivol)*conjg(p(icol,ivol)))
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
            do ivol=1,nvol_eo
               do icol=1,ncol
                  x(icol,ivol,iterm)=x(icol,ivol,iterm)-
     $                 betas(iterm)*ps(icol,ivol,iterm)
               enddo
            enddo
         endif
         enddo
         
!     calcolo
!     -r'=r+betaa*s=r+beta*Ap
!     -rfrf=(r',r')
         rfrf=0
         do ivol=1,nvol_eo
            do icol=1,ncol
               c1=r(icol,ivol)+betaa*s(icol,ivol)
               r(icol,ivol)=c1
               rfrf=rfrf+real(c1)**2+aimag(c1)**2
            enddo
         enddo

!     calcola alpha=rfrf/rr=(r',r')/(r,r)
         alpha=rfrf/rr
         
!     calcola p'=r'+alpha*p
         do ivol=1,nvol_eo
            do icol=1,ncol
               p(icol,ivol)=r(icol,ivol)+alpha*p(icol,ivol)
            enddo 
         enddo

!     calcola alphas=alpha*zfs*betas/zas*beta
         do iterm=1,nterm
         if(flag(iterm).eq.1) then
            alphas(iterm)=alpha*zfs(iterm)*betas(iterm)/
     $           (zas(iterm)*betaa)
!     calcola ps'=r'+alpha*p
            do icol=1,ncol
               do ivol=1,nvol_eo
                  ps(icol,ivol,iterm)=zfs(iterm)*r(icol,ivol)+
     $                 alphas(iterm)*ps(icol,ivol,iterm)
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

!     check sul residuo
         if(rfrf<residue) exit
      enddo

#ifdef debug1 
      write(*,*) "Iterazioni impiegate:",iter
#endif

#ifdef mf
      call rem_extf(ud)
#endif

      return
      end


c=========================================================
      subroutine multi_shift_fuffa_inverter(x,b,potc,pole,nterm,ud)
c  written by Sanfo
c     version 1.0 - 2007/2008
c  esegue lo shift termine a termine e lo schiaffain x:
c  x_i = (pole_i+MM~)^-1 b
c=========================================================
      implicit none
#include "parameters.f"

c     argomenti
      complex x(ncol,nvol_eo,nmr)
      complex b(ncol,nvol_eo)
      complex potc
      real pole(nmr)
      integer nterm
      integer ud

c     variabili interne
      integer iterm

#ifdef debug1
      write(*,*) " Uso il MULTI-FUFFA inverter"
#endif

#ifdef mf
      call add_extf(ud)
#else
      ud=ud !to avoid warning
#endif
      do iterm=1,nterm
         call singol_inverter(x(1,1,iterm),b,potc,pole(iterm))
      enddo
#ifdef mf
      call rem_extf(ud)
#endif

      return
      end

c=========================================================
      subroutine multi_shift_summed_inverter
     $     (x,b,potc,cost,pole,coef,nterm,ud)
c  written by Sanfo
c  version 1.0 - 2008
c  somma i vari termini calcolati sopra con i loro coefficienti
c=========================================================
      implicit none
#include "parameters.f"

c     argomenti
      complex x(ncol,nvol_eo),b(ncol,nvol_eo)
      complex potc
      real cost
      real pole(nmr),coef(nmr)
      integer nterm
      integer ud

c     variabili interne
      complex s(ncol,nvol_eo,nmr)
      integer ivol,icol,iterm

#ifdef debug1 
      write(*,*) " Uso il MULTI_SHIFT_SUMMED inverter"
#endif
      
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


