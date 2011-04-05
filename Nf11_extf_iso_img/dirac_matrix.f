c===================================================
      subroutine m2d(y,x,potc,mass)
c  written by Sanfo
c  calculate y=(m^dag m)*x
c===================================================
      implicit none
#include "parameters.f"

c     arguments
      complex y(ncol,nvol_eo),x(ncol,nvol_eo)
      complex potc
      real mass

c     internal variables
      integer ivol,icol
      complex h(ncol,nvolh)
      real mass2

      mass2=mass**2

c     even part
      call D(OE,    h,x,potc)
      call D(OE_DAG,y,h,potc)
#ifdef ficp
      call D(EO_P_OE_DAG,h,x(1,nvolh+1),potc)
#endif

      do ivol=1,nvolh
         do icol=1,ncol
            y(icol,ivol)=mass2*x(icol,ivol)+y(icol,ivol)
#ifdef ficp
     $           +mass*h(icol,ivol)
#endif
         enddo
      enddo


#ifdef ficp
c     odd part, to be calculated only if EO improvement off
      call D(EO,    h,x(1,nvolh+1),potc)
      call D(EO_DAG,y(1,nvolh+1),h,potc)
      call D(OE_P_EO_DAG,h,x,potc)

      do ivol=1,nvolh
         do icol=1,ncol
            y(icol,ivol+nvolh)=mass2*x(icol,ivol+nvolh)+y(icol,ivol
     $           +nvolh)+mass*h(icol,ivol)
         enddo
      enddo
#endif

      return
      end

c=========================================================
      subroutine D(eooe,w,v,potc)
c  written by Massimo D'Elia, unified by Zumbo
c  apply the even/odd part or the odd/even part of the Dirac
c  Matrix, depending on the value of the passed variable eooe
c  OE and EO variable defined in "parameters.f"
c  version 1.0 - 
c=========================================================
      implicit none
#include "parameters.f"

c     arguments
      integer eooe
      complex w(ncol,nvolh),v(ncol,nvolh)
      complex potc

c     common blocks
      integer sindeo,sindeoh,parity,cooreo,forweo,backeo
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),parity(nx,ny
     $     ,nz,nt),cooreo(nvol,4), forweo(nvol,4),backeo(nvol,4)
      integer forweo2,backeo2
      common/geo2/forweo2(nvol,4),backeo2(nvol,4)
      complex u
      common/field/u(ncol,ncol,4,nvol)

c     internal variables
      integer base_ivol,ivol,icol,idir
      complex w1(ncol)
#ifndef isotropic
      complex w2(ncol),w3(ncol)
#endif
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

c     select which part of D to apply
         if(eooe.eq.EO.or.eooe.eq.OE_DAG.or.eooe.eq.EO_P_OE_DAG) then
            ivol=base_ivol
         else
            ivol=base_ivol+nvolh
         endif

         if(eooe.ne.EO_P_OE_DAG.and.eooe.ne.OE_P_EO_DAG) then

c     reset temporary variable
            do icol=1,ncol
               w1(icol)=cmplx(0,0)
            enddo

#ifdef isotropic
            do idir=1,4
#else
            do idir=1,3
#endif
               call vmult_add(w1,u(1,1,idir,ivol),v(1,forweo(ivol
     $              ,idir)))
               call vhmult_sub(w1,u(1,1,idir,backeo2(ivol,idir)), v(1
     $              ,backeo(ivol,idir)))
            enddo

         endif
#ifndef isotropic   
         call vmult(w2,u(1,1,4,ivol),v(1,forweo(ivol,4)))
         call vhmult(w3,u(1,1,4,backeo2(ivol,4)),v(1,backeo(ivol,4)))
#endif

c     write final value
         do icol=1,ncol
            w(icol,base_ivol)=A*w1(icol)
#ifndef isotropic
     $           +B*w2(icol)+C*w3(icol)
#endif
         enddo

      enddo

      return
      end
