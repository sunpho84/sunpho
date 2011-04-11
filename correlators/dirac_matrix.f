c=========================================================
      subroutine m2d(y,x,eim,emim,mass2)
c  written by Sanfo
c  calculate y=(m^2 - DEO*DOE)*x
c  that is the even\even part of M^+M
c=========================================================
      implicit none
      include "parameters.f"

c     arguments
      complex y(ncol,nvolh),x(ncol,nvolh)
      complex eim,emim
      real mass2

c     internal variables
      integer ivol,icol
      complex h(ncol,nvolh)

      call D(OE,h(1,1),x(1,1),eim,emim)
      call D(EO,y(1,1),h(1,1),eim,emim)

c      call DOE(h(1,1),x(1,1),eim,emim)
c      call DEO(y(1,1),h(1,1),eim,emim)

      do ivol=1,nvolh
         do icol=1,ncol
            y(icol,ivol)=mass2*x(icol,ivol)-0.25*y(icol,ivol)
         enddo
      enddo   

      return
      end




c=========================================================
      subroutine D(eooe,w,v,eim,emim)
c  written by Massimo D'Elia, unified by Sanfo
c  apply the even/odd part or the odd/even part of the Dirac
c  Matrix, depending on the value of the passed variable eooe
c  OE and EO variable defined in "parameters.f"
c  version 1.0 - 
c=========================================================
      implicit none
      include "parameters.f"

c     arguments
      integer eooe
      complex w(ncol,nvolh),v(ncol,nvolh)
      complex eim,emim

c     common-block passed variables
      integer sindeo,sindeoh,parity,cooreo,forweo,backeo
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),
     $     parity(nx,ny,nz,nt),cooreo(nvol,4),
     $     forweo(nvol,4),backeo(nvol,4)
      integer forweo2,backeo2
      common/geo2/forweo2(nvol,4),backeo2(nvol,4)
      complex u
      common/field/u(ncol,ncol,4,nvol)

c     internal variables
      integer base_ivol,ivol,icol,idir
      complex w1(ncol),w2(ncol),w3(ncol)


      do base_ivol=1,nvolh

c     select which part of D to apply
         if(eooe.eq.EO) then
            ivol=base_ivol
         else
            ivol=base_ivol+nvolh
         endif

c     reset temporary variable
         do icol=1,ncol
            w1(icol)=cmplx(0,0)
         enddo

         do idir=1,3
            call vmult_add(w1,u(1,1,idir,ivol),v(1,forweo(ivol,idir)))
            call vhmult_sub(w1,u(1,1,idir,backeo2(ivol,idir)),
     $           v(1,backeo(ivol,idir)))
         enddo
         
         call vmult(w2,u(1,1,4,ivol),v(1,forweo(ivol,4)))
         call vhmult(w3,u(1,1,4,backeo2(ivol,4)),v(1,backeo(ivol,4)))

c     write final value
         do icol=1,ncol
            w(icol,base_ivol)=w1(icol)+eim*w2(icol)-emim*w3(icol) 
         enddo
         
      enddo
      
      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          Follow: 2 separate function that do the same
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c=========================================================
      subroutine DEO(w,v,eim,emim)
c  written by Massimo D'Elia
c  version 1.0 - 
c=========================================================
      implicit none
      include "parameters.f"

c     arguments
      complex w(ncol,nvolh),v(ncol,nvolh)
      complex eim,emim

c     common-block passed variables
      integer sindeo,sindeoh,parity,cooreo,forweo,backeo
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),
     $     parity(nx,ny,nz,nt),cooreo(nvol,4),
     $     forweo(nvol,4),backeo(nvol,4)
      integer forweo2,backeo2
      common/geo2/forweo2(nvol,4),backeo2(nvol,4)
      complex u
      common/field/u(ncol,ncol,4,nvol)

c     internal variables
      integer ivol,icol,idir
      complex w1(ncol),w2(ncol),w3(ncol)


      do ivol=1,nvolh
         do icol=1,ncol
            w1(icol)=cmplx(0,0)
         enddo

         do idir=1,3
            call vmult_add(w1,u(1,1,idir,ivol),v(1,forweo(ivol,idir)))
            call vhmult_sub(w1,u(1,1,idir,backeo2(ivol,idir)),
     $           v(1,backeo(ivol,idir)))
         enddo
         
         call vmult(w2,u(1,1,4,ivol),v(1,forweo(ivol,4)))
         call vhmult(w3,u(1,1,4,backeo2(ivol,4)),v(1,backeo(ivol,4)))

         do icol=1,ncol
            w(icol,ivol)=w1(icol)+eim*w2(icol)-emim*w3(icol) 
         enddo
      enddo

      return
      end



c=========================================================
      subroutine DOE(w,v,eim,emim)
c  written by Massimo D'Elia
c  version 1.0 - 
c=========================================================
      implicit none
      include "parameters.f"

c     arguments
      complex w(ncol,nvolh),v(ncol,nvolh)
      complex eim,emim

c     common-block passed variables
      integer sindeo,sindeoh,parity,cooreo,forweo,backeo
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),
     $     parity(nx,ny,nz,nt),cooreo(nvol,4),
     $     forweo(nvol,4),backeo(nvol,4)
      integer forweo2,backeo2
      common/geo2/forweo2(nvol,4),backeo2(nvol,4)
      complex u
      common/field/u(ncol,ncol,4,nvol)

c     internal variables
      integer ivol,icol,idir
      complex w1(ncol),w2(ncol),w3(ncol)


      do ivol=1+nvolh,nvol
         do icol=1,ncol
            w1(icol)=cmplx(0,0)
         enddo

         do idir=1,3
            call vmult_add(w1,u(1,1,idir,ivol),v(1,forweo(ivol,idir)))
            call vhmult_sub(w1,u(1,1,idir,backeo2(ivol,idir)),
     $           v(1,backeo(ivol,idir)))
         enddo

         call vmult(w2,u(1,1,4,ivol),v(1,forweo(ivol,4)))
         call vhmult(w3,u(1,1,4,backeo2(ivol,4)),v(1,backeo(ivol,4)))

         do icol=1,ncol
            w(icol,ivol-nvolh)=w1(icol)+eim*w2(icol)-emim*w3(icol)
         enddo
      enddo

      return
      end
