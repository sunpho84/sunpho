c=========================================================
      subroutine m2d(y_e,y_o,x_e,x_o,emu,ememu,mass2)
c  written by Sanfo
c  calculate y=(m^2 - DEO*DOE)*x
c  that is the even\even part of M^+M
c=========================================================
      implicit none
      include "parameters.f"

c     arguments
      complex y_e(ncol,nvolh),x_e(ncol,nvolh)
      complex y_o(ncol,nvolh),x_o(ncol,nvolh)
      real emu,ememu
      real mass2

c     internal variables
      integer ivol,icol
      complex h_e(ncol,nvolh),h_o(ncol,nvolh)
      real massh
      
      massh=0.5*mass2**0.5
      
      call D(OE,0,h_o(1,1),x_e(1,1),emu,ememu)
      call D(OE,1,y_e(1,1),h_o(1,1),emu,ememu)

      call D(EO,0,h_e(1,1),x_o(1,1),emu,ememu)
      call D(EO,1,y_o(1,1),h_e(1,1),emu,ememu)
      
      call DpDdag(OE,h_o(1,1),x_e(1,1),emu,ememu)
      call DpDdag(EO,h_e(1,1),x_o(1,1),emu,ememu)

      do ivol=1,nvolh
         do icol=1,ncol
            y_e(icol,ivol)=mass2*x_e(icol,ivol)-0.25*y_e(icol,ivol)
     $           +massh*h_e(icol,ivol)
            y_o(icol,ivol)=mass2*x_o(icol,ivol)-0.25*y_o(icol,ivol)
     $           +massh*h_o(icol,ivol)
         enddo
      enddo   

      return
      end




c=========================================================
      subroutine D(eooe,dag,w,v,emu,ememu)
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
      real emu,ememu
      integer dag

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
      complex w1(ncol),w2(ncol),w3(ncol)
      
      do base_ivol=1,nvolh

c     select which part of D to apply
         if((dag.eq.0.and.eooe.eq.EO).or.(dag.eq.1.and.eooe.eq.OE)) then
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
            if(dag.eq.0) then
               w(icol,base_ivol)=w1(icol)+emu*w2(icol)-ememu*w3(icol) 
            else
               w(icol,base_ivol)=w1(icol)+ememu*w2(icol)-emu*w3(icol) 
            endif
         enddo
         
      enddo
      
      return
      end

c=========================================================
      subroutine DpDdag(eooe,w,v,emu,ememu)
c     DOEDAG plus DEO: appears in Mdag M 
c=========================================================
      implicit none
      include "parameters.f"

c     arguments
      integer eooe
      complex w(ncol,nvolh),v(ncol,nvolh)
      real emu,ememu
      
c     common blocks
      integer sindeo,sindeoh,parity,cooreo,forweo,backeo
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),parity(nx,ny
     $     ,nz,nt),cooreo(nvol,4), forweo(nvol,4),backeo(nvol,4)
      integer forweo2,backeo2
      common/geo2/forweo2(nvol,4),backeo2(nvol,4)
      complex u
      common/field/u(ncol,ncol,4,nvol)

c     internal variables
      integer base_ivol,ivol,icol
      complex w2(ncol),w3(ncol)
      real twosinhmu
      
      twosinhmu=emu-ememu

      do base_ivol=1,nvolh

c     select which part of D to apply
         if(eooe.eq.EO) then
            ivol=base_ivol
         else
            ivol=base_ivol+nvolh
         endif
         
         call vmult(w2,u(1,1,4,ivol),v(1,forweo(ivol,4)))
         call vhmult(w3,u(1,1,4,backeo2(ivol,4)),v(1,backeo(ivol,4)))

        do icol=1,ncol
           w(icol,base_ivol)=twosinhmu*(w2(icol)+w3(icol))
        enddo
        
      enddo

      return
      end
