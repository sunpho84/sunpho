c=========================================================
      subroutine compute_utpdt(iepsatt)
c  written by Massimo D'Elia
c  version 1.0 - 2007/2008
c  evolve all the links for a "iepsatt" MD time
c=========================================================
      implicit none
      include "parameters.f"

c     arguments
      complex iepsatt

c     common blocks
      complex u,p_w,ipdot
      common/field/u(ncol,ncol,4,nvol)
      common/momenta/p_w(ncol,ncol,4,nvol),ipdot(ncol,ncol,4,nvol)

c     internal variables
      complex h(ncol,ncol),u1(ncol,ncol),u2(ncol,ncol)
      integer ivol,idir


      do ivol=1,nvol
         do idir=1,4
            call cmult(h,iepsatt,p_w(1,1,idir,ivol))
            call expmat(u1,h)
            call mmult(u2(1,1),u1(1,1),u(1,1,idir,ivol))
            call equal(u(1,1,idir,ivol),u2(1,1))
         enddo 
      enddo   !! end of starting step
      
      return
      end
      

c=========================================================
      subroutine compute_ptpdt(iepsatt,scale_expansion)
c  calculate p_w=p_w-iepsatt*ipdot to evolve link momenta
c  i.e calculate v(t+dt)=v(t)+a*dt
c  written by Massimo D'Elia
c  version 1.0 - 2007/2008
c=========================================================
      implicit none
      include "parameters.f"

c     arguments
      complex iepsatt
      integer scale_expansion

c     common blocks
      complex p_w,ipdot
      common/momenta/p_w(ncol,ncol,4,nvol),ipdot(ncol,ncol,4,nvol)

c     internal variables
      integer ivol,idir,icol1,icol2


      call compute_ipdot(scale_expansion)

      do ivol=1,nvol
         do idir=1,4
            do icol1=1,ncol
               do icol2=1,ncol
                  p_w(icol1,icol2,idir,ivol)=p_w(icol1,icol2,idir,ivol)-
     $                 iepsatt*ipdot(icol1,icol2,idir,ivol)
               enddo
            enddo
         enddo
      enddo

      return
      end

c=========================================================
      subroutine compute_ipdot(scale_expansion)
c  written by Massimo D'Elia
c  compute the force a(t)
c  version 1.0 - 2007/2008
c=========================================================
      implicit none
      include "parameters.f"

c     arguments
      integer scale_expansion

c     common blocks
      complex p_w,ipdot
      common/momenta/p_w(ncol,ncol,4,nvol),ipdot(ncol,ncol,4,nvol)
      complex u
      common/field/u(ncol,ncol,4,nvol)

c     internal variables
      complex ipdot_g(ncol,ncol,4,nvol),ipdot_q(ncol,ncol,4,nvol)
      integer ivol,idir,icol1,icol2
      complex u1(ncol,ncol),u2(ncol,ncol)


      if(scale_expansion.eq.1) then
         call scale_rhmc
      endif

      if(debug.ge.1) then
         write(*,*) "Starting calculation of the force"
      endif
      
      call compute_gluonic_force(ipdot_g)
      call compute_quark_force(ipdot_q)

c     sum gluonic and quark force and calculate Trace Anti-hermitian part
      do ivol=1,nvol
         do idir=1,4
            do icol1=1,ncol
               do icol2=1,ncol
                  u1(icol1,icol2)
     $                 =ipdot_g(icol1,icol2,idir,ivol)
     $                 +ipdot_q(icol1,icol2,idir,ivol)
               enddo
            enddo
            call mmult(u2(1,1),u(1,1,idir,ivol),u1(1,1))
            call TA(ipdot(1,1,idir,ivol),u2(1,1))
         enddo
      enddo

      if(debug.ge.1) then
         write(*,*) "Force calculation terminated"
      endif

      if(debug.ge.2) then
         write(*,*) "FORCE:",ipdot_g(icol1,icol2,idir,ivol),
     $        ipdot_q(icol1,icol2,idir,ivol)
      endif
      
      return
      end


c=========================================================
      subroutine compute_gluonic_force(ipdot)
c  written by Sanfo
c  calculte gluonic part of the force
c=========================================================
      implicit none
      include "parameters.f"

c     arguments
      complex ipdot(ncol,ncol,4,nvol)

c     common blocks
      integer nran,iad_flag
      real beta_f,beta_a,rho,pi
      common/param/nran,beta_f,beta_a,rho,pi,iad_flag
      complex staple
      common/staple/staple(ncol,ncol,4,nvol)

c     internal variables
      integer ivol,idir
      real r_1


      r_1=beta_f/float(ncol) 

      call compute_staples
      
      do ivol=1,nvol
         do idir=1,4
            call rmult(ipdot(1,1,idir,ivol),r_1,staple(1,1,idir,ivol))
         enddo
         
      enddo   

      return
      end

c=========================================================
      subroutine compute_quark_force(ipdot)
c  written by Sanfo
c  calculate total quark force
c=========================================================
      implicit none
      include "parameters.f"

c     arguments
      complex ipdot(ncol,ncol,4,nvol)

c     common blocks
      integer forweo,backeo,sindeo,sindeoh,cooreo,parity
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),
     $     parity(nx,ny,nz,nt),cooreo(nvol,4),
     $     forweo(nvol,4),backeo(nvol,4)
      real remu
      real emu,ememu
      common/param4/remu,emu,ememu
      complex phi_e,phi_o,fuf_e,fuf_o
      common/pseudof/phi_e(ncol,nvolh),phi_o(ncol,nvolh),fuf_e(ncol
     $     ,nvolh),fuf_o(ncol,nvolh)
      real pole,coef,cost
      integer nterm
      common/rhmc_a/nterm,cost,pole(nmr),coef(nmr)

c     internal variables
      complex w1(ncol),w2(ncol),w3(ncol),w4(ncol)
      complex v_e(ncol,nvolh,nmr),v_o(ncol,nvolh,nmr)
      integer ie,ivol,icol,idir,iterm,icol1,icol2
      complex chi_t_e(ncol,nvolh,nmr),chi_t_o(ncol,nvolh,nmr)
      complex c1,c3

      call multi_shift_inverter
     $     (chi_t_e,chi_t_o,phi_e,phi_o,
     $     emu,ememu,
     $     pole,
     $     nterm)
      
      
      do iterm=1,nterm
         call D(OE,UNDAG,v_o(1,1,iterm),chi_t_e(1,1,iterm),emu,ememu)
         call D(EO,UNDAG,v_e(1,1,iterm),chi_t_o(1,1,iterm),emu,ememu)
      enddo

      
      do ivol=1,nvol
         ie=ivol-nvolh
         do idir=1,4

            do icol1=1,ncol
               do icol2=1,ncol
                  ipdot(icol1,icol2,idir,ivol)=0
               enddo
            enddo

            do iterm=1,nterm
               do icol=1,ncol
                  if(ivol.le.nvolh) then
                     w1(icol)=v_o(icol,forweo(ivol,idir),iterm)
                     w2(icol)=CONJG(chi_t_e(icol,ivol,iterm))
                     w3(icol)=-chi_t_o(icol,forweo(ivol,idir),iterm) !minus sign 
                     w4(icol)=CONJG(v_e(icol,ivol,iterm)) !! implemented here
                  else
                     w1(icol)=-chi_t_e(icol,forweo(ivol,idir),iterm)!minus sign 
                     w2(icol)=CONJG(v_o(icol,ivol-nvolh,iterm)) !! implemented here
                     w3(icol)=v_e(icol,forweo(ivol,idir),iterm)
                     w4(icol)=CONJG(chi_t_o(icol,ivol-nvolh,iterm))
                  endif
                  if(idir.eq.4) then
                     if(ivol.le.nvolh) then
                        w1(icol)=ememu*w1(icol)
                        w3(icol)=emu*w3(icol)
                     else
                        w1(icol)=emu*w1(icol)
                        w3(icol)=ememu*w3(icol)
                     endif
                  endif
               enddo

               do icol1=1,ncol
                  c1=w1(icol1)
                  c3=w3(icol1)
                  do icol2=1,ncol
                     ipdot(icol1,icol2,idir,ivol)=ipdot(icol1,icol2,idir
     $                    ,ivol)+coef(iterm)*(c1*w2(icol2)+c3*w4(icol2))
                  enddo
               enddo
            enddo

         enddo
      enddo
      
      return
      end
