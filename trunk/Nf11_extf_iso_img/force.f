c=========================================================
      subroutine compute_utpdt(iepsatt)
c  written by Massimo D'Elia
c  version 1.0 - 2007/2008
c  evolve all the links for a "iepsatt" MD time
c=========================================================
      implicit none
#include "parameters.f"

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
            call mmult(u2,u1,u(1,1,idir,ivol))
            call equal(u(1,1,idir,ivol),u2)
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
#include "parameters.f"

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
#include "parameters.f"

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

#ifdef debug1
      write(*,*) "Starting calculation of the force"
#endif
      
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

#ifdef debug1
      write(*,*) "Force calculation terminated"
#endif
      
#ifdef debug2
      write(*,*) "FORCE:",ipdot_g(1,1,1,1), ipdot_q(1,1,1,1)
#endif      
      
      return
      end


c=========================================================
      subroutine compute_gluonic_force(ipdot)
c  written by Sanfo
c  calculte gluonic part of the force
c=========================================================
      implicit none
#include "parameters.f"

c     arguments
      complex ipdot(ncol,ncol,4,nvol)

c     common blocks
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue
      complex staple
      common/staple/staple(ncol,ncol,4,nvol)

c     internal variables
      integer ivol,idir
      real r_1

      r_1=beta/float(ncol) 

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
#include "parameters.f"

c     arguments
      complex ipdot(ncol,ncol,4,nvol)

c     common blocks
      integer forweo,backeo,sindeo,sindeoh,cooreo,parity
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),
     $     parity(nx,ny,nz,nt),cooreo(nvol,4),
     $     forweo(nvol,4),backeo(nvol,4)
      complex potc1,potc2
      common/param3/potc1,potc2
      complex phi1,phi2,fuf1,fuf2
      common/pseudof1/phi1(ncol,nvol_eo),fuf1(ncol,nvol_eo)
      common/pseudof2/phi2(ncol,nvol_eo),fuf2(ncol,nvol_eo)
      real pole,coef,cost
      integer nterm
      common/rhmc_a/nterm,cost,pole(nmr),coef(nmr)
#ifdef mf
      complex ue1,ue2
      common/fielde1/ue1(4,nvol)
      common/fielde2/ue2(4,nvol)
#endif

c     internal variables
      integer ie,ivol,icol,idir,iterm,icol1,icol2
      complex chi1_t(ncol,nvol_eo,nmr),chi2_t(ncol,nvol_eo,nmr)
      complex w11(ncol),w12(ncol)
      complex w21(ncol),w22(ncol)
      complex v1_o(ncol,nvolh,nmr),v2_o(ncol,nvolh,nmr)
      complex c1,c2
      complex exp_1,exp_2
#ifdef ficp
      complex w13(ncol),w14(ncol)
      complex w23(ncol),w24(ncol)
      complex v1_e(ncol,nvolh,nmr),v2_e(ncol,nvolh,nmr)
      complex c3,c4
      complex exp_m_c_1,exp_m_c_2
#endif

#ifdef mf
      complex c_extf1,c_extf2
#endif

      call multi_shift_inverter
     $     (chi1_t,phi1,
     $     potc1,
     $     pole,
     $     nterm,1)

      call multi_shift_inverter
     $     (chi2_t,phi2,
     $     potc2,
     $     pole,
     $     nterm,2)


#ifdef mf
      call add_extf(1)
#endif

      exp_1=exp(+potc1)
      exp_2=exp(+potc2)
#ifdef ficp
      exp_m_c_1=exp(-conjg(potc1))
      exp_m_c_2=exp(-conjg(potc2))
#endif

      do iterm=1,nterm


!     calculation of force, without or with EO improvement
         call D(OE,v1_o(1,1,iterm),chi1_t(1,1,iterm),potc1)
#ifdef ficp
         call D(EO,v1_e(1,1,iterm),chi1_t(1,1+nvolh,iterm),potc1)
#endif

      enddo


#ifdef mf
      call rem_extf(1)
      call add_extf(2)
#endif

      do iterm=1,nterm


!     calculation of force, without or with EO improvement
         call D(OE,v2_o(1,1,iterm),chi2_t(1,1,iterm),potc2)
#ifdef ficp
         call D(EO,v2_e(1,1,iterm),chi2_t(1,1+nvolh,iterm),potc2)
#endif


      enddo

#ifdef mf
      call rem_extf(2)
#endif

      do ivol=1,nvol
         ie=ivol-nvolh
         do idir=1,4

#ifdef mf
            c_extf1=ue1(idir,ivol)
            c_extf2=ue2(idir,ivol)
#endif

            do icol1=1,ncol
               do icol2=1,ncol
                  ipdot(icol1,icol2,idir,ivol)=0
               enddo
            enddo

            do iterm=1,nterm
               do icol=1,ncol
                  if(ivol.le.nvolh) then
                     w11(icol)=v1_o(icol,forweo(ivol,idir),iterm)
                     w21(icol)=v2_o(icol,forweo(ivol,idir),iterm)
                     w12(icol)=CONJG(chi1_t(icol,ivol,iterm))
                     w22(icol)=CONJG(chi2_t(icol,ivol,iterm))

#ifdef ficp
                     w13(icol)=-chi1_t(icol,forweo(ivol,idir)+nvolh
     $                    ,iterm)
                     w23(icol)=-chi2_t(icol,forweo(ivol,idir)+nvolh
     $                    ,iterm)
                     w14(icol)=CONJG(v1_e(icol,ivol,iterm))
                     w24(icol)=CONJG(v2_e(icol,ivol,iterm))
#endif

                  else
                     w11(icol)=-chi1_t(icol,forweo(ivol,idir),iterm)!minus sign 
                     w21(icol)=-chi2_t(icol,forweo(ivol,idir),iterm)
                     w12(icol)=CONJG(v1_o(icol,ie,iterm)) !! implemented here
                     w22(icol)=CONJG(v2_o(icol,ie,iterm))

#ifdef ficp
                     w13(icol)=v1_e(icol,forweo(ivol,idir),iterm)
                     w23(icol)=v2_e(icol,forweo(ivol,idir),iterm)
                     w14(icol)=CONJG(chi1_t(icol,ivol,iterm))
                     w24(icol)=CONJG(chi2_t(icol,ivol,iterm))
#endif

                  endif
               
#ifndef isotropic
                  if(idir.eq.4) then
#ifdef ficp
c     note that if re(potc)=0, -conjg(potc)=potc
                     if(ivol.le.nvolh) then
                        w11(icol)=exp_m_c_1*w11(icol)
                        w21(icol)=exp_m_c_2*w21(icol)
                        w13(icol)=exp_1*w13(icol)
                        w23(icol)=exp_2*w23(icol)
                     else
#endif
                        w11(icol)=exp_1*w11(icol)
                        w21(icol)=exp_2*w21(icol)
#ifdef ficp
                        w13(icol)=exp_m_c_1*w13(icol)
                        w23(icol)=exp_m_c_2*w23(icol)
                     endif
#endif
                  endif
#endif
               enddo

               do icol1=1,ncol

#ifdef mf                  
                  c1=w11(icol1)*c_extf1
                  c2=w21(icol1)*c_extf2

#ifdef ficp
                  c3=w13(icol1)*c_extf1
                  c4=w23(icol1)*c_extf2
#endif

#else
                  c1=w11(icol1)
                  c2=w21(icol1)

#ifdef ficp
                  c3=w13(icol1)
                  c4=w23(icol1)
#endif

#endif
                  do icol2=1,ncol
                     ipdot(icol1,icol2,idir,ivol)=ipdot(icol1,icol2,idir
     $                   ,ivol)+coef(iterm)*(c1*w12(icol2)+c2*w22(icol2)
#ifdef ficp
     $                    +c3*w14(icol2)+c4*w24(icol2)
#endif

     $                    )
                  enddo
               enddo
            enddo
            
         enddo
      enddo
      
      return
      end
