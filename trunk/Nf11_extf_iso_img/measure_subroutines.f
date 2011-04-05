c=========================================================
      subroutine measure(acc,n_traj,n_rand)
c  written by Sanfo
c  version 1.0 - 04/2008
c=========================================================
      implicit none
#include "parameters.f"

c     arguments
      integer acc
      integer n_traj
      integer n_rand

c     common blocks
      complex potc1,potc2
      common/param3/potc1,potc2

c     internal variables
      real plaq_s,plaq_t,plaq2_s,plaq2_t
      complex ploop
#ifdef chiral_meas
      complex f1,f2,f3,f4,f5,f6,f7,fv1,fv2,fv3,fv4,fv5,fv6,fv7
      complex g1,g2,g3,g4,g5,g6,g7,gv1,gv2,gv3,gv4,gv5,gv6,gv7
#endif

#ifdef mf
      complex m1
#endif

#ifdef debug1
      real timei,timef

      write(*,*)
      write(*,*)
      write(*,*) "Init measures on config: ",n_traj
      call cpu_time(timei)
#endif
      n_rand=n_rand

c     measure section
      call plaquette(plaq_s,plaq_t,plaq2_s,plaq2_t) 
      call polyakov(ploop)

#ifdef chiral_meas
      call CHIRAL(n_rand,f1,f2,f3,f4,f5,f6,f7,fv1,fv2,fv3,fv4,fv5,fv6
     $     ,fv7,potc1,1)
      call CHIRAL(n_rand,g1,g2,g3,g4,g5,g6,g7,gv1,gv2,gv3,gv4,gv5,gv6
     $     ,gv7,potc2,2)
#endif

#ifdef mf
      call MAGNET(n_rand,m1)
#endif
      
c     the file writing is collapsed in order to minimize possible failure
      write(u_meas,*) n_traj,acc,plaq_s,plaq_t,real(ploop),aimag(ploop)

#ifdef chiral_meas
      write(u_ferm_11,*) n_traj,Real(f1),Aimag(f1),Real(fv1),Aimag(fv1)
      write(u_ferm_21,*) n_traj,Real(g1),Aimag(g1),Real(gv1),Aimag(gv1)
      write(u_ferm_12,*) n_traj,Real(f2),Aimag(f2),Real(fv2),Aimag(fv2)
      write(u_ferm_22,*) n_traj,Real(g2),Aimag(g2),Real(gv2),Aimag(gv2)
      write(u_ferm_13,*) n_traj,Real(f3),Aimag(f3),Real(fv3),Aimag(fv3)
      write(u_ferm_23,*) n_traj,Real(g3),Aimag(g3),Real(gv3),Aimag(gv3)
      write(u_ferm_14,*) n_traj,Real(f4),Aimag(f4),Real(fv4),Aimag(fv4)
      write(u_ferm_24,*) n_traj,Real(g4),Aimag(g4),Real(gv4),Aimag(gv4)
#endif

#ifdef mf
      write(u_ferm_15,*) n_traj,Real(f5),imag(f5),Real(fv5),imag(fv5)
      write(u_ferm_25,*) n_traj,Real(g5),imag(g5),Real(gv5),imag(gv5)
      write(u_ferm_16,*) n_traj,Real(f6),imag(f6),Real(fv6),imag(fv6)
      write(u_ferm_26,*) n_traj,Real(g6),imag(g6),Real(gv6),imag(gv6)
      write(u_ferm_17,*) n_traj,Real(f7),imag(f7),Real(fv7),imag(fv7)
      write(u_ferm_27,*) n_traj,Real(g7),imag(g7),Real(gv7),imag(gv7)
      write(u_magnet,*) n_traj,Real(m1),Aimag(m1)
#endif

#ifdef debug1
      call cpu_time(timef)
      write(*,*) "End measures, lasted time: ",int(timef-timei)
     $     ,"seconds"
      write(*,*)
      write(*,*)
#endif

      return
      end


c=========================================================
      subroutine local_plaq(locplaq)
c     written by Massimo D'Elia
c     version 1.0 - 01/08

c     calculate the plaquette site by site. Used to avoid numerical
c     rounding during the calculation of the diference between new and
c     old gauge action (so do not change even if it seems unoptimized)
c=========================================================
      implicit none
#include "parameters.f"

c     arguments
      real locplaq(nvol)

c     common blocks
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue
      integer forweo2,backeo2
      common/geo2/forweo2(nvol,4),backeo2(nvol,4)
      complex u
      common/field/u(ncol,ncol,4,nvol)

c     internal variables
      integer i,mu,nu
      complex v1(ncol,ncol),v2(ncol,ncol)
      complex ctr
      real ctr_tot

c     parameters
      real xinvncol
      parameter(xinvncol=1.0/ncol)


      do i=1,nvol
         
         ctr_tot=0
         
         do mu=1,3

            do nu=mu+1,4        !loop on orthogonal directions

               call mmult(v1,u(1,1,nu,i),u(1,1,mu,forweo2(i,nu)))
               call mmult(v2,u(1,1,mu,i),u(1,1,nu,forweo2(i,mu)))
               call ctrace_uuh(ctr,v1,v2)

               ctr_tot=ctr_tot+real(ctr)

            enddo               !nu
            
            locplaq(i)=6+xinvncol*ctr_tot !include segno fasi stag

         enddo                  !! mu

      enddo                     !! i

      return
      end



c=========================================================
      subroutine compute_staples
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc    compute staples to be used in the updating
c=========================================================
      implicit none
#include "parameters.f"

c     common blocks
      integer sind,coor,forw,back
      common/ge/sind(nx,ny,nz,nt),coor(nvol,4),forw(nvol,4),back(nvol,4)
      integer forweo2,backeo2
      common/geo2/forweo2(nvol,4),backeo2(nvol,4)
      complex u
      common/field/u(ncol,ncol,4,nvol)
      complex staple
      common/staple/staple(ncol,ncol,4,nvol)

c     internal variables
      integer ind,idir,jdir
      complex v1(ncol,ncol)

      do ind=1,nvol
         
         do idir=1,4
            
            call zero(staple(1,1,idir,ind))
c     loop on orthogonal directions
            do jdir=1,4

               if(jdir.ne.idir) then
                  
c     staple in the forward direction
                  call mmult(v1,u(1,1,jdir,ind),u(1,1,idir,forweo2(ind
     $                 ,jdir)))
                  call mhmult_add(staple(1,1,idir,ind),u(1,1,jdir
     $                 ,forweo2(ind,idir)),v1)
c     staple in the backward direction
                  call mmult(v1,u(1,1,idir,backeo2(ind,jdir)),u(1,1
     $                 ,jdir,forweo2(backeo2(ind,jdir),idir)))
                  call hmmult_add(staple(1,1,idir,ind),v1,u(1,1,jdir
     $                 ,backeo2(ind,jdir)))
               
               endif

            enddo
            
         enddo                  !! idir
      enddo                     !! ivol

      return
      end

c=========================================================
      subroutine plaquette(plaq_s,plaq_t,plaq2_s,plaq2_t)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
c  measures fundamental and adjoint plaquette
c=========================================================

      implicit none
#include "parameters.f"

c     arguments
      real plaq_s,plaq_t,plaq2_s,plaq2_t

c     common blocks
      integer forweo2,backeo2
      common/geo2/forweo2(nvol,4),backeo2(nvol,4)
      complex u
      common/field/u(ncol,ncol,4,nvol)

c     local variables
      integer ind,idir,jdir
      complex v1(ncol,ncol),v2(ncol,ncol)
      complex ctr


      ind=0
      plaq_s=0.0                !! spatial and temporal plaquette
      plaq_t=0.0
      plaq2_s=0.0               !! spatial and temporal plaquette squared 
      plaq2_t=0.0

      do ind =1,nvol
         
         do idir=1,3

c     loop on orthogonal directions

            do jdir=idir+1,4
  
               call mmult(v1,u(1,1,jdir,ind),u(1,1,idir,forweo2(ind
     $              ,jdir)))
               call mmult(v2,u(1,1,idir,ind),u(1,1,jdir,forweo2(ind
     $              ,idir)))
               call ctrace_uuh(ctr,v1,v2)

               if (jdir.eq.4) then
                  plaq_t= plaq_t+Real(ctr)
                  plaq2_t= plaq2_t+Real(ctr)**2+aimag(ctr)**2
               else
                  plaq_s=plaq_s+Real(ctr)
                  plaq2_s=plaq2_s+Real(ctr)**2+aimag(ctr)**2
               endif

            enddo               !! jdir: closing do over orthogonal directions

         enddo                  !! idir: closing do over directions

      enddo                     !! ind

      plaq_s=-plaq_s/(ncol*3.0*nvol)     !! minus sign due to 
      plaq_t=-plaq_t/(ncol*3.0*nvol)     !! staggered phases
      plaq2_s=plaq2_s/(ncol*ncol*3.0*nvol)
      plaq2_t=plaq2_t/(ncol*ncol*3.0*nvol)

      return
      end

c=========================================================
      subroutine polyakov(poly_c)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
c=========================================================

      implicit none
#include "parameters.f"

c     arguments
      complex poly_c

c     common blocks
      integer sindeo,sindeoh,parity,cooreo,forweo,backeo
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),
     $  parity(nx,ny,nz,nt),cooreo(nvol,4),forweo(nvol,4),backeo(nvol,4)
      complex u
      common/field/u(ncol,ncol,4,nvol)

c     internal variables
      integer ix,iy,iz,it,idir
      complex v1(ncol,ncol),poly_mat(ncol,ncol)
      complex ctr


      poly_c=cmplx(0.0,0.0)

      do iz=1,nz
         do iy=1,ny
            do ix=1,nx

               call one(poly_mat)
               idir=4
               
               do it=1,nt
                  call mmult(v1,poly_mat,u(1,1,idir,sindeo(ix,iy,iz
     $                 ,it)))
                  call equal(poly_mat,v1)
               enddo
               
               call ctrace(ctr,poly_mat)
               poly_c=poly_c+ctr 
               
            enddo               !! ix
         enddo                  !! iy
      enddo                     !! iz

c     minus due to staggered phase
      poly_c=-poly_c/float(ncol*nx*ny*nz)

      return
      end


#ifdef chiral_meas
c=========================================================
      subroutine CHIRAL(n_rand,c1,c2,c3,c4,c5,c6,c7,
     $     cv1,cv2,cv3,cv4,cv5,cv6,cv7,potc,ud)
c  written by Massimo D'Elia
c  version 1.0 - 2007/2008
c  c1 -> chiral condensate
c  c2 -> energy density
c  c3 -> baryon number
c  c4 -> pressure
c  c5 -> baryon current direction x
c  c6 -> baryon current direction y
c  c7 -> baryon current direction z
c=========================================================

      implicit none
#include "parameters.f"
      
c     arguments
      integer n_rand
      complex c1,c2,c3,c4,c5,c6,c7,cv1,cv2,cv3,cv4,cv5,cv6,cv7
      complex potc
      integer ud

c     common blocks
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue
      integer sindeo,sindeoh,parity,cooreo,forweo,backeo
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),
     $     parity(nx,ny,nz,nt),cooreo(nvol,4),
     $     forweo(nvol,4),backeo(nvol,4)
      complex u
      common/field/u(ncol,ncol,4,nvol)
      integer nterm
      real pole,coef,cost
      common/rhmc_a/nterm,cost,pole(nmr),coef(nmr)

c     internal variables
      integer iaux,i_rand,ivol,idir,icol,iodd
      complex v1(ncol),v2(ncol)
      complex phi(ncol,nvol)
      complex chi(ncol,nvol)
      complex rnd(ncol,nvol)
      complex b1,b2,b3,b4,a1,a2
      complex ac1(3),ac2(3)
      complex bc(3)

c     parameters
      real sigma,nfqv
      parameter(sigma=1,nfqv=n_quark_for_flavour/4.0/(nvol))


#ifdef mf
      call add_extf(ud)
#else
      ud=ud !to avoid warning
#endif

      c1=(0,0)
      c2=(0,0)
      c3=(0,0)
      c4=(0,0)
      c5=(0,0)
      c6=(0,0)
      c7=(0,0)
      cv1=(0,0)
      cv2=(0,0)
      cv3=(0,0)
      cv4=(0,0)
      cv5=(0,0)
      cv6=(0,0)
      cv7=(0,0)

      do i_rand=1,n_rand  !! loop over random vectors

#ifdef debug1
         write(*,*) "Calculation of chiral observables, vector:",i_rand
#endif

         b1=(0,0)
         b2=(0,0)
         b3=(0,0)
         b4=(0,0)
         a1=(0,0)
         a2=(0,0)
         do iaux=1,3
            bc(iaux)=(0,0)
            ac1(iaux)=(0,0)
            ac2(iaux)=(0,0)
         enddo   
         
         do ivol = 1,nvol
            call gauss_vector(rnd(1,ivol),sigma)
         enddo 

#ifdef ficp
         call D(OE_DAG,phi(1,1),rnd(1,1+nvolh),potc)
         call D(EO_DAG,phi(1,1+nvolh),rnd(1,1),potc)
#else
         call D(EO,phi(1,1),rnd(1,1+nvolh),potc)
#endif

         do ivol = 1,nvol
            do icol = 1,ncol
               phi(icol,ivol)=mass*rnd(icol,ivol)+phi(icol,ivol)
            enddo
         enddo 

         call singol_inverter(chi,phi,potc,mass2)
#ifndef ficp
         call D(OE,phi(1,1+nvolh),chi(1,1),potc)
#endif

!!=========COMPUTE CHIRAL CONDENSATE=====================================      
         do ivol=1,nvol_eo
            do icol=1,ncol
               b1=b1+conjg(rnd(icol,ivol))*chi(icol,ivol)
            enddo
         enddo 

!!=========COMPUTE ENERGY AND BARYON DENSITY=============================      
         do ivol = 1,nvolh
            iodd = ivol + nvolh
            
            call vmult(v1,u(1,1,4,ivol),chi(1,forweo(ivol,4)+nvolh))
            call vmult(v2,u(1,1,4,ivol),rnd(1,forweo(ivol,4)+nvolh))
c     v1 should be multiplied by eim -> a1 by eim
c     v2 should be multiplied by eim -> a2 by emim
            do icol = 1,ncol
               a1=a1+conjg(rnd(icol,ivol))*v1(icol)
               a2=a2+conjg(v2(icol))*chi(icol,ivol)
            enddo
            
            call vmult(v1,u(1,1,4,iodd),chi(1,forweo(iodd,4)))
            call vmult(v2,u(1,1,4,iodd),rnd(1,forweo(iodd,4)))
c     v1 should be multiplied by eim -> a1 by eim
c     v2 should be multiplied by eim -> a2 by emim
            do icol = 1,ncol
               a1=a1+conjg(rnd(icol,iodd))*v1(icol)
               a2=a2+conjg(v2(icol))*chi(icol,iodd)
            enddo
            
         enddo 
         
         b2=exp(potc)*a1-exp(-potc)*a2  !! energy density
         b3=exp(potc)*a1+exp(-potc)*a2  !! baryon density
         
!!=========COMPUTE PRESSURE DENSITY=====================================      
         do ivol=1,nvolh            
            iodd=ivol+nvolh !! pressure 
            
            do idir=1,3
               
               call vmult(v1,u(1,1,idir,ivol),chi(1,forweo(ivol,idir)
     $              +nvolh))
               call vmult(v2,u(1,1,idir,ivol),rnd(1,forweo(ivol,idir)
     $              +nvolh))
              do icol = 1,ncol
                 ac1(idir)=ac1(idir)+CONJG(rnd(icol,ivol))*v1(icol)
                 ac2(idir)=ac2(idir)+CONJG(v2(icol))*chi(icol,ivol)
              enddo
              
              call vmult(v1,u(1,1,idir,iodd),chi(1,forweo(iodd,idir)))
              call vmult(v2,u(1,1,idir,iodd),rnd(1,forweo(iodd,idir)))
              do icol = 1,ncol
                 ac1(idir)=ac1(idir)+CONJG(rnd(icol,iodd))*v1(icol)
                 ac2(idir)=ac2(idir)+CONJG(v2(icol))*chi(icol,iodd)
              enddo
              
           enddo                !! spatial direction
           
        enddo 
        
        do iaux=1,3
           bc(iaux)=ac1(iaux)+ac2(iaux) 
           b4=b4+ac1(iaux)-ac2(iaux) 
        enddo

        b1=b1*nfqv
        b2=0.5*b2*nfqv
        b3=0.5*b3*nfqv
        b4=0.5*b4*nfqv

        do iaux=1,3
           bc(iaux)=0.5*bc(iaux)/float(nvol)
        enddo
       
        c1=c1+b1
        c2=c2+b2
        c3=c3+b3
        c4=c4+b4
        c5=c5+bc(1)
        c6=c6+bc(2)
        c7=c7+bc(3)
        cv1=cv1+cmplx(real(b1)**2,aimag(b1)**2)
        cv2=cv2+cmplx(real(b2)**2,aimag(b2)**2)
        cv3=cv3+cmplx(real(b3)**2,aimag(b3)**2)
        cv4=cv4+cmplx(real(b4)**2,aimag(b4)**2)
        cv5=cv5+cmplx(real(bc(1))**2,aimag(bc(1))**2)
        cv6=cv6+cmplx(real(bc(2))**2,aimag(bc(2))**2)
        cv7=cv7+cmplx(real(bc(3))**2,aimag(bc(3))**2)

      enddo                     !! random vectors
      
      c1=c1/float(n_rand)
      c2=c2/float(n_rand)
      c3=c3/float(n_rand)
      c4=c4/float(n_rand)
      c5=c5/float(n_rand)
      c6=c6/float(n_rand)
      c7=c7/float(n_rand)
      cv1=cv1/float(n_rand)
      cv2=cv2/float(n_rand)
      cv3=cv3/float(n_rand)
      cv4=cv4/float(n_rand)
      cv5=cv5/float(n_rand)
      cv6=cv6/float(n_rand)
      cv7=cv7/float(n_rand)

#ifdef mf
      call rem_extf(ud)
#endif

      return
      end
#endif


#ifdef mf
c=========================================================
      subroutine magnet_flav(rnd,b3,potc,ud)
c     written by Massimo D'Elia
c     version 1.0 - 2007/2008
c=========================================================

      implicit none
#include "parameters.f"

c     arguments
      complex rnd(ncol,nvol)
      integer ud
      complex b3,potc

c     variabili passate da blocchi common
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue
      integer sindeo,sindeoh,parity,cooreo,forweo,backeo
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),parity(nx,ny
     $     ,nz,nt),cooreo(nvol,4),forweo(nvol,4),backeo(nvol,4)
      complex u
      common/field/u(ncol,ncol,4,nvol)
      real u1phase
      common/expu1/u1phase(4,nvol)

c     variabili interne
      integer ivol,iodd,icol,idir
      complex phi(ncol,nvol)
      complex chi(ncol,nvol)
      complex v1(ncol),v2(ncol)
      complex t3
      

      call add_extf(ud)

      
#ifdef ficp
      call D(OE_DAG,phi,rnd(1,1+nvolh),potc)
      call D(EO_DAG,phi(1,1+nvolh),rnd,potc)
#else
      call D(EO,phi(1,1),rnd(1,1+nvolh),potc)
#endif
      
      do ivol=1,nvol_eo
         do icol=1,ncol
            phi(icol,ivol)=mass*rnd(icol,ivol)+phi(icol,ivol)
         enddo
      enddo 
      

      call singol_inverter(chi,phi,potc,mass2)
#ifndef ficp
         call D(OE,phi(1,1+nvolh),chi(1,1),potc)
#endif
      
!     =========================COMPUTE MAGN_SUSC=========================
      do ivol=1,nvolh

         iodd=ivol+nvolh
         
         do idir=1,4
            call vmult(v1,u(1,1,idir,ivol),chi(1,forweo(ivol,idir)
     $           +nvolh))
            call vmult(v2,u(1,1,idir,ivol),rnd(1,forweo(ivol,idir)
     $           +nvolh))

c     v1,v2 should be multiplied by +,-u1phase -> t3 by +u1phase
            t3=(0,0)
            do icol=1,ncol
               t3=t3+CONJG(rnd(icol,ivol))*v1(icol)
               t3=t3+CONJG(v2(icol))*chi(icol,ivol)
            enddo
            b3=b3+t3*u1phase(idir,ivol)

            call vmult(v1,u(1,1,idir,iodd),chi(1,forweo(iodd,idir)))
            call vmult(v2,u(1,1,idir,iodd),rnd(1,forweo(iodd,idir)))

c     v1,v2 should be multiplied by +,-u1phase -> t3 by +u1phase
            t3=(0,0)
            do icol=1,ncol
               t3=t3+CONJG(rnd(icol,iodd))*v1(icol)
               t3=t3+CONJG(v2(icol))*chi(icol,iodd)
            enddo
            b3=b3+t3*u1phase(idir,iodd)

         enddo

      enddo 
      
      call rem_extf(ud)

      
      return
      end

c=========================================================
      subroutine magnet(n_rand,c3)
c     written by Massimo D'Elia
c     version 1.0 - 2007/2008
c     calculate the derivative of F/T=lnZ with respect to n=B/B0, where
c     B0=2pi/Nx or 2pi/Nx^2 according to the quantization condition. The
c     result is divided by Nvol for historical reason
c=========================================================

      implicit none
#include "parameters.f"

c     arguments
      integer n_rand
      complex c3

c     variabili interne
      complex c1,c2
      integer i_rand,ivol
      complex rnd(ncol,nvol)

c     variabili passate da blocchi common
      complex potc1,potc2
      common/param3/potc1,potc2

c     parametri
      real nfqv
      parameter(nfqv=n_quark_for_flavour/(4.0*nvol))
      real sigma
      parameter(sigma=1)


      c1=(0,0)
      c2=(0,0)

      do i_rand=1,n_rand        !! loop over random vectors
         
#ifdef debug1
         write(*,*) "Calculation of magnetic observables, vector:"
     $        ,i_rand
#endif

         do ivol=1,nvol
            call gauss_vector(rnd(1,ivol),sigma)
         enddo 
         
         call magnet_flav(rnd,c1,potc1,1)
         call magnet_flav(rnd,c2,potc2,2)
         
      enddo

      c3=(2*c1-c2)*nfqv/float(n_rand)
      
      
      return
      end
#endif
