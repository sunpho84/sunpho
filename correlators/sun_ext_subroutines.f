c=========================================================
      subroutine expmat(u,h)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc exponential of su(n) subroutine up to sixth order
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex u(ncol,ncol),h(ncol,ncol),hp1(ncol,ncol),hp2(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
            u(i,j) = h(i,j)
         enddo
         u(i,i) = u(i,i) + (1.0,0.0)
      enddo

c (1 + x) done

      call mmult(hp1(1,1),h(1,1),h(1,1))

      do i = 1,ncol
         do j = 1,ncol
               u(i,j) = u(i,j) + 0.50000000000*hp1(i,j)
         enddo
      enddo

c second order done

      call mmult(hp2(1,1),hp1(1,1),h(1,1))

      do i = 1,ncol
         do j = 1,ncol
               u(i,j) = u(i,j) + 0.1666666667*hp2(i,j)
         enddo
      enddo

c third order done

      call mmult(hp1(1,1),hp2(1,1),h(1,1))

      do i = 1,ncol
         do j = 1,ncol
               u(i,j) = u(i,j) + 0.0416666667*hp1(i,j)
         enddo
      enddo

c fourth order done

      call mmult(hp2(1,1),hp1(1,1),h(1,1))

      do i = 1,ncol
         do j = 1,ncol
               u(i,j) = u(i,j) + 0.0083333333*hp2(i,j)
         enddo
      enddo

c fifth order done

      call mmult(hp1(1,1),hp2(1,1),h(1,1))

      do i = 1,ncol
         do j = 1,ncol
               u(i,j) = u(i,j) + 0.0013888889*hp1(i,j)
         enddo
      enddo

c sixth order done

      call unitarize1(u(1,1))
   
      return
      end  
c=========================================================

c=========================================================
      subroutine sun_random(u_ran)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc create a random su(n) matrix by multiplying 
cc n(n-1)/2 su(2) random matrices
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      common/param/nran,beta_f,beta_a,rho,pi,iad_flag 
      complex u_l(ncol,ncol),u_ran(ncol,ncol),u_aux(ncol,ncol)

CC    generate i1,i2 indices of SU(2) subgroup
      
      call one(u_ran)
      do i1 = 1,ncol
         do i2 = i1+1,ncol
            call equal(u_aux,u_ran)
            
CC    generate u0,u1,u2,u3 random on the four dim. sphere
            
            u0 = 1.0 - 2.0*ran2()
            alpha = sqrt(1 - u0**2)
            phi = 2.0*pigr*ran2()      
            costheta = 1.0 - 2.0*ran2()      
            sintheta = sqrt(1.0 - costheta**2)
            u3 = alpha*costheta      
            u1 = alpha*sintheta*cos(phi)
            u2 = alpha*sintheta*sin(phi)
            
cc define u_l as unit matrix ...
            call one(u_l)
            
cc ... and then modify the elements in the chosen su(2) subgroup
            u_l(i1,i1) = CMPLX(u0,u3)
            u_l(i1,i2) = CMPLX(u2,u1)
            u_l(i2,i1) = CMPLX(-u2,u1)
            u_l(i2,i2) = CMPLX(u0,-u3)
            
            call mmult(u_ran,u_l,u_aux)
            
         enddo
      enddo

      return
      end

c=========================================================
      subroutine generate_su2rhotables
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc create a table of su2 matrices  in the form of 4
cc real numbers to be used in the metropolis algorithm 
cc the matrices are in a rho neighborhood of the identity 
cc also simple 3-vectors are created to be used in heatbath
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      common/param/nran,beta_f,beta_a,rho,pi,iad_flag
      parameter (MTAB = 1000000) 
      common/su2table/stable(4,MTAB)
      common/su2table_h/stable_h(4,MTAB)

         do itab = 1,MTAB
C     generate u0,u1,u2,u3
            
            u0 = 1.0 - rho*ran2()
            alpha = sqrt(1 - u0**2)
            
            phi = 2.0*pigr*ran2()      
            costheta = 1.0 - 2.0*ran2()      
            sintheta = sqrt(1.0 - costheta**2)
            stable_h(1,itab) = 1.0
            stable_h(2,itab) = sintheta*cos(phi)
            stable_h(3,itab) = sintheta*sin(phi)
            stable_h(4,itab) = costheta      
            u3 = alpha*costheta      
            u1 = alpha*sintheta*cos(phi)
            u2 = alpha*sintheta*sin(phi)
            stable(1,itab) = u0
            stable(2,itab) = u1
            stable(3,itab) = u2
            stable(4,itab) = u3
         enddo

      return
      end
CC=============================================

c=========================================================
      subroutine mmult_su2(ua,i1,i2,a0,a1,a2,a3)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc multiply an sun matrix ua by an su2 matrix extended
cc to sun in the i1,i2 subgroup. the su2 matrix is given
cc in the four-vector form as taken from su2rhotables.
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      common/param/nran,beta_f,beta_a,rho,pi,iad_flag
      parameter (MTAB = 1000000) 
      common/su2table/stable(4,MTAB)
      complex ua(ncol,ncol),col1(ncol),col2(ncol)
      complex c11,c12,c21,c22
      real a0,a1,a2,a3
      integer i1,i2

!! creo gli elementi della matrice su2

      c11 = cmplx(a0,a3)
      c22 = conjg(c11)
      c12 = cmplx(a2,a1)
      c21 = - conjg(c12)

!! creo le due righe della matrice che vengono modificate
      do ic = 1,ncol
         col1(ic) = c11 * ua(i1,ic) + c12 * ua(i2,ic)
         col2(ic) = c21 * ua(i1,ic) + c22 * ua(i2,ic)
      enddo

!! sostituisco le due righe nuove nella matrice originaria
      do ic = 1,ncol
         ua(i1,ic) = col1(ic) 
         ua(i2,ic) = col2(ic) 
      enddo

      return
      end
CC=============================================




