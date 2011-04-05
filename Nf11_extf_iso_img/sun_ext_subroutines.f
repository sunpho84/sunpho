c=========================================================
      subroutine expmat(u,h)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc exponential of su(n) subroutine up to sixth order
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
#include "parameters.f"
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
      implicit none
#include "parameters.f"
      complex u_l(ncol,ncol),u_ran(ncol,ncol),u_aux(ncol,ncol)

      integer i1,i2
      real u0,alpha,phi,ran2,sintheta,costheta,u1,u2,u3

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
