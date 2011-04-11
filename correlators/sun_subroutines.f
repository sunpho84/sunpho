c=========================================================
      subroutine mmult(u,u1,u2)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc multiplication of two su(n) matrices: u = u1 * u2
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex u(ncol,ncol),u1(ncol,ncol),u2(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
            u(i,j) = (0.,0.)
            do l = 1,ncol
               u(i,j) = u(i,j) + u1(i,l)*u2(l,j)
            enddo
         enddo
      enddo
   
      return
      end  
c=========================================================
c=========================================================
      subroutine mmult_add(u,u1,u2)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc multiplication of two su(n) matrices: u = u1 * u2
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex u(ncol,ncol),u1(ncol,ncol),u2(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
c            u(i,j) = (0.,0.)
            do l = 1,ncol
               u(i,j) = u(i,j) + u1(i,l)*u2(l,j)
            enddo
         enddo
      enddo
   
      return
      end  
c=========================================================
c=========================================================
      subroutine hmmult(u,u1,u2)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc multiplication of two su(n) matrices: u = u1~ * u2
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex u(ncol,ncol),u1(ncol,ncol),u2(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
            u(i,j) = (0.,0.)
            do l = 1,ncol
               u(i,j) = u(i,j) + CONJG(u1(l,i))*u2(l,j)
            enddo
         enddo
      enddo
   
      return
      end  
c=========================================================
c=========================================================
      subroutine hmmult_add(u,u1,u2)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc addition of multiplication of two su(n) matrices: u = u1~ * u2
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex u(ncol,ncol),u1(ncol,ncol),u2(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
            do l = 1,ncol
               u(i,j) = u(i,j) + CONJG(u1(l,i))*u2(l,j)
            enddo
         enddo
      enddo
   
      return
      end  
c=========================================================
c=========================================================
      subroutine mhmult(u,u1,u2)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc multiplication of two su(n) matrices: u = u1 * u2~
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex u(ncol,ncol),u1(ncol,ncol),u2(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
            u(i,j) = (0.,0.)
            do l = 1,ncol
               u(i,j) = u(i,j) + u1(i,l)*CONJG(u2(j,l))
            enddo
         enddo
      enddo
   
      return
      end  
c=========================================================

c=========================================================
      subroutine mhmult_add(u,u1,u2)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc addition of multiplication of two su(n) matrices: u = u + u1 * u2~
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex u(ncol,ncol),u1(ncol,ncol),u2(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
            do l = 1,ncol
               u(i,j) = u(i,j) + u1(i,l)*CONJG(u2(j,l))
            enddo
         enddo
      enddo
   
      return
      end  
c=========================================================

c=========================================================
      subroutine madd(u,u1,u2)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc addition of two su(n) matrices: u = u1 + u2 
cc can be called with identical arguments
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex u(ncol,ncol),u1(ncol,ncol),u2(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
               u(i,j) = u1(i,j) + u2(i,j)
         enddo
      enddo
   
      return
      end  
c=========================================================
      

c=========================================================
      subroutine lincomb(u,a,u1,b,u2)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc addition of two su(n) matrices: u = a*u1 + b*u2 
cc can be called with identical arguments
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex u(ncol,ncol),u1(ncol,ncol),u2(ncol,ncol)
      real a,b   
  
      do i = 1,ncol
         do j = 1,ncol
               u(i,j) = a*u1(i,j) + b*u2(i,j)
         enddo
      enddo
   
      return
      end  
c=========================================================
c=========================================================
      subroutine rmult(u,a,u1)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc multiplication of su(n) matrix by real : u = a*u1 
cc can be called with identical arguments
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex u(ncol,ncol),u1(ncol,ncol)
      real a
  
      do i = 1,ncol
         do j = 1,ncol
               u(i,j) = a*u1(i,j)
         enddo
      enddo
   
      return
      end  
c=========================================================
      
c=========================================================
      subroutine cmult(u,a,u1)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc multiplication of su(n) matrix by real : u = a*u1 
cc can be called with identical arguments
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex u(ncol,ncol),u1(ncol,ncol)
      complex a
  
      do i = 1,ncol
         do j = 1,ncol
               u(i,j) = a*u1(i,j)
         enddo
      enddo
   
      return
      end  
c=========================================================
      
c=========================================================
      subroutine unitarize(u)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc project a (ncol,ncol) matrix onto SU(n)
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex d,c(ncol),u(ncol,ncol)
      real rnorm

CC ORTONORMALIZATION ROW BY ROW

      do i = 1,ncol

c compute the scalar product with previous rows
         do j = 1,i-1
            c(j) = (0.,0.)
            do k = 1,ncol
               c(j) = c(j) + CONJG(u(j,k))*u(i,k)
            enddo
         enddo

c ortogonalize with respect to previous rows
         do j = 1,i-1
            do k = 1,ncol
               u(i,k) = u(i,k) - c(j)*u(j,k)
            enddo
         enddo
   
c normalize
         rnorm = 0.0
         do k = 1,ncol
            rnorm = rnorm + Real(u(i,k))**2 + aimag(u(i,k))**2  
         enddo
         rnorm = 1./sqrt(rnorm)  
         do k = 1,ncol
            u(i,k) = rnorm*u(i,k)
         enddo
         
         
      enddo !! close the loop on rows

c compute the determinant
      call det(d,u)
c correct last row to have determinant = 1
      do k = 1,ncol
         u(ncol,k) = conjg(d)*u(ncol,k)
      enddo
              
      return
      end  
c=========================================================

c=========================================================
      subroutine unitarize1(u)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc project a (ncol,ncol) matrix onto SU(n)
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex c(ncol),u(ncol,ncol)
      real rnorm

CC ORTONORMALIZATION ROW BY ROW

      do i = 1,ncol

c compute the scalar product with previous rows
         do j = 1,i-1
            c(j) = (0.,0.)
            do k = 1,ncol
               c(j) = c(j) + CONJG(u(j,k))*u(i,k)
            enddo
         enddo

c ortogonalize with respect to previous rows
         do j = 1,i-1
            do k = 1,ncol
               u(i,k) = u(i,k) - c(j)*u(j,k)
            enddo
         enddo
   
c normalize
         rnorm = 0.0
         do k = 1,ncol
            rnorm = rnorm + Real(u(i,k))**2 + aimag(u(i,k))**2  
         enddo
         rnorm = 1./sqrt(rnorm)  
         do k = 1,ncol
            u(i,k) = rnorm*u(i,k)
         enddo
         
         
      enddo !! close the loop on rows

c compute the determinant
c      call det(d,u)
c correct last row to have determinant = 1
c      do k = 1,ncol
c         u(ncol,k) = conjg(d)*u(ncol,k)
c      enddo
              
      return
      end  
c=========================================================

c=========================================================
      subroutine det(d,u)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc compute the determinant of an (ncol,ncol) matrix
cc using the sum over permutations 
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      integer perm
      common/perm/ncfact,perm(nperm,ncol),sign(nperm)
      complex u(ncol,ncol),d,c

      d = cmplx(0.,0.)
      
      do i = 1,ncfact
        c = sign(i)
        do j = 1,ncol
           c = c * u(j,perm(i,j))
        enddo
        d = d + c
      enddo   
      
      return
      end  
c=========================================================

c=========================================================
      subroutine one(u)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc define the unit matrix
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex u(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
            u(i,j) = cmplx(0.0,0.0)
         enddo
         u(i,i) = (1.0,0.0)
      enddo

      return
      end

c=========================================================
      subroutine zero(u)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc define the zero matrix
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex u(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
            u(i,j) = cmplx(0.0,0.0)
         enddo
      enddo

      return
      end
c=========================================================
      subroutine equal(ua,ub)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc ua = ub for su(N) matrices
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex ua(ncol,ncol),ub(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
            ua(i,j) = ub(i,j) 
         enddo
      enddo

      return
      end

c=========================================================
      subroutine equalh(ua,ub)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc ua = ub~ for su(N) matrices
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex ua(ncol,ncol),ub(ncol,ncol)

      do i = 1,ncol
         do j = 1,ncol
            ua(i,j) = CONJG(ub(j,i)) 
         enddo
      enddo

      return
      end



c=========================================================
      subroutine ctrace(trace,ua)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc trace = complex trace of ua 
c=========================================================

      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex ua(ncol,ncol)
      complex trace

      trace = (0.,0.)
      do i = 1, ncol
         trace = trace + ua(i,i)
      enddo

      return
      end
c=========================================================
      subroutine rtrace(trace,ua)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc trace = real trace of ua
c=========================================================

      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex ua(ncol,ncol)
      real trace

      trace = 0.
      do i = 1, ncol
         trace = trace + Real(ua(i,i))
      enddo

      return
      end

c=========================================================
      subroutine ctrace_uu(ctrace,ua,ub)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc ctrace = complex trace of     ua * ub
c=========================================================

      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex ua(ncol,ncol),ub(ncol,ncol)
      complex ctrace
      integer i,j

      ctrace = (0.,0.)
      do i = 1, ncol
         do j = 1,ncol
            ctrace = ctrace + ua(i,j)*ub(j,i)
         enddo
      enddo

      return
      end
c=========================================================
      subroutine ctrace_uuh(ctrace,ua,ub)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc ctrace = complex trace of     ua * ub~
c=========================================================

      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex ua(ncol,ncol),ub(ncol,ncol)
      complex ctrace
      integer i,j

      ctrace = (0.,0.)
      do i = 1, ncol
         do j = 1,ncol
            ctrace = ctrace + ua(i,j)*CONJG(ub(i,j))
         enddo
      enddo

      return
      end
c=========================================================
      subroutine rtrace_uu(rtrace,ua,ub)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc ctrace = real trace of     ua * ub~
c=========================================================

      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex ua(ncol,ncol),ub(ncol,ncol)
      complex ctrace
      real rtrace

      ctrace = (0.,0.)
      do i = 1, ncol
         do j = 1,ncol
            ctrace = ctrace + ua(i,j)*ub(j,i)
         enddo
      enddo
      rtrace = real(ctrace)
      return
      end

c=========================================================
      subroutine rtrace_uuh(rtrace,ua,ub)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc ctrace = real trace of     ua * ub~
c=========================================================

      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex ua(ncol,ncol),ub(ncol,ncol)
      complex ctrace
      real rtrace

      ctrace = (0.,0.)
      do i = 1, ncol
         do j = 1,ncol
            ctrace = ctrace + ua(i,j)*CONJG(ub(i,j))
         enddo
      enddo
      rtrace = real(ctrace)
      return
      end

c=========================================================
      subroutine TA(ua,ub)
c  written by Massimo D'Elia
c  version 1.0 - 01/2008
cc ua is the traceless antihermitean part of ub
c=========================================================

      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"

c     argomenti
      complex ua(ncol,ncol),ub(ncol,ncol)

c     variabili interne
      complex ctrace

c     parametri
      parameter(xinvncol=1.0/ncol)


      ctrace=cmplx(0,0)
      do i=1,ncol
         do j=1,ncol
            ua(i,j)=0.5*(ub(i,j)-conjg(ub(j,i)))
         enddo
         ctrace=ctrace+ua(i,i)
      enddo
      do i=1,ncol
         ua(i,i)=ua(i,i)-xinvncol*ctrace
      enddo

      return
      end

