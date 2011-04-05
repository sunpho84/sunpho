c =========================================================================
      subroutine gauss_matrix(h,sigma)
c =========================================================================
c IT'S GOOD ONLY FOR SU(3) !!!
!!===========================================================================!!
!! A gauss_matrix is a linear combination of eight LAMBDA matrices. The      !!
!! coefficients are real numbers with gaussian distribution. We have         !!
!!          f(z) = 1/(sigma pi) * exp[ - z z* / sigma]  ,  z complex,        !!
!! so that the integral f(z)dz dz* is one.                                   !!
!! Further  we have < z z* > = sigma .                                       !!
!!                                                                           !!
!! We first create 8 gaussian numbers                                        !!
!! with sigma=1 and then we use the eight real components of this vector     !!
!! to get a linear combination of the eight LAMDA matrices. Note, that the   !!
!! expectation value of each coefficient is 1/2 because we created complex   !!
!! numbers with expectation value 1.                                         !! 
!!                                                                           !!
!! Our choice of the LAMBDA matrices is so, that the trace of the product    !!
!! of two LAMBDAS A and B is 2 if A = B and 0 else:                          !!
!!                                                                           !!
!!           |  0  1  0 |             |  0 -i  0 |             |  1  0  0 |  !!
!! lambda1 = |  1  0  0 | , lambda2 = |  i  0  0 | , lambda3 = |  0 -1  0 |  !!
!!           |  0  0  0 |             |  0  0  0 |             |  0  0  0 |  !!
!!                                                                           !!
!!           |  0  0  1 |             |  0  0 -i |             |  0  0  0 |  !!
!! lambda4 = |  0  0  0 | , lambda5 = |  0  0  0 | , lambda6 = |  0  0  1 |  !!
!!           |  1  0  0 |             |  i  0  0 |             |  0  1  0 |  !!
!!                                                                           !!
!!           |  0  0  0 |                       |  1  0  0 |                 !!
!! lambda7 = |  0  0 -i | , lambda8 = 1/sqrt(3) |  0  1  0 | .               !!
!!           |  0  i  0 |                       |  0  0 -2 |                 !!
!!                                                                           !!
!!===========================================================================!!

      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      parameter(twopi = 6.28318530718)
      parameter(one_by_sqrt3 = 0.57735026919)
      parameter(two_by_sqrt3 = 1.154700538379)
      complex h(ncol,ncol)
      real temp(ncol2),phi,radius,raux
      real sigma

      do igen = 1,ncol2/2
         phi = twopi * ran2()
         raux = - sigma * log(ran2())
         radius = sqrt(raux)
         temp(2*igen-1) =  radius * cos(phi)
         temp(2*igen) = radius * sin(phi)
      enddo

      x1 = temp(3) + one_by_sqrt3 * temp(8)
      x2 = 0.0 
      h(1,1) = cmplx(x1,x2)  

      h(1,2) = cmplx(temp(1),-temp(2))

      h(1,3) = cmplx(temp(4),-temp(5))

      h(2,1) = cmplx(temp(1),temp(2))

      x1 = - temp(3) + one_by_sqrt3 * temp(8)
      x2 = 0.0 
      h(2,2) = cmplx(x1,x2)  

      h(2,3) = cmplx(temp(6),-temp(7))

      h(3,1) = cmplx(temp(4),temp(5))

      h(3,2) = cmplx(temp(6),temp(7))

      x1 = - two_by_sqrt3 * temp(8)
      x2 = 0.0 
      h(3,3) = cmplx(x1,x2)  

      return
      end

c =========================================================================
      subroutine gauss_vector(v,sigma)
c =========================================================================
!!===========================================================================!!
!! A gauss vector has complex components z which are gaussian distributed    !!
!! with <z~ z> = sigma                                                       !!
!!===========================================================================!!


      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      parameter(twopi = 6.28318530718)
c      parameter(sigma = 1.0000000) 
      complex v(ncol)
      real temp1,temp2,phi,radius,raux
      real sigma

      do icol = 1,ncol
         phi = twopi * ran2()
         raux = - sigma * log(ran2())
         radius = sqrt(raux)
         temp1 =  radius * cos(phi)
         temp2 = radius * sin(phi)
         v(icol) = cmplx(temp1,temp2)
      enddo
      
      return
      end
