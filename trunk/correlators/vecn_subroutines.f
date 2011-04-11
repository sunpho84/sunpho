c=========================================================
      subroutine vmult(v,u,w)
c  written by Massimo D'Elia
c  version 1.0 - 
cc multiplication of an su(n) matrices by a su(n) vector: v = u * w
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex u(ncol,ncol),v(ncol),w(ncol)

      do i = 1,ncol
         v(i) = (0.,0.)
         do l = 1,ncol
            v(i) = v(i) + u(i,l)*w(l)
         enddo
      enddo
   
      return
      end  
c=========================================================
c=========================================================
      subroutine vhmult(v,u,w)
c  written by Massimo D'Elia
cc multiplication of an su(n) matrices by a su(n) vector: v = u~ * w
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex u(ncol,ncol),v(ncol),w(ncol)

      do i = 1,ncol
         v(i) = (0.,0.)
         do l = 1,ncol
            v(i) = v(i) + CONJG(u(l,i))*w(l)
         enddo
      enddo
   
      return
      end  
c=========================================================
c=========================================================
      subroutine vmult_add(v,u,w)
c  written by Massimo D'Elia
c  version 1.0 - 
cc add the multiplication of an su(n) matrices by a su(n) vector: v = u * w
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex u(ncol,ncol),v(ncol),w(ncol)

      do i = 1,ncol
         do l = 1,ncol
            v(i) = v(i) + u(i,l)*w(l)
         enddo
      enddo
   
      return
      end  
c=========================================================
c=========================================================
      subroutine vhmult_add(v,u,w)
c  written by Massimo D'Elia
cc add the multiplication of an su(n) matrices by a su(n) vector: v = u~ * w
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex u(ncol,ncol),v(ncol),w(ncol)

      do i = 1,ncol
         do l = 1,ncol
            v(i) = v(i) + CONJG(u(l,i))*w(l)
         enddo
      enddo
   
      return
      end  
c=========================================================
c=========================================================
      subroutine vmult_sub(v,u,w)
c  written by Massimo D'Elia
c  version 1.0 - 
cc add the multiplication of an su(n) matrices by a su(n) vector: v = u * w
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex u(ncol,ncol),v(ncol),w(ncol)

      do i = 1,ncol
         do l = 1,ncol
            v(i) = v(i) - u(i,l)*w(l)
         enddo
      enddo
   
      return
      end  
c=========================================================
c=========================================================
      subroutine vhmult_sub(v,u,w)
c  written by Massimo D'Elia
cc add the multiplication of an su(n) matrices by a su(n) vector: v = u~ * w
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex u(ncol,ncol),v(ncol),w(ncol)

      do i = 1,ncol
         do l = 1,ncol
            v(i) = v(i) - CONJG(u(l,i))*w(l)
         enddo
      enddo
   
      return
      end  
c=========================================================
c=========================================================
      subroutine tensn(u,v,w)
c  written by Massimo D'Elia
cc tensor product of two sun vectors 
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex u(ncol,ncol),v(ncol),w(ncol)

      do i = 1,ncol
         do j = 1,ncol
            u(i,j) = v(i) * w(j)
         enddo
      enddo
   
      return
      end  
c=========================================================
c=========================================================
      subroutine tensn_add(u,v,w)
c  written by Massimo D'Elia
cc tensor product of two sun vectors 
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      complex u(ncol,ncol),v(ncol),w(ncol)

      do i = 1,ncol
         do j = 1,ncol
            u(i,j) = u(i,j) + v(i) * w(j)
         enddo
      enddo
   
      return
      end  
c=========================================================
