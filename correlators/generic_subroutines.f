c============================================================================
c  Control if close
c  written by Sanfo
c============================================================================
      function imposed_close()

      integer imposed_close

      open(7,file='control_file',status='unknown')
      read(7,*,end=312) imposed_close
 312  close(7)

      return
      end



c============================================================================
      subroutine close_all(n_traj)
c  Call ending subroutines
c  written by Massimo D'Elia
c============================================================================
      implicit none
      include "parameters.f"

c     arguments
      integer n_traj


      write(*,*) "ENDING JOB ..."
      call addrem_stagphase

c     save lattice
c      call write_lattice(n_traj) !binary version
      call write_lattice2(n_traj) !text version, intersystem

c     save random number generator status
      call ranfinish 

c     create an empty file "allright"
      open(11,file='allright',status='unknown')
      close(11)

      write(*,*) "JOB ENDED"

      return
      end


c============================================================================
      subroutine write_lattice(n_traj)
c  Save link-array in a binary file
c  written by Massimo D'Elia
c============================================================================

      implicit none
      include "parameters.f"

c     arguments
      integer n_traj

c     common-block passed variables
      complex u
      common/field/u(ncol,ncol,4,nvol)

      open(2,file='lattice',status='unknown',form='unformatted')
      write(2) u
      write(2) n_traj
      close(2)

      return
      end

c============================================================================
      subroutine write_lattice2(n_traj)
c  Save link-array in a text file that can be 
c  easily opened by computer based on different architecure
c  written by Sanfo
c============================================================================

      implicit none
      include "parameters.f"

c     arguments
      integer n_traj

c     common-block passed variables
      complex u
      common/field/u(ncol,ncol,4,nvol)

c     internal variables
      integer icol1,icol2,idir,ivol
      complex c1

      open(8,file='lattice2',status='unknown')
      write(8,*) n_traj
      do icol1=1,ncol
         do icol2=1,ncol
            do idir=1,4
               do ivol=1,nvol
                  c1=u(icol1,icol2,idir,ivol)
                  write(8,*) real(c1),aimag(c1)
               enddo
            enddo
         enddo
      enddo
      close(8)

      return
      end

c=========================================================
      subroutine geometry
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
c  set lattice geometry
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      integer perm
      common/perm/ncfact,perm(nperm,ncol),sign(nperm)
      common/param/nran,beta_f,beta_a,rho,pi,iad_flag
      integer sind,coor,forw,back,sum,parity,par
      integer forweo,backeo,forweo2,backeo2,sindeo,sindeoh,cooreo
      common/geo2/forweo2(nvol,4),backeo2(nvol,4)
      common/ge/sind(nx,ny,nz,nt),coor(nvol,4),forw(nvol,4),back(nvol,4)
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),
     $  parity(nx,ny,nz,nt),cooreo(nvol,4),forweo(nvol,4),backeo(nvol,4)
      common/gegeo/eotolex(nvol),lextoeo(nvol)
      common/stagphase/eta(4,nvol)
      complex u,usave,staple
      common/field/u(ncol,ncol,4,nvol)
      common/field1/usave(ncol,ncol,4,nvol)
      common/staple/staple(ncol,ncol,4,nvol)
      real gen_fact
      common/cartan/gen_fact(ncol)

c     to begin let us define the factors multiplying
c     the diagonal cartan generators, used in gauss.f

      do icol = 1,ncol-1
         x = float(icol*(icol + 1)/2)
         gen_fact(icol) = 1.0/sqrt(x)
      enddo


      do it = 1,nt
         do iz = 1,nz
            do iy = 1,ny
               do ix = 1,nx
                  sum = ix + iy + iz + it
                  par = sum - 2*(sum/2)
                  i =  ix + (iy-1)*nx + (iz-1)*nx*ny + (it-1)*nx*ny*nz
                  j = (i-1)/2 + 1 + par*nvolh
                  sind(ix,iy,iz,it) = i 
                  sindeo(ix,iy,iz,it) = j
                  parity(ix,iy,iz,it) = par 
                  sindeoh(ix,iy,iz,it) = (i-1)/2 + 1
                  eotolex(sindeo(ix,iy,iz,it)) = sind(ix,iy,iz,it)
                  lextoeo(sind(ix,iy,iz,it)) = sindeo(ix,iy,iz,it)
                  coor(i,1) = ix
                  coor(i,2) = iy
                  coor(i,3) = iz
                  coor(i,4) = it
                  cooreo(j,1) = ix
                  cooreo(j,2) = iy
                  cooreo(j,3) = iz
                  cooreo(j,4) = it
               enddo
            enddo
         enddo
      enddo 

      do it = 1,nt
         do iz = 1,nz
            do iy = 1,ny
               do ix = 1,nx
                  ixp = ix + 1 
                  iyp = iy + 1 
                  izp = iz + 1 
                  itp = it + 1
                  if(ix.eq.nx) ixp = 1
                  if(iy.eq.ny) iyp = 1 
                  if(iz.eq.nz) izp = 1
                  if(it.eq.nt) itp = 1
                  ixm = ix - 1 
                  iym = iy - 1
                  izm = iz - 1
                  itm = it - 1
                  if(ixm.eq.0) ixm = nx
                  if(iym.eq.0) iym = ny
                  if(izm.eq.0) izm = nz
                  if(itm.eq.0) itm = nt
c funzioni spostamento per i siti cartesiani
                  forw(sind(ix,iy,iz,it),1) = sind(ixp,iy,iz,it)
                  forw(sind(ix,iy,iz,it),2) = sind(ix,iyp,iz,it)
                  forw(sind(ix,iy,iz,it),3) = sind(ix,iy,izp,it)
                  forw(sind(ix,iy,iz,it),4) = sind(ix,iy,iz,itp)
                  back(sind(ix,iy,iz,it),1) = sind(ixm,iy,iz,it)
                  back(sind(ix,iy,iz,it),2) = sind(ix,iym,iz,it)
                  back(sind(ix,iy,iz,it),3) = sind(ix,iy,izm,it)
                  back(sind(ix,iy,iz,it),4) = sind(ix,iy,iz,itm)
c funzioni spostamento per i siti even-odd
                  forweo(sindeo(ix,iy,iz,it),1) = sindeoh(ixp,iy,iz,it)
                  forweo(sindeo(ix,iy,iz,it),2) = sindeoh(ix,iyp,iz,it)
                  forweo(sindeo(ix,iy,iz,it),3) = sindeoh(ix,iy,izp,it)
                  forweo(sindeo(ix,iy,iz,it),4) = sindeoh(ix,iy,iz,itp)
                  backeo(sindeo(ix,iy,iz,it),1) = sindeoh(ixm,iy,iz,it)
                  backeo(sindeo(ix,iy,iz,it),2) = sindeoh(ix,iym,iz,it)
                  backeo(sindeo(ix,iy,iz,it),3) = sindeoh(ix,iy,izm,it)
                  backeo(sindeo(ix,iy,iz,it),4) = sindeoh(ix,iy,iz,itm)

                  forweo2(sindeo(ix,iy,iz,it),1) = sindeo(ixp,iy,iz,it)
                  forweo2(sindeo(ix,iy,iz,it),2) = sindeo(ix,iyp,iz,it)
                  forweo2(sindeo(ix,iy,iz,it),3) = sindeo(ix,iy,izp,it)
                  forweo2(sindeo(ix,iy,iz,it),4) = sindeo(ix,iy,iz,itp)
                  backeo2(sindeo(ix,iy,iz,it),1) = sindeo(ixm,iy,iz,it)
                  backeo2(sindeo(ix,iy,iz,it),2) = sindeo(ix,iym,iz,it)
                  backeo2(sindeo(ix,iy,iz,it),3) = sindeo(ix,iy,izm,it)
                  backeo2(sindeo(ix,iy,iz,it),4) = sindeo(ix,iy,iz,itm)

c                 definition of staggered phases                  
                  eta(1,sindeo(ix,iy,iz,it)) = 1.0
                  sum = ix - 1
                  eta(2,sindeo(ix,iy,iz,it)) = (-1.0)**sum
                  sum = ix - 1 + iy - 1
                  eta(3,sindeo(ix,iy,iz,it)) = (-1.0)**sum
                  sum = ix - 1 + iy - 1 + iz - 1 
                  eta(4,sindeo(ix,iy,iz,it)) = (-1.0)**sum
                  if(it.eq.nt) then  !! antiperiodic boundary conditions
                eta(4,sindeo(ix,iy,iz,it)) = -eta(4,sindeo(ix,iy,iz,it))
                  endif
               enddo
            enddo
         enddo
      enddo 

      return
      end

c=========================================================
      subroutine initialize_lattice(init_flag,n_traj)
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc 0 -> cold start; 1 -> hot start; 
cc 2 -> bin stored config; 3 -> text stored config
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      common/param/nran,beta_f,beta_a,rho,pi,iad_flag 
      integer sind,coor,forw,back,parity 
      integer forweo,backeo,forweo2,backeo2,sindeo,sindeoh,cooreo
      common/geo2/forweo2(nvol,4),backeo2(nvol,4)
      common/ge/sind(nx,ny,nz,nt),coor(nvol,4),forw(nvol,4),back(nvol,4)
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),
     $     parity(nx,ny,nz,nt),cooreo(nvol,4),
     $     forweo(nvol,4),backeo(nvol,4)
      common/gegeo/eotolex(nvol),lextoeo(nvol)
      common/stagphase/eta(4,nvol)
      complex u
      common/field/u(ncol,ncol,4,nvol)
      integer idir,ind,ix,iy,iz,it,n_traj

      if(init_flag.eq.0) then
         ind = 0
         do it = 1,nt
            do iz = 1,nz
               do iy = 1,ny
                  do ix = 1,nx
                     ind = ind + 1
                     do idir = 1,4
                        call one(u(1,1,idir,ind))
                     enddo
                  enddo
               enddo
            enddo
         enddo
         n_traj = 0
      endif
      
      if(init_flag.eq.1) then
         ind = 0
         do it = 1,nt
            do iz = 1,nz
               do iy = 1,ny
                  do ix = 1,nx
                     ind = ind + 1
                     do idir = 1,4
                        call sun_random(u(1,1,idir,ind))
                     enddo
                  enddo
               enddo
            enddo
         enddo
         n_traj = 0
      endif
      
      if(init_flag.eq.2) then
         open(2,file='lattice',status='old',form='unformatted')
         read(2) u
         read(2) n_traj
         close(2)
      endif
      
      if(init_flag.eq.3) then
         open(2,file='lattice2',status='old')
         read(2,*) n_traj
         do icol1=1,ncol
            do icol2=1,ncol
               do idir=1,4
                  do ivol=1,nvol
                     read(2,*) ur,ui
                     u(icol1,icol2,idir,ivol)=cmplx(ur,ui)
                  enddo
               enddo
            enddo
         enddo
         close(2)
      endif

      if(init_flag.gt.3) then
         write(*,*) "BAD INIT_FLAG! use 0,1,2,3 (cold,hot,stored)"
      endif

      return
      end

c=========================================================
      subroutine initialize_extf(val_ext_f)
c  written by Massimo D'Elia
c  version 1.0 - 2008
c=========================================================
      implicit none
      include "parameters.f"

c     arguments
      real val_ext_f

c     common-block passed variables
      complex ue1,ue2
      real u1phase
      common/fielde1/ue1(4,nvol)
      common/fielde2/ue2(4,nvol)
      common/expu1/u1phase(4,nvol)
      integer sindeo,sindeoh,parity,cooreo,forweo,backeo
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),
     $     parity(nx,ny,nz,nt),cooreo(nvol,4),
     $     forweo(nvol,4),backeo(nvol,4)

c     internal variables
      integer ind,ix,iy,iz,it,idir
      integer relaxed
      real bphase
      complex caux
      
      relaxed=1
      
      if(relaxed.eq.0) then
         bphase=2*val_ext_f*pigr/float(nx)
      else
         bphase=2*val_ext_f*pigr/float(nx**2)
      endif

      do it=1,nt
         do iz=1,nz
            do iy=1,ny
               do ix=1,nx
                  
                  ind=sindeo(ix,iy,iz,it)
!     defintions of the phases relative to the x and y directions
!     according to the appropriate quantization condition
                  if((relaxed.eq.1).and.(ix.eq.nx)) then
                     u1phase(1,ind)=-float(iy)*nx
                  else
                     u1phase(1,ind)=0
                  endif
                  u1phase(2,ind)=float(ix)
!     whatever the quantization condition, the phase in the z,t
!     directions are always 0
                  u1phase(3,ind)=0
                  u1phase(4,ind)=0
                     
                  do idir=1,4
                     
                     caux=cmplx(cos(u1phase(idir,ind)*bphase)
     $                    ,sin(u1phase(idir,ind)*bphase))
 
                     ue1(idir,ind)=caux**2
                     ue2(idir,ind)=conjg(caux)
                  enddo
               enddo
            enddo
         enddo
      enddo
      
      return
      end




c=========================================================
      subroutine add_rem_extf(act,ud)
c  written by Massimo D'Elia
c  version 1.0 - 2007/2008
c  add or remove external fields 1 or 2 according 
c  to passed arguments
c=========================================================
      implicit none
      include "parameters.f"
      
c     arguments
      integer act,ud

c     common-block passed variables
      complex ue1,ue2
      common/fielde1/ue1(4,nvol)
      common/fielde2/ue2(4,nvol)
      integer sindeo,sindeoh,parity,cooreo,forweo,backeo
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),
     $     parity(nx,ny,nz,nt),cooreo(nvol,4),
     $     forweo(nvol,4),backeo(nvol,4)
      complex u
      common/field/u(ncol,ncol,4,nvol)

c     internal variables
      integer eos
      integer ix,iy,iz,it,idir,icol1,icol2
      complex caux


      do it=1,nt
         do iz=1,nz
            do iy=1,ny
               do ix=1,nx

                  eos=sindeo(ix,iy,iz,it) 

                  do idir=1,4
                     do icol2=1,ncol
                        do icol1=1,ncol

                           if(ud.eq.1) then
                              caux=ue1(idir,eos)
                           else
                              caux=ue2(idir,eos)
                           endif

                           if(act.eq.2) then
                              caux=conjg(caux)
                           endif

                           u(icol1,icol2,idir,eos)=
     $                          caux*u(icol1,icol2,idir,eos)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
      end

c=========================================================
      subroutine add_extf(ud)
c  written by Sanfo
c=========================================================
      implicit none

c     arguments
      integer ud

      call add_rem_extf(1,ud)

      return
      end

c=========================================================
      subroutine rem_extf(ud)
c  written by Sanfo
c=========================================================
      implicit none

c     arguments
      integer ud

      call add_rem_extf(2,ud)

      return
      end

c=========================================================
      subroutine addrem_stagphase()
c  written by Massimo D'Elia
c  version 1.0 - 2007/2008
c  add or remove staggered phases
c=========================================================
      implicit none
      include "parameters.f"

c     common-block passed variables
      integer sindeo,sindeoh,parity,cooreo,forweo,backeo
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),
     $     parity(nx,ny,nz,nt),cooreo(nvol,4),
     $     forweo(nvol,4),backeo(nvol,4)
      real eta
      common/stagphase/eta(4,nvol)
      complex u
      common/field/u(ncol,ncol,4,nvol)

c     internal variables
      integer idir,ix,iy,iz,it,icol1,icol2,eos

      do it=1,nt
         do iz=1,nz
            do iy=1,ny
               do ix=1,nx
                  eos=sindeo(ix,iy,iz,it) 
                  do idir=1,4
                     do icol2=1,ncol
                        do icol1=1,ncol
                           u(icol1,icol2,idir,eos)=
     $                          eta(idir,eos)*u(icol1,icol2,idir,eos)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      
      return
      end


c=========================================================
      subroutine normalize_lattice()
c  written by Massimo D'Elia
c  version 1.0 - 2007/2008
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      common/param/nran,beta_f,beta_a,rho,pi,iad_flag 
      integer sind,coor,forw,back,parity
      integer forweo,backeo,forweo2,backeo2,sindeo,sindeoh,cooreo
      common/geo2/forweo2(nvol,4),backeo2(nvol,4)
      common/ge/sind(nx,ny,nz,nt),coor(nvol,4),forw(nvol,4),back(nvol,4)
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),
     $     parity(nx,ny,nz,nt),cooreo(nvol,4),
     $     forweo(nvol,4),backeo(nvol,4)
      common/gegeo/eotolex(nvol),lextoeo(nvol)
      common/stagphase/eta(4,nvol)
      complex u
      common/field/u(ncol,ncol,4,nvol)
      integer idir

         do ivol = 1,nvol
            do idir = 1,4
               call unitarize(u(1,1,idir,ivol))
            enddo
         enddo

      return
      end


c=========================================================
      subroutine normalize1_lattice()
c  written by Massimo D'Elia
c  version 1.0 - 2007/2008
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      common/param/nran,beta_f,beta_a,rho,pi,iad_flag 
      integer sind,coor,forw,back,parity
      integer forweo,backeo,forweo2,backeo2,sindeo,sindeoh,cooreo
      common/geo2/forweo2(nvol,4),backeo2(nvol,4)
      common/ge/sind(nx,ny,nz,nt),coor(nvol,4),forw(nvol,4),back(nvol,4)
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),
     $     parity(nx,ny,nz,nt),cooreo(nvol,4),
     $     forweo(nvol,4),backeo(nvol,4)
      common/gegeo/eotolex(nvol),lextoeo(nvol)
      common/stagphase/eta(4,nvol)
      complex u
      common/field/u(ncol,ncol,4,nvol)
      integer idir

         do ivol = 1,nvol
            do idir = 1,4
               call unitarize1(u(1,1,idir,ivol))
            enddo
         enddo

      return
      end



c=========================================================
      subroutine generate_permutations
c  written by Massimo D'Elia
c  version 1.0 - 11/07/2004
cc generate all permutations of the first ncol numbers and
cc their signs using an awful stupid algorithm 
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "parameters.f"
      integer perm
      common/perm/ncfact,perm(nperm,ncol),sign(nperm)
      
      integer perm1(ncol),perm11(ncol,ncol)
      integer lfact(ncol),lst(ncol),jref(ncol)
      real sign1(ncol)      
      
      do i1 = 1,ncol
         perm1(i1) = i1
         do i2 = 1,ncol
            perm11(i1,i2) = i2
         enddo
         lst(i1) = 1
         sign1(i1) = 1.
         jref(i1) = 1
      enddo
      
      k2 = ncfact
      do j = 1,ncol
         lfact(j) = k2
c         write(*,*) j,k2
         k2 = k2/(ncol - j + 1)
      enddo
      
      iperm = 1
      
 500  continue
      
      signaux = 1.0   
      do k = 1,ncol
         perm(iperm,k) = perm1(k)
         signaux = signaux*sign1(k)
c         write(*,*) iperm, signaux, perm1(k)
      enddo
      sign(iperm) = signaux
      
      
      do j = 2,ncol
         
         k2 = (iperm/lfact(j))*lfact(j) - iperm
         if(k2.eq.0) then
            do k = 1,ncol
               perm11(j,k) = perm11(jref(j),k)
            enddo
            sign1(j) = -1.0
            do k = j+1,ncol
               jref(k) = j
               sign1(k) = 1.0 
            enddo
            k3 = perm11(j,j-1+lst(j))
            perm11(j,j-1+lst(j)) = perm11(j,j-1)
            perm11(j,j-1) = k3
            lst(j) = lst(j) + 1
            do k1 = j+1,ncol
               lst(k1) = 1
            enddo
            do k = 1,ncol
               perm1(k) = perm11(j,k)
            enddo
            goto 510
         endif
         
         
      enddo   
      
      
      
 510  iperm = iperm + 1
c      write(*,*) "----------------"
      if (iperm.lt.ncfact) goto 500
      
      signaux = 1.0
      do k = 1,ncol
         perm(iperm,k) = perm1(k)
         signaux = signaux * sign1(k)
c         write(*,*) iperm, signaux, perm1(k)
      enddo
      sign(iperm) = signaux
      
      
      return
      end
c=========================================================

c============================================================================
c  RANDOM NUMBER GENERATOR: standard ran2 from numerical recipes
c============================================================================
      function ran2()
      implicit real (a-h,o-z)
      implicit integer (i-n)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real ran2,am,eps,rnmx
      parameter(im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     &          ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,
     &          ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,
     &          rnmx=1.-eps)
      integer idum2,j,k,iv,iy
      common /dasav/ idum,idum2,iv(ntab),iy
c      save iv,iy,idum2
c      data idum2/123456789/, iv/NTAB*0/, iy/0/

      if(idum.le.0) then
         idum=max0(-idum,1)
         idum2=idum
         do j=ntab+8,1,-1
            k=idum/iq1
            idum=ia1*(idum-k*iq1)-k*ir1
            if(idum.lt.0) idum=idum+im1
            if(j.le.ntab) iv(j)=idum
         enddo
         iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if(idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1) iy=iy+imm1
      ran2=min(am*iy,rnmx)

      return
      end

c=============================================================================
      subroutine ranstart
      implicit real (a-h,o-z)
      implicit integer (i-n)
      common /dasav/ idum,idum2,iv(32),iy

      open(unit=23, file='randomseed', status='unknown')
      read(23,*) idum
      read(23,*,end=117) idum2
      do i=1,32
         read(23,*) iv(i)
      enddo
      read(23,*) iy
      close(23)
      goto 118                          !!takes account of the first start
 117  if(idum.ge.0) idum = -idum -1     !!
      close(23)
 118  continue                          !!

      return
      end

c=============================================================================
      subroutine ranfinish
      implicit real (a-h,o-z)
      implicit integer (i-n)
      common /dasav/ idum,idum2,iv(32),iy

      open(unit=23, file='randomseed', status='unknown')
      write(23,*) idum
      write(23,*) idum2
      do i=1,32
         write(23,*) iv(i)
      enddo
      write(23,*) iy
      close(23)

      return
      end
c=============================================================================
