*=======================================================================
      function rani(lim)
*=======================================================================
      implicit none
      real ran2
      integer rani,lim
      
      rani=ran2()*lim

      return
      end
      
*=======================================================================
      function z2()
*=======================================================================
      implicit none
      real ran2
      real z2
      
      if(ran2().le.0.5) then
         z2=-1
      else
         z2=+1
      endif
      
      return
      end
      
*=======================================================================
      function z4()
*=======================================================================
      implicit none
      complex z4
      real z2
      
      z4=complex(z2(),z2())/sqrt(2.0)
      
      return
      end
      
*=======================================================================
      program meson_mass
*=======================================================================
*  Program to compute correlators of local staggered meson operators
*  along---
*                        t-direction
*                     ----------------  
*                 m              NO/O states
*                ----           -------------
*                 1             pi_2 / a0
*                 2             pi / --
*                 3             rho / a1
*                 4             rho / a1 
*                 5             rho / a1
*                 6             rho / b1
*                 7             rho / b1
*                 8             rho / b1
*
*  for flavor combination: (ubar u), (dbar d), (ubar d) & (dbar u)
*
*  Output files---
*
*   t-direction: meson_corr_t_ubu.dat,
*                meson_corr_t_dbd.dat,
*                meson_corr_t_ubd.dat.
*                meson_corr_t_dbu.dat.
*                (columns: 1 ==> t; 2-9 ==> results m1-m8)
*
*  Last modified : 05/11/2010
*-----------------------------------------------------------------------
*  Variable
*-----------------------------------------------------------------------
        implicit none
        include "parameters.f"
* common variable: geometry
        integer sind, coor, forw, back
        common/ge/ sind(nx,ny,nz,nt), coor(nvol,4),
     $             forw(nvol,4), back(nvol,4)     
      integer sindeo,sindeoh,parity,cooreo,forweo,backeo
      common/geo/sindeo(nx,ny,nz,nt),sindeoh(nx,ny,nz,nt),
     $     parity(nx,ny,nz,nt),cooreo(nvol,4),
     $     forweo(nvol,4),backeo(nvol,4)
* internal variable
        integer init_flag, termalizza, n_meas, n_rand, n_traj
        real immu_quark, immu_iso, val_ext_f
        character(LEN=30) file_rhmc

        complex prop_u(ncol,ncol,nvol), prop_d(ncol,ncol,nvol)
        real tr_ubu(nvol), tr_dbd(nvol), tr_ubd(nvol), tr_dbu(nvol)
        integer meson_phase(8,nvol)
        real meson_corr_t(4,8,nt)
        integer ix, iy, iz, it, ivol,icol
        complex u
        common/field/u(ncol,ncol,4,nvol)
        complex ori_source_e(nvolh),ori_source_o(nvolh)
        real normq
        integer site,only_vertexes,vertex,i_rand
        integer ivol_dest,tshift,rani
        complex z4
        parameter(only_vertexes=0)
        complex utemp(ncol,ncol,4,nvol)
        real plaq_s,plaq_t,plaq2_s,plaq2_t
*-----------------------------------------------------------------------
*  Initialize
*-----------------------------------------------------------------------
        call load_par( init_flag, termalizza, n_meas, n_rand, 
     $                 immu_quark, immu_iso, file_rhmc, val_ext_f )

        call set_const( immu_quark, immu_iso )

        call geometry

        call initialize_extf( val_ext_f )
        
        call initialize_lattice( init_flag, n_traj )
        call ranstart

        meson_corr_t = 0.0

        write(*,*)
        write(*,*) "---------Executing: meson_mass---------------"
        write(*,*) " Trajectory Number: ", n_traj
        write(*,*)

*-----------------------------------------------------------------------

        call plaquette(plaq_s,plaq_t,plaq2_s,plaq2_t) 

        do i_rand=1,n_rand

c     shift the configuration arbitrarily in the t direction
           call addrem_stagphase
           tshift=rani(nt)+1
           do ix=1,nx
              do iy=1,ny
                 do iz=1,nz
                    do it=1,nt
                       ivol=sindeo(ix,iy,iz,it)
                       ivol_dest=sindeo(ix,iy,iz,mod(it+tshift,nt)+1)
                       utemp(:,:,:,ivol_dest)=u(:,:,:,ivol)
                    enddo
                 enddo
              enddo
           enddo
           u=utemp
           call addrem_stagphase
           
c     create the source
           normq=0
           ori_source_e=0
           ori_source_o=0
           it=1        
           do ix=1,nx
              do iy=1,ny
                 do iz=1,nz
                
                    site=sindeo(ix,iy,iz,it)
                    vertex=mod(ix-1,2)+mod(iy-1,2)+mod(iz-1,2)
                    if(only_vertexes.ne.1.or.vertex.eq.0) then

                       if(site.le.nvolh) then
                          ori_source_e(site)=z4()
                       else
                          ori_source_o(site-nvolh)=z4()
                       endif

                       normq=normq+1

                    endif
                    
                 enddo
              enddo
           enddo

           ori_source_e=ori_source_e/sqrt(normq)
           ori_source_o=ori_source_o/sqrt(normq)

           write(*,*) "Source",i_rand,"/",n_rand
           
c           ori_source_e=0
c           ori_source_o=0
c           ori_source_e(1)=1
*-----------------------------------------------------------------------
*  Compute lexically ordered quark propagators : 
*                       M^(-1) ( x,y,z,t ; 1,1,1,1 )_[i,color]
*  for all color
*-----------------------------------------------------------------------
* u-quark
           
           call add_extf( 1 )
           
           do icol = 1, ncol
              call quark_propagator( icol, prop_u(icol,:,:),ori_source_e
     $             ,ori_source_o )
           enddo
           
           call rem_extf( 1 )
           
* d-quark
           
           call add_extf( 2 )
           
           do icol = 1, ncol
              call quark_propagator( icol, prop_d(icol,:,:),ori_source_e
     $             ,ori_source_o)
           enddo
           
           call rem_extf( 2 )
*-----------------------------------------------------------------------
*  Compute :
*     Tr_c { [ M^(-1)(x,y,z,t;1,1,1,1) ]+ * M^(-1)(x,y,z,t;1,1,1,1) }
*  for each (x,y,z,t)
*-----------------------------------------------------------------------
* ubar u
           
           tr_ubu = 0.0
           
           do icol = 1, ncol
              tr_ubu(:) = tr_ubu(:) + 
     $             abs( prop_u(icol,1,:) ) ** 2 + 
     $             abs( prop_u(icol,2,:) ) ** 2 + 
     $             abs( prop_u(icol,3,:) ) ** 2
           enddo
           
* dbar d
           
           tr_dbd = 0.0

           do icol = 1, ncol
              tr_dbd(:) = tr_dbd(:) + 
     $             abs( prop_d(icol,1,:) ) ** 2 + 
     $             abs( prop_d(icol,2,:) ) ** 2 + 
     $             abs( prop_d(icol,3,:) ) ** 2
           enddo
           
* ubar d
           
           tr_ubd = 0.0
           
           do icol = 1, ncol
              tr_ubd(:) = tr_ubd(:) + 
     $             real( conjg(prop_u(icol,1,:)) * prop_d(icol,1,:) ) + 
     $             real( conjg(prop_u(icol,2,:)) * prop_d(icol,2,:) ) + 
     $             real( conjg(prop_u(icol,3,:)) * prop_d(icol,3,:) )
           enddo

* dbar u
           
           tr_dbu = 0.0
           
           do icol = 1, ncol
              tr_dbu(:) = tr_dbu(:) + 
     $             real( conjg(prop_d(icol,1,:)) * prop_u(icol,1,:) ) + 
     $             real( conjg(prop_d(icol,2,:)) * prop_u(icol,2,:) ) + 
     $             real( conjg(prop_d(icol,3,:)) * prop_u(icol,3,:) )
           enddo
*-----------------------------------------------------------------------
*  Build meson correnxrs along t-direction       
*-----------------------------------------------------------------------
           call meson_phase_t( meson_phase )
           
           do it = 1, nt        !! t - loop
              
              do iz = 1, nz     !! loop over othogonal 3-volume
                 do iy = 1, ny
                    do ix = 1, nx
                       
                       ivol = sind(ix,iy,iz,it)
                       meson_corr_t(1,:,it) = meson_corr_t(1,:,it) + 
     $                      meson_phase(:,ivol) * tr_ubu(ivol) 
                       meson_corr_t(2,:,it) = meson_corr_t(2,:,it) + 
     $                      meson_phase(:,ivol) * tr_dbd(ivol) 
                       meson_corr_t(3,:,it) = meson_corr_t(3,:,it) + 
     $                      meson_phase(:,ivol) * tr_ubd(ivol) 
                       meson_corr_t(4,:,it) = meson_corr_t(4,:,it) + 
     $                      meson_phase(:,ivol) * tr_dbu(ivol) 
                       
                    enddo
                 enddo  
              enddo             !! loop over othogonal 3-volume
              
           enddo                !! t - loop
           
        enddo                   !! loop over sources
           
*-----------------------------------------------------------------------
*  Write-out results
*-----------------------------------------------------------------------
        open( unit=101, file='meson_corr_t_ubu.dat', status='unknown' ) 
        open( unit=102, file='meson_corr_t_dbd.dat', status='unknown' ) 
        open( unit=103, file='meson_corr_t_ubd.dat', status='unknown' ) 
        open( unit=104, file='meson_corr_t_dbu.dat', status='unknown' ) 

        write(101,*) "# n_traj = ", n_traj
        write(102,*) "# n_traj = ", n_traj
        write(103,*) "# n_traj = ", n_traj
        write(104,*) "# n_traj = ", n_traj
        do it = 1, nt
          write(101,10) it, meson_corr_t(1,:,it)/n_rand
          write(102,10) it, meson_corr_t(2,:,it)/n_rand
          write(103,10) it, meson_corr_t(3,:,it)/n_rand
          write(104,10) it, meson_corr_t(4,:,it)/n_rand
        enddo  

10      format( i2, 8(5x,e12.5) )          
        close(101); close(102); close(103); close(104)
*-----------------------------------------------------------------------
*  End program
*-----------------------------------------------------------------------
        stop
        end        
*-----------------------------------------------------------------------
*  Include all required subroutines
*-----------------------------------------------------------------------
        include "measure_subroutines.f"
        include "vecn_subroutines.f"
        include "sun_subroutines.f"
        include "sun_ext_subroutines.f"
        include "init_modified.f"
        include "generic_subroutines.f"
        include "dirac_matrix.f"
        include "inverter.f"
*=======================================================================


*=======================================================================
        subroutine quark_propagator( color, prop , ori_source_e,
     $       ori_source_o )
*=======================================================================
*  Computes lexically ordered quark propagators :  
*                       M^(-1) ( x,y,z,t ; 1,1,1,1 )_[i,color]
*  for color indices [i,color], for all i.
*
*  Last modified : 01/22/2009
*-----------------------------------------------------------------------
*  Variable
*-----------------------------------------------------------------------
        implicit none 
        include "parameters.f"
* argument
        integer color                   !! input
        complex prop(ncol,nvol)         !! output
        complex ori_source_e(nvolh)
        complex ori_source_o(nvolh)
* common variable: parameters
        real mass, mass2, residue
        common/param2/ mass, mass2, residue
        real immu1, immu2
        complex eim1, emim1, eim2, emim2
        common/param4/ immu1, immu2, eim1, emim1, eim2, emim2
* common variable: geometry
        integer sind(nx,ny,nz,nt), coor(nvol,4),

     $          forw(nvol,4), back(nvol,4)     
        common/ge/ sind, coor, forw, back
        integer eotolex(nvol), lextoeo(nvol)
        common/gegeo/ eotolex, lextoeo
* common variable : source
        complex source_e(ncol,nvolh), source_o(ncol,nvolh)
        common/source/ source_e, source_o
* internal variable
        integer ivolh, ivolh_pl_nvolh, ix, iy, iz, it, ivol, ieo
        real inv_mass
        complex inv_e(ncol,nvolh), inv_o(ncol,nvolh), chi(ncol,nvolh),
     $          inv_eo(ncol,nvol)
        
*-----------------------------------------------------------------------
*  Set point source S at (x=1,y=1,z=1,t=1) for color_index = color
*-----------------------------------------------------------------------
        call make_source( color, ori_source_e,ori_source_o )
*-----------------------------------------------------------------------
*  Compute :          I = M^(-1) * S 
*  for all even (e) and odd (o) sites. 
*  I = ( Ie  Io ) , inverse vector
*  S = ( Se  So ) , source vector
*  M = (  m    Deo/2 ) 
*      ( Doe/2    m  ) , with Deo = -Doe+   
*-----------------------------------------------------------------------
*  chi = M+ * S, on even sites
        
        call D( EO, chi, source_o, eim1, emim1 )

        chi = mass * source_e - 0.5 * chi

*  I = (M+ * M)^(-1) * chi = M^(-1) * S, on even sites

        call singol_inverter( inv_e, chi, eim1, emim1, mass2 )

*  Io using Ie

        call D( OE, chi, inv_e, eim1, emim1 )

        inv_mass = 1.0 / mass
        inv_o = inv_mass * ( source_o - 0.5 * chi )
*-----------------------------------------------------------------------
*  Reorder I=(Ie  Io) into lexical ordering   
*-----------------------------------------------------------------------
*  I from Ie (comes first) and Io

        do ivolh = 1, nvolh
          ivolh_pl_nvolh = ivolh + nvolh
          inv_eo(:,ivolh) = inv_e(:,ivolh) 
          inv_eo(:,ivolh_pl_nvolh) = inv_o(:,ivolh) 
        enddo

*  I to lexical ordering 

        do it = 1, nt
          do iz = 1, nz
            do iy = 1, ny
              do ix = 1, nx
                ivol = sind(ix,iy,iz,it)
                ieo = lextoeo(ivol)
                prop(:,ivol) = inv_eo(:,ieo)
              enddo
            enddo
          enddo
        enddo
*-----------------------------------------------------------------------
*  End subroutine
*-----------------------------------------------------------------------
        return
        end
*=======================================================================


*=======================================================================
        subroutine make_source( color , ori_source_e,ori_source_o )
*=======================================================================
*  Sets point sources at (x=1,y=1,z=1,t=1) for color_index = color
*
*  Last modified : 01/22/2009
*-----------------------------------------------------------------------
*  Variable
*-----------------------------------------------------------------------
        implicit none
        include "parameters.f"
* argument
        integer color           !! input
        complex ori_source_e(nvolh)
        complex ori_source_o(nvolh)
* common variable : source      !! output
        complex source_e(ncol,nvolh), source_o(ncol,nvolh)
        common/source/ source_e, source_o
* internal variable
        integer ivolh
*-----------------------------------------------------------------------
        do ivolh=1,nvolh
           source_e(color,ivolh) = ori_source_e(ivolh)
           source_o(color,ivolh) = ori_source_o(ivolh)
        enddo
*-----------------------------------------------------------------------
*  End subroutine
*-----------------------------------------------------------------------
        return
        end        
*=======================================================================


*=======================================================================
        subroutine meson_phase_t( meson_phase )
*=======================================================================
*  Computes phases for 8 local staggered meson operators for
*  correnxrs along :      t-direction
*
*  Last modified : 06/11/2010
*-----------------------------------------------------------------------
*  Variable
*-----------------------------------------------------------------------
        implicit none
        include "parameters.f"
* argument:
        integer meson_phase(8,nvol)     !! output
* common variable: geometry
        integer sind(nx,ny,nz,nt), coor(nvol,4),
     $          forw(nvol,4), back(nvol,4)     
        common/ge/ sind, coor, forw, back
* internal variable
        integer ix, iy, iz, it, ivol, jx ,jy, jz
*-----------------------------------------------------------------------
*  Set meson phases :
*
*  m              phase = eps(n)*phi(n)*(-)^t           NO/O states
* ----           ------------------------------        -------------
*  1              (-)^x+y+z                            pi_2 / a0
*  2              +1                                   pi / --
*  3              (-)^y+z                              rho / a1
*  4              (-)^x+z                              rho / a1 
*  5              (-)^x+y                              rho / a1
*  6              (-)^x                                rho / b1
*  7              (-)^y                                rho / b1
*  8              (-)^z                                rho / b1
*
*  where x, y, z, t should be treated as the relative distance form the
*  origin. Since in our case the origin is at (x_o=1,y_o=1,z_o=1,t_o=1)
*  in practice our x=x-1, y=y-1, z=z-1, t=t-1. 
*-----------------------------------------------------------------------
        meson_phase = 1

        do it = 1, nt

          do iz = 1, nz
            jz = - (-1) ** iz

            do iy = 1, ny
              jy = - (-1) ** iy

              do ix = 1, nx
                jx = - (-1) ** ix
                ivol = sind(ix,iy,iz,it)

                meson_phase(1,ivol) = jx * jy * jz
                meson_phase(2,ivol) = 1 
                meson_phase(3,ivol) = jy * jz 
                meson_phase(4,ivol) = jx * jz 
                meson_phase(5,ivol) = jx * jy
                meson_phase(6,ivol) = jx
                meson_phase(7,ivol) = jy
                meson_phase(8,ivol) = jz
              enddo

            enddo

          enddo

        enddo  
*-----------------------------------------------------------------------
*  End subroutine
*-----------------------------------------------------------------------
        return
        end        
*=======================================================================
