*=======================================================================
        program meson_corr
*=======================================================================
*  Program to compute correlators of local staggered meson operators
*  along---
*                        z-direction
*                     ----------------  
*                 m              NO/O states
*                ----           -------------
*                 1             pi_2 / a0
*                 2             pi / --
*                 3             rho_T / a1_T
*                 4             rho_T / a1_T 
*                 5             rho_L / a1_L
*                 6             rho_T / b1_T
*                 7             rho_T / b1_T
*                 8             rho_L / b1_L
*
*                        x-direction
*                     -----------------  
*                 m              NO/O states
*                ----           -------------
*                 1             pi_2 / a0
*                 2             pi / --
*                 3             rho_L / a1_L
*                 4             rho_T / a1_T 
*                 5             rho_T / a1_T
*                 6             rho_L / b1_L
*                 7             rho_T / b1_T
*                 8             rho_T / b1_T
*
*  for flavor combination: (ubar u), (dbar d) & (ubar d).
*
*  Output files---
*
*   z-direction: meson_corr_z_ubu.dat,
*                meson_corr_z_dbd.dat,
*                meson_corr_z_ubd.dat.
*                (columns: 1 ==> z; 2-9 ==> results m1-m8)
*
*   x-direction: meson_corr_x_ubu.dat,
*                meson_corr_x_dbd.dat,
*                meson_corr_x_ubd.dat.
*                (columns: 1 ==> x; 2-9 ==> results m1-m8)
*
*  Last modified : 01/22/2009
*-----------------------------------------------------------------------
*  Variable
*-----------------------------------------------------------------------
        implicit none
        include "parameters.f"
* common variable: geometry
        integer sind, coor, forw, back
        common/ge/ sind(nx,ny,nz,nt), coor(nvol,4),
     $             forw(nvol,4), back(nvol,4)     
* internal variable
        integer init_flag, termalizza, n_meas, n_rand, n_traj
        real immu_quark, immu_iso, val_ext_f
        character(LEN=30) file_rhmc

        complex prop_u(ncol,ncol,nvol), prop_d(ncol,ncol,nvol)
        real tr_ubu(nvol), tr_dbd(nvol), tr_ubd(nvol)
        integer meson_phase(8,nvol)
        real meson_corr_z(3,8,nz), meson_corr_x(3,8,nx)

        integer ix, iy, iz, it, ivol, icol
*-----------------------------------------------------------------------
*  Initialize
*-----------------------------------------------------------------------
        call load_par( init_flag, termalizza, n_meas, n_rand, 
     $                 immu_quark, immu_iso, file_rhmc, val_ext_f )

        call set_const( immu_quark, immu_iso )

        call initialize_extf( val_ext_f )

        call geometry

        call initialize_lattice( init_flag, n_traj )

        call addrem_stagphase
*-----------------------------------------------------------------------
        write(*,*)
        write(*,*) "---------Executing: meson_corr---------------"
        write(*,*) " Trajectory Number: ", n_traj
        write(*,*)
*-----------------------------------------------------------------------
*  Compute lexically ordered quark propagators : 
*                       M^(-1) ( x,y,z,t ; 1,1,1,1 )_[i,color]
*  for all color
*-----------------------------------------------------------------------
* u-quark

        call add_extf( 1 )

        do icol = 1, ncol
          call quark_propagator( icol, prop_u(icol,:,:) )
        enddo

        call rem_extf( 1 )

* d-quark

        call add_extf( 2 )

        do icol = 1, ncol
          call quark_propagator( icol, prop_d(icol,:,:) )
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
     $                abs( prop_u(icol,1,:) ) ** 2 + 
     $                abs( prop_u(icol,2,:) ) ** 2 + 
     $                abs( prop_u(icol,3,:) ) ** 2
        enddo

* dbar d

        tr_dbd = 0.0

        do icol = 1, ncol
          tr_dbd(:) = tr_dbd(:) + 
     $                abs( prop_d(icol,1,:) ) ** 2 + 
     $                abs( prop_d(icol,2,:) ) ** 2 + 
     $                abs( prop_d(icol,3,:) ) ** 2
        enddo

* ubar d

        tr_ubd = 0.0

        do icol = 1, ncol
          tr_ubd(:) = tr_ubd(:) + 
     $             real( conjg(prop_u(icol,1,:)) * prop_d(icol,1,:) ) + 
     $             real( conjg(prop_u(icol,2,:)) * prop_d(icol,2,:) ) + 
     $             real( conjg(prop_u(icol,3,:)) * prop_d(icol,3,:) )
        enddo
*-----------------------------------------------------------------------
*  Build meson correlators along z-direction       
*-----------------------------------------------------------------------
        call meson_phase_z( meson_phase )

        meson_corr_z = 0.0

        do iz = 1, nz                   !! z - loop

          do it = 1, nt                 !! loop over othogonal 3-volume
            do iy = 1, ny
              do ix = 1, nx

                ivol = sind(ix,iy,iz,it)
                meson_corr_z(1,:,iz) = meson_corr_z(1,:,iz) + 
     $                              meson_phase(:,ivol) * tr_ubu(ivol) 
                meson_corr_z(2,:,iz) = meson_corr_z(2,:,iz) + 
     $                              meson_phase(:,ivol) * tr_dbd(ivol) 
                meson_corr_z(3,:,iz) = meson_corr_z(3,:,iz) + 
     $                              meson_phase(:,ivol) * tr_ubd(ivol) 

              enddo
            enddo  
          enddo                         !! loop over othogonal 3-volume

        enddo                           !! z - loop
*-----------------------------------------------------------------------
*  Build meson correlators along x-direction       
*-----------------------------------------------------------------------
        call meson_phase_x( meson_phase )

        meson_corr_x = 0.0

        do ix = 1, nx                   !! x - loop

          do it = 1, nt                 !! loop over othogonal 3-volume
            do iz = 1, nz
              do iy = 1, ny

                ivol = sind(ix,iy,iz,it)
                meson_corr_x(1,:,ix) = meson_corr_x(1,:,ix) + 
     $                              meson_phase(:,ivol) * tr_ubu(ivol) 
                meson_corr_x(2,:,ix) = meson_corr_x(2,:,ix) + 
     $                              meson_phase(:,ivol) * tr_dbd(ivol) 
                meson_corr_x(3,:,ix) = meson_corr_x(3,:,ix) + 
     $                              meson_phase(:,ivol) * tr_ubd(ivol) 

              enddo
            enddo  
          enddo                         !! loop over othogonal 3-volume

        enddo                           !! x - loop
*-----------------------------------------------------------------------
*  Write-out results
*-----------------------------------------------------------------------
        open( unit=101, file='meson_corr_z_ubu.dat', status='unknown' ) 
        open( unit=102, file='meson_corr_z_dbd.dat', status='unknown' ) 
        open( unit=103, file='meson_corr_z_ubd.dat', status='unknown' ) 

        open( unit=201, file='meson_corr_x_ubu.dat', status='unknown' ) 
        open( unit=202, file='meson_corr_x_dbd.dat', status='unknown' ) 
        open( unit=203, file='meson_corr_x_ubd.dat', status='unknown' ) 

        write(101,*) "# n_traj = ", n_traj
        write(102,*) "# n_traj = ", n_traj
        write(103,*) "# n_traj = ", n_traj
        do iz = 1, nz
          write(101,10) iz, meson_corr_z(1,:,iz)
          write(102,10) iz, meson_corr_z(2,:,iz)
          write(103,10) iz, meson_corr_z(3,:,iz)
        enddo  

        write(201,*) "# n_traj = ", n_traj
        write(202,*) "# n_traj = ", n_traj
        write(203,*) "# n_traj = ", n_traj
        do ix = 1, nx
          write(201,10) ix, meson_corr_x(1,:,ix)
          write(202,10) ix, meson_corr_x(2,:,ix)
          write(203,10) ix, meson_corr_x(3,:,ix)
        enddo  

10      format( i2, 8(5x,e12.5) )          
        close(101); close(102); close(103)
        close(201); close(202); close(203)
*-----------------------------------------------------------------------
*  End program
*-----------------------------------------------------------------------
        stop
        end        
*-----------------------------------------------------------------------
*  Include all required subroutines
*-----------------------------------------------------------------------
        include "vecn_subroutines.f"
        include "sun_subroutines.f"
        include "sun_ext_subroutines.f"
        include "init_modified.f"
        include "generic_subroutines.f"
        include "dirac_matrix.f"
        include "inverter.f"
*=======================================================================


*=======================================================================
        subroutine quark_propagator( color, prop )
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
        call pt_source( color )
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
        subroutine pt_source( color )
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
* common variable : source      !! output
        complex source_e(ncol,nvolh), source_o(ncol,nvolh)
        common/source/ source_e, source_o
*-----------------------------------------------------------------------
*  Set point source at : (x=1,y=1,z=1,t=1) for color_index = color
*-----------------------------------------------------------------------
        source_e = cmplx(0.0,0.0) 
        source_o = cmplx(0.0,0.0) 

        source_e(color,1) = cmplx(1.0,0.0) 
*-----------------------------------------------------------------------
*  End subroutine
*-----------------------------------------------------------------------
        return
        end        
*=======================================================================


*=======================================================================
        subroutine meson_phase_z( meson_phase )
*=======================================================================
*  Computes phases for 8 local staggered meson operators for
*  correlators along :      z-direction
*
*  Last modified : 01/22/2009
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
        integer ix, iy, iz, it, ivol, jx ,jy, jt
*-----------------------------------------------------------------------
*  Set meson phases :
*
*  m              phase = eps(n)*phi(n)*(-)^z           NO/O states
* ----           ------------------------------        -------------
*  1              (-)^x+y+t                            pi_2 / a0
*  2              +1                                   pi / --
*  3              (-)^y+t                              rho_T / a1_T
*  4              (-)^x+t                              rho_T / a1_T 
*  5              (-)^x+y                              rho_L / a1_L
*  6              (-)^x                                rho_T / b1_T
*  7              (-)^y                                rho_T / b1_T
*  8              (-)^t                                rho_L / b1_L
*
*  where x, y, z, t should be treated as the relative distance form the
*  origin. Since in our case the origin is at (x_o=1,y_o=1,z_o=1,t_o=1)
*  in practice our x=x-1, y=y-1, z=z-1, t=t-1. 
*-----------------------------------------------------------------------
        meson_phase = 1

        do it = 1, nt
          jt = - (-1) ** it

          do iz = 1, nz

            do iy = 1, ny
              jy = - (-1) ** iy

              do ix = 1, nx
                jx = - (-1) ** ix
                ivol = sind(ix,iy,iz,it)

                meson_phase(1,ivol) = jx * jy * jt 
                meson_phase(2,ivol) = 1 
                meson_phase(3,ivol) = jy * jt 
                meson_phase(4,ivol) = jx * jt 
                meson_phase(5,ivol) = jx * jy
                meson_phase(6,ivol) = jx
                meson_phase(7,ivol) = jy
                meson_phase(8,ivol) = jt
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
        subroutine meson_phase_x( meson_phase )
*=======================================================================
*  Computes phases for 8 local staggered meson operators for
*  correlators along :      x-direction
*
*  Last modified : 01/22/2009
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
        integer ix, iy, iz, it, ivol, jy, jz, jt
*-----------------------------------------------------------------------
*  Set meson phases :
*
*  m              phase = eps(n)*phi(n)*(-)^x           NO/O states
* ----           ------------------------------        -------------
*  1              (-)^y+z+t                            pi_2 / a0
*  2              +1                                   pi / --
*  3              (-)^y+z                              rho_L / a1_L
*  4              (-)^z+t                              rho_T / a1_T 
*  5              (-)^y+t                              rho_T / a1_T
*  6              (-)^t                                rho_L / b1_L
*  7              (-)^y                                rho_T / b1_T
*  8              (-)^z                                rho_T / b1_T
*
*  where x, y, z, t should be treated as the relative distance form the
*  origin. Since in our case the origin is at (x_o=1,y_o=1,z_o=1,t_o=1)
*  in practice our x=x-1, y=y-1, z=z-1, t=t-1. 
*-----------------------------------------------------------------------
        meson_phase = 1

        do it = 1, nt
          jt = - (-1) ** it

          do iz = 1, nz
            jz = - (-1) ** iz

            do iy = 1, ny
              jy = - (-1) ** iy

              do ix = 1, nx
                ivol = sind(ix,iy,iz,it)

                meson_phase(1,ivol) = jy * jz * jt 
                meson_phase(2,ivol) = 1 
                meson_phase(3,ivol) = jy * jz 
                meson_phase(4,ivol) = jz * jt 
                meson_phase(5,ivol) = jy * jt
                meson_phase(6,ivol) = jt
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
