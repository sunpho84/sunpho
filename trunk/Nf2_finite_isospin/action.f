c=========================================================
      subroutine uaction(action,beta)
c  written by Massimo D'Elia
c  version 1.0 - 01/08
c  local plaquette for molecular dynamics updating
c=========================================================
      implicit none
      include "parameters.f"

c     arguments
      real action(nvol),beta

c     internal variables
      integer ivol
      real locplaq(nvol)


      call local_plaq(locplaq)

      do ivol=1,nvol
         action(ivol)=action(ivol)+locplaq(ivol)*beta
      enddo


      return
      end


c=========================================================
      subroutine qaction_flav(action,phi_e,phi_o,emu,ememu)
c  written by Massimo D'Elia
c  version 1.0 - 01/08
c  local action for ud flavour
c=========================================================

      implicit none
      include "parameters.f"

c     arguments
      real action(nvol)
      complex phi_e(ncol,nvolh),phi_o(ncol,nvolh)
      real emu,ememu
      
c     common blocks
      integer nterm
      real cost,pole,coef
      common/rhmc_a/nterm,cost,pole(nmr),coef(nmr)
      
c     internal variables
      integer ivol,icol
      complex chi_e(ncol,nvolh),chi_o(ncol,nvolh)


      call multi_shift_summed_inverter(chi_e,chi_o,phi_e,phi_o,emu,ememu
     $     ,cost,pole,coef,nterm)

      do ivol=1,nvol

         do icol=1,ncol
            if(ivol.le.nvolh) then
               action(ivol)=action(ivol)+real(chi_e(icol,ivol)
     $              *conjg(phi_e(icol,ivol)))
            else
               action(ivol)=action(ivol)+real(chi_o(icol,ivol-nvolh)
     $              *conjg(phi_o(icol,ivol-nvolh)))
            endif
         enddo

      enddo

      return
      end


c=========================================================
      subroutine haction(action,h)
c  written by Massimo D'Elia
c  version 1.0 - 01/08
c  link momentum - part system action
c=========================================================

      implicit none
      include "parameters.f"

c     arguments
      real action(nvol)
      complex h(ncol,ncol,4,nvol)

c     internal variables
      integer ivol,idir
      real rtr
      real loc_action


      do ivol=1,nvol

         loc_action=0

         do idir=1,4
            call rtrace_uu(rtr,h(1,1,idir,ivol),h(1,1,idir,ivol))
            loc_action=loc_action+rtr
         enddo

         action(ivol)=action(ivol)+loc_action*0.5

      enddo


      return
      end




c=========================================================
      subroutine action(az)
c  written by Massimo D'Elia
c  version 1.0 - 01/08
c  calculate locally the total action of the system
c  IMPORTANT: the action have to be calculated locally
c=========================================================

      implicit none
      include "parameters.f"

c     arguments
      real az(nvol)

c     common blocks
      integer nran,iad_flag
      real beta_f,beta_a,rho,pi
      common/param/nran,beta_f,beta_a,rho,pi,iad_flag
      complex p_w,ipdot
      common/momenta/p_w(ncol,ncol,4,nvol),ipdot(ncol,ncol,4,nvol)
      complex phi_e,phi_o,chi_e,chi_o
      common/pseudof/phi_e(ncol,nvolh),chi_e(ncol,nvolh),phi_o(ncol
     $     ,nvolh),chi_o(ncol,nvolh)
      real remu
      real emu,ememu
      common/param4/remu,emu,ememu

c     internal variables
      integer ivol


      do ivol=1,nvol
         az(ivol)=0
      enddo

      call uaction(az,beta_f)
      call haction(az,p_w)
      call qaction_flav(az,phi_e,phi_o,emu,ememu)


      return
      end



c=========================================================
      subroutine dif_action(dif,az_2,az_1)
c  written by Massimo D'Elia
c  version 1.0 - 01/08
c  calculate the global differenxe between two locally-calculated action
c=========================================================

      implicit none
      include "parameters.f"

c     arguments
      real dif
      real az_2(nvol),az_1(nvol)

c     internal variables
      integer ivol


      dif=0
      do ivol=1,nvol
         dif=dif+az_2(ivol)-az_1(ivol)
      enddo


      return
      end
