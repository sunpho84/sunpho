c=========================================================
      subroutine uaction(action,beta)
c  written by Massimo D'Elia
c  version 1.0 - 01/08
c  local plaquette for molecular dynamics updating
c=========================================================
      implicit none
#include "parameters.f"

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
      subroutine qaction_flav(action,phi,potc,ud)
c  written by Massimo D'Elia
c  version 1.0 - 01/08
c  local action for ud flavour
c=========================================================
      implicit none
#include "parameters.f"

c     arguments
      real action(nvol)
      complex phi(ncol,nvol_eo)
      complex potc
      integer ud
      
c     common-block passed variables
      integer nterm
      real cost,pole,coef
      common/rhmc_a/nterm,cost,pole(nmr),coef(nmr)
      
c     internal variables
      integer ivol,icol
      complex chi(ncol,nvol_eo)
      real loc_action

      call multi_shift_summed_inverter
     $     (chi,phi,
     $     potc,
     $     cost,pole,coef,
     $     nterm,ud)

      do ivol=1,nvol_eo
         loc_action=0
         do icol=1,ncol
            loc_action=
     $           loc_action+real(chi(icol,ivol)*conjg(phi(icol,ivol)))
         enddo
         action(ivol)=action(ivol)+loc_action
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
#include "parameters.f"

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
#include "parameters.f"

c     arguments
      real az(nvol)

c     common-block passed variables
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue
      complex p_w,ipdot
      common/momenta/p_w(ncol,ncol,4,nvol),ipdot(ncol,ncol,4,nvol)
      complex phi1,phi2,chi1,chi2
      common/pseudof1/phi1(ncol,nvol_eo),chi1(ncol,nvol_eo)
      common/pseudof2/phi2(ncol,nvol_eo),chi2(ncol,nvol_eo)
      complex potc1,potc2
      common/param3/potc1,potc2

      integer ivol

      do ivol=1,nvol
         az(ivol)=0
      enddo

      call uaction(az,beta)
      call haction(az,p_w)
      call qaction_flav(az,phi1,potc1,1)
      call qaction_flav(az,phi2,potc2,2)


      return
      end



c=========================================================
      subroutine dif_action(dif,az_2,az_1)
c  written by Massimo D'Elia
c  version 1.0 - 01/08
c  calculate the global differenxe between two locally-calculated action
c=========================================================
      implicit none
#include "parameters.f"

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
