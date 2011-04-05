c=========================================================
      subroutine test_reversibility(scale_each_step)
c  written by Sanfo
c  testa la reversibilit delle equazioni usate per la dinamica
c=========================================================
      implicit none
#include "parameters.f"

c     argomenti
      integer scale_each_step

c     variabili passate da blocchi common
      complex u
      common/field/u(ncol,ncol,4,nvol)
      complex p_w,ipdot
      common/momenta/p_w(ncol,ncol,4,nvol),ipdot(ncol,ncol,4,nvol)
      
c     variabili interne
      real dif_act01,dif_act02
      real normp_w,normu
      integer nudof,npdof,nfdof,nhdof
      parameter(nudof=4*8*nvol,npdof=nudof,nfdof=3*nvol_eo,
     $     nhdof=nudof+npdof+nfdof)
      complex ur(ncol,ncol,4,nvol)
      complex p_wr(ncol,ncol,4,nvol)

c     VARIABILI CHE CONTENGONO L'AZIONE
      real az_0(nvol),az_1(nvol),az_2(nvol)

c     crea la configurazione iniziale
      call create_momenta
      call create_phi
      call action(az_0)
      
c     copia la configurazione iniziale
      call copyconf(p_wr,p_w)
      call copyconf(ur,u)

c     va avanti,calcola l'azione,torna indietro e ricalcola l'azione
      call dinamica(scale_each_step)
      call action(az_1)
      
      call revert_momenta

      call dinamica(scale_each_step)
      call revert_momenta
      call action(az_2)

c     calcola la differenza tra le norme dei campi u e p_w e delle azioni
      call diff_norm(normp_w,p_wr,p_w)
      call diff_norm(normu,ur,u)
      call dif_action(dif_act01,az_0,az_1)
      call dif_action(dif_act02,az_0,az_2)

      write(*,*) "Azione1-azione0:",dif_act01
      write(*,*) "Azione2-azione0:",dif_act02

      write(*,*) "|P''-P|/npdof^0.5:",normp_w/real(npdof)**0.5
      write(*,*) "|U''-U|/nudof^0.5:",normu/real(nudof)**0.5
      write(*,*) "|H''-H|/nhdof^0.5:",dif_act02/real(nhdof)**0.5
      write(*,*) "|H''-H|/|H'-H|:",dif_act02/dif_act01

      return
      end



c=========================================================
      subroutine revert_momenta()
c  written by Sanfo
c  inverte il segno agli H in modo da tornare indietro
c=========================================================
      implicit real (a-h,o-z)
      implicit integer (i-n)
#include "parameters.f"

c     variabili passate con blocchi common
      complex p_w,ipdot
      common/momenta/p_w(ncol,ncol,4,nvol),ipdot(ncol,ncol,4,nvol)

      do ivol=1,nvol
         do idir=1,4
            do icol1=1,ncol
               do icol2=1,ncol
                  p_w(icol1,icol2,idir,ivol)=-p_w(icol1,icol2,idir,ivol)
               enddo
            enddo
         enddo
      enddo 
      
      return
      end


c=========================================================
      subroutine diff_norm(norm,mr,m)
c  written by Sanfo
c  calcola una fuffa norma della differenza tra il campo ottenuto
c  ritornando indietro e quello originale per la coppia di matrici
c  passate da argomenti
c=========================================================

      implicit real (a-h,o-z)
      implicit integer (i-n)
#include "parameters.f"

c     argomenti
      real norm
      complex mr(ncol,ncol,4,nvol)
      complex m(ncol,ncol,4,nvol)
      
c     variabili interne
      complex c1


      norm=0

      do ivol=1,nvol
         do idir=1,4
            do icol1=1,ncol
               do icol2=1,ncol
                  c1=mr(icol1,icol2,idir,ivol)-m(icol1,icol2,idir,ivol)
                  norm=norm+real(c1*conjg(c1))
               enddo
            enddo
         enddo
      enddo 

      norm=norm**0.5
      
      return
      end


c=========================================================
      subroutine back_campo(mr,m)
c  written by Sanfo
c  copia il campo m in mr
c=========================================================

      implicit real (a-h,o-z)
      implicit integer (i-n)
#include "parameters.f"

c     argomenti
      complex mr(ncol,ncol,4,nvol)
      complex m(ncol,ncol,4,nvol)
     

      do ivol=1,nvol
         do idir=1,4
            do icol1=1,ncol
               do icol2=1,ncol
                  mr(icol1,icol2,idir,ivol)=m(icol1,icol2,idir,ivol)
               enddo
            enddo
         enddo
      enddo 

      return
      end



c=========================================================
      subroutine test_invertibility()
c  written by Sanfo
c  testa la reversibilit delle equazioni usate per la dinamica
c=========================================================
      implicit none
#include "parameters.f"

c     variabili passate da blocchi common
      complex u
      common/field/u(ncol,ncol,4,nvol)
      complex p_w,ipdot
      common/momenta/p_w(ncol,ncol,4,nvol),ipdot(ncol,ncol,4,nvol)
      real beta,mass,mass2,residue
      common/param/beta,mass,mass2,residue
      complex potc1,potc2
      common/param3/potc1,potc2
      integer nterm
      real cost,pole,coef
      common/rhmc_a/nterm,cost,pole(nmr),coef(nmr)

c     variabili interne
      complex temp(ncol,nvol_eo)
      complex origi(ncol,nvol_eo)
      complex final(ncol,nvol_eo)
      integer i,icol,ivol
      complex c1
      real norm

c     parametri
      real sigma
      parameter(sigma=1)

      call scale_rhmc

c     estrae un vettore a caso
      do ivol=1,nvol_eo
         call gauss_vector(temp(1,ivol),sigma)
      enddo
      
c     copia il vettore estratto in quello originale
      do icol=1,ncol
         do ivol=1,nvol_eo
            origi(icol,ivol)=temp(icol,ivol)
         enddo
      enddo

c     applica l'inverter 4 volte, così ottiene temp=m^-1orig
      do i=1,4
         call multi_shift_summed_inverter
     $        (final,temp,
     $        potc1,cost,pole,coef,
     $        nterm,1)

         do icol=1,ncol
            do ivol=1,nvol_eo
               temp(icol,ivol)=final(icol,ivol)
            enddo
         enddo
      enddo
      
c     riapplica M così in teoria avrei quello originale
#ifdef mf
      call add_extf(1)
#endif
      call m2d(final,temp,potc1,mass)
#ifdef mf
      call rem_extf(1)
#endif
c     stampa il vettore originale, quello finale e quello intermedio
c     intanto calcola la norma della differenza tra l'originale ed il finale
      norm=0
      do icol=1,ncol
         do ivol=1,nvol_eo
c            write(*,*) real(origi(icol,ivol)),real(final(icol,ivol)),
c     $           real(temp(icol,ivol))
            c1=origi(icol,ivol)-final(icol,ivol)
            norm=norm+real(c1*conjg(c1))
         enddo
      enddo
      norm=norm**0.5
      norm=norm/nvol_eo/ncol
      write(*,*) "norma della differenza:",norm

      return
      end



#ifdef mf
c=====================================================================
      subroutine magnetic_field(ud)
c  written by Sanfo 
c  calculate the value of the magnetic field defined by the ue1 or ue2
c  fields for each plaquette, to check the correct definition
c=====================================================================
      implicit none
#include "parameters.f"

c     arguments
      integer ud

c     common-block passed variables
      complex ue1,ue2
      real u1phase
      common/fielde1/ue1(4,nvol)
      common/fielde2/ue2(4,nvol)
      common/expu1/u1phase(4,nvol)
      integer eotolex,lextoeo
      common/gegeo/eotolex(nvol),lextoeo(nvol)
      integer sind,coor,forw,back
      common/ge/sind(nx,ny,nz,nt),coor(nvol,4),forw(nvol,4),back(nvol,4)

c     internal variables
      integer ix,iy,iz,it
      complex caux
      real br,bi
      real mr,mi
      
      integer i1,i2,i3,i4,ivol

      mr=0
      mi=0

      do it=1,nt
         do iz=1,nz
            do iy=1,ny
               do ix=1,nx                  
                  ivol=sind(ix,iy,iz,it)
                  i1=lextoeo(ivol)
                  i2=lextoeo(forw(ivol,1))
                  i3=lextoeo(forw(forw(ivol,1),2))
                  i4=lextoeo(forw(ivol,2))

c                  write(*,*) ivol,": ",i1,i2,i3,i4

                  if(ud.eq.1) then
                     caux=ue1(1,i1)*ue1(2,i2)*conjg(ue1(1,i4)*ue1(2,i1))
                  else
                     caux=ue2(1,i1)*ue2(2,i2)*conjg(ue2(1,i4)*ue2(2,i1))
                  endif

                  br=acos(real(caux))/(2*pigr/float(nx**2))
                  bi=asin(aimag(caux))/(2*pigr/float(nx**2))
                  
                  mr=mr+br
                  mi=mi+bi
   
                  write(*,*) ivol,br,bi

               enddo
            enddo
         enddo
      enddo
      
      write(*,*) "Mean value: ",mr/nvol,mi/nvol

      return
      end
#endif
