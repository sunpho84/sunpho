      program rab

cc************************************************cc
cc Program that reads a gauge configuration       cc
cc written in binary mode by a quadrics machine   cc
cc and converts it                                cc
cc                                                cc
cc**********************************************************cc


      parameter(nodes1_x=4,nodes1_y=8,nodes1_z=8,nodes1_t=1)
      parameter(nodes2_x=2,nodes2_y=2,nodes2_z=2,nodes2_t=1)
      parameter(ltot_x =16,ltot_y =16,ltot_z =16,ltot_t =16)

      parameter(lloc1_x = ltot_x/nodes1_y,lloc1_y =ltot_y/nodes1_x,
     $          lloc1_z =ltot_z/nodes1_z,lloc1_t =ltot_t/nodes1_t)
      parameter(lloc2_x = ltot_x/nodes2_x,lloc2_y =ltot_y/nodes2_y,
     $          lloc2_z =ltot_z/nodes2_z,lloc2_t =ltot_t/nodes2_t)

      parameter (NOL1 = lloc1_x*lloc1_y*lloc1_z*lloc1_t*4)
      parameter (NOL2 = lloc2_x*lloc2_y*lloc2_z*lloc2_t*4)

      parameter(nvol_tot = ltot_x*ltot_y*ltot_z*ltot_t)
      parameter(nlink4 = 4*4*nvol_tot)
      parameter(nvols1_loc = lloc1_x*lloc1_y*lloc1_z)
      parameter(nvols2_loc = lloc2_x*lloc2_y*lloc2_z)

      complex*16 us(3,3,4,nvol_tot)
      complex*16 useo(3,3,4,nvol_tot)

      complex*16 uape1 (3,3,NOL1)
      complex*16 uape2 (3,3,NOL2)
      complex ul(3,3)
      complex uld(3,3)
      complex*16 ele
      integer par



      irecl1 = 8*3*3*2*4*nvols1_loc*lloc1_t
      irecl2 = 4*2*2*2*4*nvols2_loc*lloc2_t
      nvols1_loc_half = nvols1_loc/2
      nvols2_loc_half = nvols2_loc/2
      nvoloc = nol1/4
      nvoloch = nvoloc/2
      nvol_toth = nvol_tot/2
      open(1,file='lattice',form='unformatted',
     $     access='direct',
     $     recl=irecl1)

      open(2,file='lattice2',status='unknown')

CC READ THE CONFIGURATION AND TRANSLATE INTO SCALAR
      nrecord = 1

      do inodo_x = 0,nodes1_x - 1
         do inodo_y = 0,nodes1_y - 1
            do inodo_z = 0,nodes1_z - 1

               read(1,rec=nrecord) uape1
c               write(*,*) nrecord

               ie = 0
               io = 0
               do i_t = 1,lloc1_t
                  do i_z = 1,lloc1_z
                     do i_y = 1,lloc1_y
                        do i_x = 1,lloc1_x

                           isum = (i_x + i_y + i_z + i_t) - 
     $                            ((i_x + i_y + i_z + i_t)/2)*2
                           if(isum.eq.0) ie = ie + 1
                           if(isum.eq.1) io = io + 1
c                           write(*,*) isum,ie,io,nvoloc,nvoloch

                           i_cart = inodo_y*lloc1_x + i_x + 
     $                      (inodo_x*lloc1_y + i_y -1)*ltot_x +  
     $                      (inodo_z*lloc1_z + i_z -1)*ltot_x*ltot_y +  
     $                         (i_t-1)*ltot_x*ltot_y*ltot_z 
               
                           
                           do i_dir = 1,4
                           if(isum.eq.0) iape = ie
                           if(isum.eq.1) iape = io + nvoloch
                             
                           iape = iape + (i_dir - 1)*nvoloc
                             
                           i_cart_loc = i_dir + 
     $                        (i_x - 1)*4 + (i_y - 1)*lloc1_x*4 +
     $                         (i_z - 1)*lloc1_y*lloc1_x*4 + 
     $                         (i_t - 1)*lloc1_z*lloc1_y*lloc1_x*4      

                              do icol1 = 1,3
                                 do icol2 = 1,3
                                    us(icol1,icol2,i_dir,i_cart)  
     $                 = uape1 (icol1,icol2,iape) 
                                 enddo  !! col2
                              enddo   !! col1
                              
                              

                           enddo  !! i_dir

                        enddo !! i_x
                     enddo !! i_y
                  enddo !! i_z
               enddo  !! i_t
               
               nrecord = nrecord + 1
            enddo
         enddo
      enddo


	       ie = 0
	       io = 0
               do i_t = 1,ltot_t
                  do i_z = 1,ltot_z
                     do i_y = 1,ltot_y
                        do i_x = 1,ltot_x
                           isum = (i_x + i_y + i_z + i_t) - 
     $                            ((i_x + i_y + i_z + i_t)/2)*2
                           if(isum.eq.0) ie = ie + 1
                           if(isum.eq.1) io = io + 1
                           if(isum.eq.0) i_eo = ie
                           if(isum.eq.1) i_eo = io + nvol_toth
                           
                           i_cart = (i_x) + (i_y -1)*ltot_x + 
     $                      (i_z -1)*ltot_x*ltot_y +  
     $                          (i_t-1)*ltot_x*ltot_y*ltot_z    

	                 do icol1 = 1,3
	                 do icol2 = 1,3
	                 do idir = 1,4


	                 useo(icol1,icol2,idir,i_eo) =                   
     $	                                     us(icol1,icol2,idir,i_cart)

         	         enddo
	                 enddo
	                 enddo

	                enddo
	              enddo
                    enddo	
                  enddo	
          

                  write(2,*) 0
	 do icol1 = 1,3
            do icol2 = 1,3
               do idir = 1,4
                  do ivol = 1,nvol_tot
                     ele  = useo(icol2,icol1,idir,ivol)
                     write(2,*) dreal(ele),dimag(ele)
                  enddo
               enddo
            enddo
         enddo   



      stop
      end

      
