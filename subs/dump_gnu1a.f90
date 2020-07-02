
       WRITE(which,'(I6)') 100000 + icount

       string1 = 'data/wek'  // '_' // trim(which)//'.gplt'
       string2 = 'data/psurf'  // '_' // trim(which)//'.gplt'
       string3 = 'data/psurf_1st'  // '_' // trim(which)//'.gplt'
       string4 = 'data/psurf_2nd'  // '_' // trim(which)//'.gplt'
       string5 = 'data/eta_qg'  // '_' // trim(which)//'.gplt'

       open(11, file = string1, status = 'unknown')
       open(12, file = string2, status = 'unknown')
       open(13, file = string3, status = 'unknown')
       open(14, file = string4, status = 'unknown')
       open(15, file = string5, status = 'unknown')


       icount = icount + 1


       do j = 11,14  !write header on files
       rewind(j)
       write(j,*) 'set pm3d map'
       write(j,*) 'splot "-"'
       enddo


       do j = 1,ny
       y =  -Ly/2. + (j-1)*dx

       do i = 1,nx
       x =  -Lx/2. + (i-1)*dx
       tmp_out(1) = div_ek(i,j)
       tmp_out(2) = Psurf_csum(i,j)       
       tmp_out(3) = Psurf_csum(i,j)-Psurf(i,j)
       tmp_out(4) = Psurf(i,j)
       tmp_out(5) = eta_qg(i,j)
       write(11,*) x, y,  tmp_out(1)
       write(12,*) x, y,  tmp_out(2)       
       write(13,*) x, y,  tmp_out(3)
       write(14,*) x, y,  tmp_out(4)
       write(15,*) x, y,  tmp_out(5)
       enddo  
       write(11,*) 
       write(12,*)        
       write(13,*) 
       write(14,*) 
       write(15,*) 
       enddo  

       do j = 11,15  !forces write from buffer
       flush(j) 
       enddo

       do j = 11,15
       close(j)
       enddo

