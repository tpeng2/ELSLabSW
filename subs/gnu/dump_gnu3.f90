


      WRITE(which,'(I3)') 100 + icount -1
       string1 = 'data/zeta1'  // '_' // which(1:3)
       open(11, file = string1, status = 'unknown')
      WRITE(which,'(I3)') 100 + icount -1
       string2 = 'data/wek'  // '_' // which(1:3)
       open(12, file = string2, status = 'unknown')


!      icount = icount + 1


       do j = 11,12  !write header on files
       rewind(j)
       write(j,*) 'set pm3d map'
       write(j,*) 'splot "-"'
       enddo


       do j = 1,ny
       y =  -Ly/2. + (j-1)*dx

       do i = 1,nx
       x =  -Lx/2. + (i-1)*dx

       tmp_out(1) =  zeta1(i,j)/f0 
       tmp_out(2) =  div_ek(i,j)
       write(11,*) x, y,  tmp_out(1)
       write(12,*) x, y,  tmp_out(2)
       enddo  ! end of row


       write(11,*) ! write blank line between rows 
       write(12,*) 
       enddo  ! end of j loop

       do j = 11,12  !forces write from buffer
       flush(j) 
       enddo
       do j = 11,12
       close(j)
       enddo

