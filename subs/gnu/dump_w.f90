
       WRITE(which,'(I3)') 100 + icount
       string1 = 'data/w_total' !  // '_' // which(1:3)
       open(11, file = string1, status = 'unknown')
       string2 = 'data/w_lp'   !// '_' // which(1:3)
       open(12, file = string2, status = 'unknown')
       string3 = 'data/w_hp'   !// '_' // which(1:3)
       open(13, file = string3, status = 'unknown')


!      icount = icount + 1


       do j = 11,13  !write header on files
       rewind(j)
       write(j,*) 'set pm3d map'
       write(j,*) 'splot "-"'
       enddo


       do j = 1,ny
       y =  -Ly/2. + (j-1)*dx

       do i = 1,nx
       x =  -Lx/2. + (i-1)*dx

       tmp_out(1) = wek_snap(i,j)
       tmp_out(2) = w_filtered(i,j)
       tmp_out(3) = wek_snap(i,j)  &
           &      - w_filtered(i,j)
       write(11,*) x, y,  tmp_out(1)
       write(12,*) x, y,  tmp_out(2)
       write(13,*) x, y,  tmp_out(3)
       enddo  ! end of row


       write(11,*) ! write blank line between rows 
       write(12,*) ! write blank line between rows 
       write(13,*) ! write blank line between rows 
       enddo  ! end of j loop

       do j = 11,13  !forces write from buffer
       flush(j) 
       enddo
       do j = 11,13
       close(j)
       enddo

