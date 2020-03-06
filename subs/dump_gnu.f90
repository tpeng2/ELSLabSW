
       WRITE(which,'(I3)') 100 + icount1
       string1 = 'data/w_avg'  // '_' // which(1:3)
       string2 = 'data/p_avg'  // '_' // which(1:3)
       string3 = 'data/wp_avg'  // '_' // which(1:3)
       open(14, file = string1, status = 'unknown')
       open(15, file = string2, status = 'unknown')
       open(16, file = string3, status = 'unknown')
       icount1 = icount1 + 1


       print*, 'dump gnu', its
       do j = 14,16  !write header on files
       rewind(j)
       write(j,*) 'set pm3d map'
       write(j,*) 'splot "-"'
       enddo


       do j = 1,ny
       y =  -Ly/2. + (j-1)*dx

       do i = 1,nx
       x =  -Lx/2. + (i-1)*dx
!      tmp_out(1) = Psurf(i,j)/9.8/2/dt
       tmp_out(1) = w_avg(i,j)/icount_wp
       tmp_out(2) = p_avg(i,j)/icount_wp
       tmp_out(3) = wp_avg(i,j)/icount_wp - &
           &        tmp_out(1)*tmp_out(2)
       write(14,*) x, y,  tmp_out(1)
       write(15,*) x, y,  tmp_out(2)
       write(16,*) x, y,  tmp_out(3)
       enddo  ! end of row


       write(14,*) ! write blank line between rows 
       write(15,*) 
       write(16,*) 
       enddo  ! end of j loop

       do j = 14,16  !forces write from buffer
       flush(j) 
       enddo
       do j = 14,16
       close(j)
       enddo

