
       WRITE(which,'(I3)') 100 + icount1
       string1 = 'data/Psurf'  // '_' // which(1:3)
       open(14, file = string1, status = 'unknown')


       print*, 'dump gnu', its
       do j = 14,14  !write header on files
       rewind(j)
       write(j,*) 'set pm3d map'
       write(j,*) 'splot "-"'
       enddo


       do j = 1,ny
       y =  -Ly/2. + (j-1)*dx

       do i = 1,nx
       x =  -Lx/2. + (i-1)*dx
       tmp_out(1) = Psurf(i,j)
       write(14,*) x, y,  tmp_out(1)
       enddo  ! end of row


       write(14,*) ! write blank line between rows 
       enddo  ! end of j loop

       do j = 14,14  !forces write from buffer
       flush(j) 
       enddo
       do j = 14,14
       close(j)
       enddo

