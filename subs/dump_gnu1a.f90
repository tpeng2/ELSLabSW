
       WRITE(which,'(I6)') 100000 + ignucount
       string1 = 'data/wek'  // '_' // which(1:6)
       open(11, file = string1, status = 'unknown')


       ignucount = ignucount + 1


       do j = 11,11  !write header on files
       rewind(j)
       write(j,*) 'set pm3d map'
       write(j,*) 'splot "-"'
       enddo


       do j = 1,ny
       y =  -Ly/2. + (j-1)*dx

       do i = 1,nx
       x =  -Lx/2. + (i-1)*dx
       tmp_out(1) = div_ek(i,j)
       write(11,*) x, y,  tmp_out(1)
       enddo  
       write(11,*) 
       enddo  

       do j = 11,11  !forces write from buffer
       flush(j) 
       enddo
       do j = 11,11
       close(j)
       enddo

