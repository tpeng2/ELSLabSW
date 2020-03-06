
       WRITE(which,'(I3)') 100 + icount -1
       string1 = 'data/p_lp' !  // '_' // which(1:3)
       open(11, file = string1, status = 'unknown')

       string2 = 'data/p_hp'  ! // '_' // which(1:3)
       open(12, file = string2, status = 'unknown')

       string3 = 'data/p_A1'   !// '_' // which(1:3)
       open(13, file = string3, status = 'unknown')

       string4 = 'data/p_-G'   !// '_' // which(1:3)
       open(14, file = string4, status = 'unknown')

       string5 = 'data/p_hp-A1'  ! // '_' // which(1:3)
       open(15, file = string5, status = 'unknown')

       string6 = 'data/p_hp-A1+G'  ! // '_' // which(1:3)
       open(16, file = string6, status = 'unknown')


!      icount = icount + 1


       do j = 11,16  !write header on files
       rewind(j)
       write(j,*) 'set pm3d map'
       write(j,*) 'splot "-"'
       enddo

       do j = 1,ny
       y =  -Ly/2. + (j-1)*dx

       do i = 1,nx
       x =  -Lx/2. + (i-1)*dx

       tmp_out(1) = p_filtered(i,j) 
       tmp_out(2) = Psurf_hp(i,j) 
       tmp_out(3) = -H(2)*Psurf_A(i,j)/Htot
       tmp_out(4) = -Psurf_G(i,j)
       tmp_out(5) = Psurf_hp(i,j) + H(2)*Psurf_A(i,j)/Htot
       tmp_out(6) = tmp_out(5) - tmp_out(4)

       write(11,*) x, y,  tmp_out(1)
       write(12,*) x, y,  tmp_out(2)
       write(13,*) x, y,  tmp_out(3)
       write(14,*) x, y,  tmp_out(4)
       write(15,*) x, y,  tmp_out(5)
       write(16,*) x, y,  tmp_out(6)
       enddo  ! end of row


       write(11,*) ! write blank line between rows 
       write(12,*) ! write blank line between rows 
       write(13,*) ! write blank line between rows 
       write(14,*) ! write blank line between rows 
       write(15,*) ! write blank line between rows 
       write(16,*) ! write blank line between rows 
       enddo  ! end of j loop

       do j = 11,16  !forces write from buffer
       flush(j) 
       enddo
       do j = 11,16
       close(j)
       enddo

