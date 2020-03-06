
       WRITE(which,'(I3)') 100 + icount -1
       string1 = 'data/u1_hp'  // '_' // which(1:3)
       open(11, file = string1, status = 'unknown')

       string2 = 'data/Uek_hp'  // '_' // which(1:3)
       open(12, file = string2, status = 'unknown')

       string3 = 'data/Uek'  // '_' // which(1:3)
       open(13, file = string3, status = 'unknown')

       string4 = 'data/u1_lp'  // '_' // which(1:3)
       open(14, file = string4, status = 'unknown')

       string5 = 'data/Psurf_bt'  // '_' // which(1:3)
       open(15, file = string5, status = 'unknown')

       string6 = 'data/Psurf_A'  // '_' // which(1:3)
       open(16, file = string6, status = 'unknown')

       string7 = 'data/Psurf_G'  // '_' // which(1:3)
       open(17, file = string7, status = 'unknown')
       
       string8 = 'data/Psurf_G1'  // '_' // which(1:3)
       open(18, file = string8, status = 'unknown')

       string9 = 'data/Phi_bt'  // '_' // which(1:3)
       open(19, file = string9, status = 'unknown')

       string10 = 'data/u1_G1'  // '_' // which(1:3)
       open(20, file = string10, status = 'unknown')

       string11 = 'data/u1_A'  // '_' // which(1:3)
       open(21, file = string11, status = 'unknown')

       string12 = 'data/u_sum'  // '_' // which(1:3)
       open(22, file = string12, status = 'unknown')

       string13 = 'data/u_diff'  // '_' // which(1:3)
       open(23, file = string13, status = 'unknown')

       string14 = 'data/Psurf_hp'  // '_' // which(1:3)
       open(24, file = string14, status = 'unknown')


!      icount = icount + 1

       array = (u1_snap - u1_filtered)*H(1)/Htot &
           & + (u2_snap - u2_filtered)*H(2)/Htot


       do j = 11,24  !write header on files
       rewind(j)
       write(j,*) 'set pm3d map'
       write(j,*) 'splot "-"'
       enddo

       do j = 1,ny
       y =  -Ly/2. + (j-1)*dx

       do i = 1,nx
       x =  -Lx/2. + (i-1)*dx

       tmp_out(1) = u1_snap(i,j) &
           &      - u1_filtered(i,j)
       tmp_out(2) = Uek_snap(i,j) &
           &      - Uek_filtered(i,j)
       tmp_out(3) = Uek_snap(i,j)
       tmp_out(4) = u1_filtered(i,j)
       tmp_out(5) = Psurf_bt(i,j)
       tmp_out(6) = Psurf_A(i,j)
       tmp_out(7) = Psurf_G(i,j)
       tmp_out(8) = Htot*tmp_out(7)/H(2)
       tmp_out(9) = tmp_out(5) - tmp_out(8) 
       tmp_out(10) = array(i,j) - H(2)*uu_G(i,j)/Htot
       tmp_out(11) = -H(2)*uu_A(i,j)/Htot
       tmp_out(12) = tmp(10) + tmp(11)
       tmp_out(13) = tmp(1) - tmp(12)
       tmp_out(14) = Psurf_hp(i,j)

       write(11,*) x, y,  tmp_out(1)
       write(12,*) x, y,  tmp_out(2)/hek
       write(13,*) x, y,  tmp_out(3)/hek
       write(14,*) x, y,  tmp_out(4)
       write(15,*) x, y,  tmp_out(5)
       write(16,*) x, y,  tmp_out(6)
       write(17,*) x, y,  tmp_out(7)
       write(18,*) x, y,  tmp_out(8)
       write(19,*) x, y,  tmp_out(9)
       write(20,*) x, y,  tmp_out(10)
       write(21,*) x, y,  tmp_out(11)
       write(22,*) x, y,  tmp_out(12)
       write(23,*) x, y,  tmp_out(13)
       write(24,*) x, y,  tmp_out(14)
       enddo  ! end of row


       write(11,*) ! write blank line between rows 
       write(12,*) ! write blank line between rows 
       write(13,*) ! write blank line between rows 
       write(14,*) ! write blank line between rows 
       write(15,*) ! write blank line between rows 
       write(16,*) ! write blank line between rows 
       write(17,*) ! write blank line between rows 
       write(18,*) ! write blank line between rows 
       write(19,*) ! write blank line between rows 
       write(20,*) ! write blank line between rows 
       write(21,*) ! write blank line between rows 
       write(22,*) ! write blank line between rows 
       write(23,*) ! write blank line between rows 
       write(24,*) ! write blank line between rows 
       enddo  ! end of j loop

       do j = 11,24  !forces write from buffer
       flush(j) 
       enddo
       do j = 11,24
       close(j)
       enddo

