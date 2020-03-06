!      Filter kx = 0, ky = 1 mode

        call dfftw_execute_dft_r2c(pr2c,datr,datc)
        datc = datc/nx/ny        
        do j = 1,ny
          do i = 1,nx/2+1
           if(nkx(i,j).eq.0) then
             if(nky(i,j).eq.1) then
              datc(i,j) = 0
             endif
             if(nky(i,j).eq.-1) then
              datc(i,j) = 0
             endif
           endif
         enddo
        enddo
        call dfftw_execute_dft_c2r(pc2r,datc,datr)

