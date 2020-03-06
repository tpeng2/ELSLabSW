!      COMPUTE THE LAPLACIAN

        call dfftw_execute_dft_r2c(pr2c,datr,datc)
        datc = datc/nx/ny        
        do j = 1,ny
         do i = 1,nx/2+1
          datc(i,j) = - datc(i,j)*k2(i,j)
         enddo
        enddo
        call dfftw_execute_dft_c2r(pc2r,datc,datr)

