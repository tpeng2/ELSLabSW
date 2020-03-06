!      COMPUTE THE INVERSE LAPLACIAN

        call dfftw_execute_dft_r2c(pr2c,datr,datc)
        datc = datc/nx/ny        
        do j = 1,ny
         do i = 1,nx/2+1
          if ( k2(i,j) /= 0 )  then
           datc(i,j) = - datc(i,j)/(k2(i,j)+Fmode(k))
          else
           datc(i,j) = 0       
          endif
         enddo
        enddo
        call dfftw_execute_dft_c2r(pc2r,datc,datr)
