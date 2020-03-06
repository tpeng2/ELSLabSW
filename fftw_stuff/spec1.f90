
    CALL dfftw_execute_dft_r2c(pr2c,datr,datc)

    spectrum(:) = 0.
    do j = 1,ny
    do i = 1,nx/2 + 1
       tmp(1) = datc(i,j)*conjg(datc(i,j))
       spectrum(nk(i,j)) = spectrum(nk(i,j))+tmp(1)
    enddo
    enddo

