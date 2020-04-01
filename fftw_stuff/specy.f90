
    
    !1-D spectrum for x-direction
    CALL dfftw_execute_dft_r2c(vr2cy,datry,datcy)

    spectrum(:) = 0.
    do j = 1,ny/2+1
       tmp(1) = datcy(j)*conjg(datcy(j))
       spectrumy(nky1d(j)) = spectrumy(nky1d(j))+tmp(1)
    enddo

