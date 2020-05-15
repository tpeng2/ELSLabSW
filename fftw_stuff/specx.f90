
    
    !1-D spectrum for x-direction
    CALL dfftw_execute_dft_r2c(vr2cx,datrx,datcx)

    spectrum(:) = 0.
    do i = 1,nx/2+1
       tmp(1) = datcx(i)*conjg(datcx(i))
       spectrumx(nkx1d(i)) = spectrumx(nkx1d(i))+tmp(1)
    enddo

