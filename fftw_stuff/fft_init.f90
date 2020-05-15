    include 'fftw_stuff/fftw3.f'

    datr = 0; datc = 0
    CALL dfftw_plan_dft_r2c_2d(pr2c,Nx,Ny,datr,datc,FFTW_ESTIMATE)
    CALL dfftw_plan_dft_c2r_2d(pc2r,Nx,Ny,datc,datr,FFTW_ESTIMATE)
    write(*,*) 'Check Nx,Ny',Nx,Ny
    ! 1-D FFT
    CALL dfftw_plan_dft_r2c_1d(vr2cx,Nx,datrx,datcx,FFTW_ESTIMATE)
    CALL dfftw_plan_dft_r2c_1d(vr2cy,Ny,datry,datcy,FFTW_ESTIMATE)
    CALL dfftw_plan_dft_c2r_1d(vc2rx,Nx,datcx,datrx,FFTW_ESTIMATE)
    CALL dfftw_plan_dft_c2r_1d(vc2ry,Ny,datcy,datry,FFTW_ESTIMATE)

    include 'fftw_stuff/index.f90'
