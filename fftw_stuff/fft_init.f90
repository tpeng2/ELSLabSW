
    INCLUDE '/usr/include/fftw3.f'


    datr = 0; datc = 0
    CALL dfftw_plan_dft_r2c_2d(pr2c,Nx,Ny,datr,datc,FFTW_ESTIMATE)
    CALL dfftw_plan_dft_c2r_2d(pc2r,Nx,Ny,datc,datr,FFTW_ESTIMATE)

    include 'fftw_stuff/index.f90'
