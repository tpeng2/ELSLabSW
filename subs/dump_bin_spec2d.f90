    WRITE(which3,'(I6)') 100000 + iftcount
    
    ! U FFT
    string1 = 'specdata/ufft2d_r'  // '_' // trim(which3)
    string1i = 'specdata/ufft2d_i'  // '_' // trim(which3)
    open(unit=9999,file=string1,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
    write(9999,REC=1) ((real(ufft(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,nz)
    close(9999)
    open(unit=9998,file=string1i,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
    write(9998,REC=1) ((aimag(ufft(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,nz)
    close(9998)
    
    ! V FFT
    string2 = 'specdata/vfft2d_r'  // '_' // trim(which3)
    string2i = 'specdata/vfft2d_i'  // '_' // trim(which3)
    open(unit=9997,file=string2,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
    write(9997,REC=1) ((real(vfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,nz)
    close(9997)
    open(unit=9996,file=string2i,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
    write(9996,REC=1) ((aimag(vfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,nz)
    close(9996)

    ! Forcing FFT
    string3 = 'specdata/ftotalfft_r'  // '_' // trim(which3)
    string3i = 'specdata/ftotalfft_i'  // '_' // trim(which3)
    open(unit=9995,file=string3,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol)
    write(9995,REC=1) (real(ftotalfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2))),inkrow=1,szftrdrow)
    close(9995)
    open(unit=9994,file=string3i,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol)
    write(9994,REC=1) (aimag(ftotalfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2))),inkrow=1,szftrdrow)
    close(9994)

    ! Forcing-AG FFT
    string4 = 'specdata/f_agfft_r'  // '_' // trim(which3)
    string4i = 'specdata/f_agfft_i'  // '_' // trim(which3)
    open(unit=9993,file=string4,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol)
    write(9993,REC=1) (real(fagfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2))),inkrow=1,szftrdrow)
    close(9993)
    open(unit=9992,file=string4i,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol)
    write(9992,REC=1) (aimag(fagfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2))),inkrow=1,szftrdrow)
    close(9992)

    ! U-AG (BC) FFT
    string5 = 'specdata/u_agfft_bc_r'  // '_' // trim(which3)
    string5i = 'specdata/u_agfft_bc_i'  // '_' // trim(which3)
    open(unit=9991,file=string5,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol)
    write(9991,REC=1) (real(u_agfft_bc(iftsubkl(inkrow,1),iftsubkl(inkrow,2))),inkrow=1,szftrdrow)
    close(9991)
    open(unit=9990,file=string5i,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol)
    write(9990,REC=1) (aimag(u_agfft_bc(iftsubkl(inkrow,1),iftsubkl(inkrow,2))),inkrow=1,szftrdrow)
    close(9990)

    ! U-AG (BC) FFT
    string6 = 'specdata/v_agfft_bc_r'  // '_' // trim(which3)
    string6i = 'specdata/v_agfft_bc_i'  // '_' // trim(which3)
    open(unit=9989,file=string6,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol)
    write(9989,REC=1) (real(v_agfft_bc(iftsubkl(inkrow,1),iftsubkl(inkrow,2))),inkrow=1,szftrdrow)
    close(9989)
    open(unit=9988,file=string6i,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol)
    write(9988,REC=1) (aimag(v_agfft_bc(iftsubkl(inkrow,1),iftsubkl(inkrow,2))),inkrow=1,szftrdrow)
    close(9988)

    ! U-AG (POINCARE) FFT
    string7 = 'specdata/u_agfft_p_r'  // '_' // trim(which3)
    string7i = 'specdata/u_agfft_p_i'  // '_' // trim(which3)
    open(unit=9987,file=string7,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
    write(9987,REC=1) ((real(u_agfft_p(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,2)
    close(9987)
    open(unit=9986,file=string7i,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol)
    write(9986,REC=1) ((aimag(u_agfft_p(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,2)
    close(9986)

    ! V-AG (POINCARE) FFT
    string8 = 'specdata/v_agfft_p_r'  // '_' // trim(which3)
    string8i = 'specdata/v_agfft_p_i'  // '_' // trim(which3)
    open(unit=9985,file=string8,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
    write(9985,REC=1) ((real(v_agfft_p(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,2)
    close(9985)
    open(unit=9984,file=string8i,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
    write(9984,REC=1) ((aimag(v_agfft_p(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,2)
    close(9984)

    ! U-G FFT
    string9 = 'specdata/u_qgfft_r'  // '_' // trim(which3)
    string9i = 'specdata/u_qgfft_i'  // '_' // trim(which3)
    open(unit=9983,file=string9,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
    write(9983,REC=1) ((real(u_qgfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,nz)
    close(9983)
    open(unit=9982,file=string9i,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
    write(9982,REC=1) ((aimag(u_qgfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,nz)
    close (9982)

    ! V-G FFT
    string10 = 'specdata/v_qgfft_r'  // '_' // trim(which3)
    string10i = 'specdata/v_qgfft_i'  // '_' // trim(which3)
    open(unit=9981,file=string10,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
    write(9981,REC=1) ((real(v_qgfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,nz)
    close(9981)
    open(unit=9980,file=string10i,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
    write(9980,REC=1) ((aimag(v_qgfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,nz)
    close (9980)

    !U-AG FFT
    string11 = 'specdata/u_agfft_r'  // '_' // trim(which3)
    string11i = 'specdata/u_agfft_i'  // '_' // trim(which3)
    open(unit=9979,file=string11,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
    write(9979,REC=1) ((real(u_agfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,nz)
    close(9979)
    open(unit=9978,file=string11i,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
    write(9978,REC=1) ((aimag(u_agfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,nz)
    close (9978)

    !V-AG FFT
    string12 = 'specdata/v_agfft_r'  // '_' // trim(which3)
    string12i = 'specdata/v_agfft_i'  // '_' // trim(which3)

    open(unit=9977,file=string12,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
    write(9977,REC=1) ((real(v_agfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,nz)
    close(9977)
    open(unit=9976,file=string12i,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
    write(9976,REC=1) ((aimag(v_agfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,nz)
    close (9976)

    ! ETA FFT (Total potential energy)
    string13 = 'specdata/eta_fft_r'  // '_' // trim(which3)
    string13i = 'specdata/eta_fft_i'  // '_' // trim(which3)
    open(unit=9979,file=string13,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
    write(9979,REC=1) (real(etafft(iftsubkl(inkrow,1),iftsubkl(inkrow,2))),inkrow=1,szftrdrow)
    close(9979)
    open(unit=9978,file=string13i,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
    write(9978,REC=1)  (aimag(etafft(iftsubkl(inkrow,1),iftsubkl(inkrow,2))),inkrow=1,szftrdrow)
    close (9978)

    ! ETA-G FFT (Geostrophic part, for potential energy)
    string14 = 'specdata/eta_qgfft_r'  // '_' // trim(which3)
    string14i = 'specdata/eta_qgfft_i'  // '_' // trim(which3)
    open(unit=9977,file=string14,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
    write(9977,REC=1) (real(eta_qgfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2))),inkrow=1,szftrdrow)
    close(9977)
    open(unit=9976,file=string14i,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
    write(9976,REC=1) (aimag(eta_qgfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2))),inkrow=1,szftrdrow)
    close (9976)
