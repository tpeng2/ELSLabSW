        WRITE(which3,'(I6)') 100000 + iftcount
        string1 = 'specdata/ufft2d_r'  // '_' // trim(which3)
        string1i = 'specdata/ufft2d_i'  // '_' // trim(which3)
        string2 = 'specdata/vfft2d_r'  // '_' // trim(which3)
        string2i = 'specdata/vfft2d_i'  // '_' // trim(which3)
        string3 = 'specdata/ftotalfft_r'  // '_' // trim(which3)
        string3i = 'specdata/ftotalfft_i'  // '_' // trim(which3)
        string4 = 'specdata/f_agfft_r'  // '_' // trim(which3)
        string4i = 'specdata/f_agfft_i'  // '_' // trim(which3)
        string5 = 'specdata/u_agfft_r'  // '_' // trim(which3)
        string5i = 'specdata/u_agfft_i'  // '_' // trim(which3)
        string6 = 'specdata/v_agfft_r'  // '_' // trim(which3)
        string6i = 'specdata/v_agfft_i'  // '_' // trim(which3)
        !ufft,vfft,ftotalfft,f_agfft

        open(unit=9999,file=string1,access='DIRECT',&
        & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
        write(9999,REC=1) ((real(ufft(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,nz)
        close(9999)
        open(unit=9998,file=string1i,access='DIRECT',&
        & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
        write(9998,REC=1) ((aimag(ufft(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,nz)
        close(9998)

        open(unit=9997,file=string2,access='DIRECT',&
        & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
        write(9997,REC=1) ((real(vfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,nz)
        close(9997)
        open(unit=9996,file=string2i,access='DIRECT',&
        & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
        write(9996,REC=1) ((aimag(vfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,nz)
        close(9996)

        ! open(unit=9995,file=string3,access='DIRECT',&
        ! & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol)
        ! write(9995,REC=1) (real(ftotalfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2))),inkrow=1,szftrdrow)
        ! close(9995)
        ! open(unit=9994,file=string3i,access='DIRECT',&
        ! & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol)
        ! write(9994,REC=1) (aimag(ftotalfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2))),inkrow=1,szftrdrow)
        ! close(9994)

        ! open(unit=9993,file=string4,access='DIRECT',&
        ! & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol)
        ! write(9993,REC=1) (real(f_agfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2))),inkrow=1,szftrdrow)
        ! close(9993)
        ! open(unit=9992,file=string4i,access='DIRECT',&
        ! & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol)
        ! write(9992,REC=1) (aimag(f_agfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2))),inkrow=1,szftrdrow)
        ! close(9992)

        open(unit=9991,file=string5,access='DIRECT',&
        & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
        write(9991,REC=1) ((real(u_agfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,nz)
        close(9991)
        open(unit=9990,file=string5i,access='DIRECT',&
        & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
        write(9990,REC=1) ((aimag(u_agfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,nz)
        close(9990)

        open(unit=9989,file=string6,access='DIRECT',&
        & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
        write(9989,REC=1) ((real(v_agfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,nz)
        close(9989)
        open(unit=9988,file=string6i,access='DIRECT',&
        & form='BINARY',status='UNKNOWN',RECL=4*szftrdrow*szftrdcol*nz)
        write(9988,REC=1) ((aimag(v_agfft(iftsubkl(inkrow,1),iftsubkl(inkrow,2),k)),inkrow=1,szftrdrow),k=1,nz)
        close(9988)
