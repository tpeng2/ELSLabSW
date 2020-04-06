
        WRITE(which3,'(I5)') 10000 + iftcount
        string1 = 'specdata/ufft2d_r'  // '_' // which3(1:5)
        string1i = 'specdata/ufft2d_i'  // '_' // which3(1:5)
        string2 = 'specdata/vfft2d_r'  // '_' // which3(1:5)
        string2i = 'specdata/vfft2d_i'  // '_' // which3(1:5)
        string3 = 'specdata/ftotalfft_r'  // '_' // which3(1:5)
        string3i = 'specdata/ftotalfft_i'  // '_' // which3(1:5)
        string4 = 'specdata/fagfft_i'  // '_' // which3(1:5)
        string4i = 'specdata/fagfft_i'  // '_' // which3(1:5)
        !ufft,vfft,ftotalfft,fagfft

        open(unit=9999,file=string1,access='DIRECT',&
        & form='BINARY',status='UNKNOWN',RECL=4*((nx/2+1)*ny*nz))
        write(9999,REC=1) (((real(ufft(i,j,k)),i=1,(nx/2+1)/ftsubsmprto),j=1,ny/ftsubsmprto),k=1,nz)
        close(9999)
        open(unit=9998,file=string1i,access='DIRECT',&
        & form='BINARY',status='UNKNOWN',RECL=4*((nx/2+1)*ny*nz))
        write(9998,REC=1) (((aimag(ufft(i,j,k)),i=1,(nx/2+1)/ftsubsmprto),j=1,ny/ftsubsmprto),k=1,nz)
        close(9998)

        open(unit=9997,file=string2,access='DIRECT',&
        & form='BINARY',status='UNKNOWN',RECL=4*((nx/2+1)*ny*nz))
        write(9997,REC=1) (((real(vfft(i,j,k)),i=1,(nx/2+1)/ftsubsmprto),j=1,ny/ftsubsmprto),k=1,nz)
        close(9997)
        open(unit=9996,file=string2i,access='DIRECT',&
        & form='BINARY',status='UNKNOWN',RECL=4*((nx/2+1)*ny*nz))
        write(9996,REC=1) (((aimag(vfft(i,j,k)),i=1,(nx/2+1)/ftsubsmprto),j=1,ny/ftsubsmprto),k=1,nz)
        close(9996)

        open(unit=9995,file=string3,access='DIRECT',&
        & form='BINARY',status='UNKNOWN',RECL=4*((nx/2+1)*ny))
        write(9995,REC=1) ((real(ftotalfft(i,j)),i=1,(nx/2+1)/ftsubsmprto),j=1,ny/ftsubsmprto)
        close(9995)
        open(unit=9994,file=string3i,access='DIRECT',&
        & form='BINARY',status='UNKNOWN',RECL=4*((nx/2+1)*ny))
        write(9994,REC=1) ((aimag(ftotalfft(i,j)),i=1,(nx/2+1)/ftsubsmprto),j=1,ny/ftsubsmprto)
        close(9994)

        open(unit=9993,file=string4,access='DIRECT',&
        & form='BINARY',status='UNKNOWN',RECL=4*((nx/2+1)*ny))
        write(9993,REC=1) ((real(fagfft(i,j)),i=1,(nx/2+1)/ftsubsmprto),j=1,ny/ftsubsmprto)
        close(9993)
        open(unit=9992,file=string4i,access='DIRECT',&
        & form='BINARY',status='UNKNOWN',RECL=4*((nx/2+1)*ny))
        write(9992,REC=1) ((aimag(fagfft(i,j)),i=1,(nx/2+1)/ftsubsmprto),j=1,ny/ftsubsmprto)
        close(9992)

