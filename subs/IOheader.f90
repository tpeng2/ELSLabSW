!     setups for writing spectra
               !write 1D table heat
          WRITE(spectbszx,'(I3)') size(kex1_spec)+1
          specform='('//trim(spectbszx)//'e13.5)'
          write(*,*) 'spectra-x time series length+1=',spectbszx
          open (unit=999, file = 'specdata/ke1x.dat', form = 'formatted',access = 'sequential')
          write(999,specform) 0.0,kx1d
          close(999)

          WRITE(spectbszy,'(I3)') size(key1_spec)+1
          specform='('//trim(spectbszy)//'e13.5)'
          write(*,*) 'spectra-y time series length+1=',spectbszy
          open (unit=998, file = 'specdata/ke1y.dat', form = 'formatted',access = 'sequential')
          write(998,specform) 0.0,ky1d
          close(998)

          WRITE(spectbszx,'(I3)') size(kex1_spec)+1
          specform='('//trim(spectbszx)//'e13.5)'
          write(*,*) 'spectra-x time series length+1=',spectbszx
          open (unit=997, file = 'specdata/ke2x.dat', form = 'formatted',access = 'sequential')
          write(997,specform) 0.0,kx1d
          close(997)

          WRITE(spectbszy,'(I3)') size(key1_spec)+1
          specform='('//trim(spectbszy)//'e13.5)'
          write(*,*) 'spectra-y time series length+1=',spectbszy
          open (unit=996, file = 'specdata/ke2y.dat', form = 'formatted',access = 'sequential')
          write(996,specform) 0.0,ky1d
          close(996)

          WRITE(spectbszx,'(I3)') size(kex1_spec)+1
          specform='('//trim(spectbszx)//'e13.5)'
          write(*,*) 'spectra-x time series length+1=',spectbszx
          open (unit=995, file = 'specdata/ke_ekx.dat', form = 'formatted',access = 'sequential')
          write(995,specform) 0.0,kx1d
          close(995)

          WRITE(spectbszy,'(I3)') size(key1_spec)+1
          specform='('//trim(spectbszy)//'e13.5)'
          write(*,*) 'spectra-y time series length+1=',spectbszy
          open (unit=994, file = 'specdata/ke_eky.dat', form = 'formatted',access = 'sequential')
          write(994,specform) 0.0,ky1d
          close(994)
