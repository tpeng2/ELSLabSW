      WRITE(which2,'(I5)') 10000 + nspecfile

      write(*,*) 'its,itlocal,itsrow,nspecfile',its,itlocal,itsrow,nspecfile
      pathspects='specdata/ke1x_'//trim(which2)//'.dat'
      open (unit=999, file = pathspects, form = 'formatted', access = 'sequential')
         do i=1,ntsrow
            write(999,specform) tstime(i),kex1_spec_tb(:,i)/nx/ny
         end do
      close (999)

      pathspects='specdata/ke1y_'//trim(which2)//'.dat'
      open (unit=998, file = pathspects, form = 'formatted', access = 'sequential')
         do i=1,ntsrow
            write(998,specform) tstime(i),key1_spec_tb(:,i)/nx/ny
         end do
      close (998)

      pathspects='specdata/ke2x_'//trim(which2)//'.dat'
      open (unit=997, file = pathspects, form = 'formatted', access = 'sequential')
         do i=1,ntsrow
            write(997,specform) tstime(i),kex2_spec_tb(:,i)/nx/ny
         end do
      close (997)

      pathspects='specdata/ke2y_'//trim(which2)//'.dat'
      open (unit=996, file = pathspects, form = 'formatted', access = 'sequential')
         do i=1,ntsrow
            write(996,specform) tstime(i),key2_spec_tb(:,i)/nx/ny
         end do
      close (996)

      pathspects='specdata/ke_ekx_'//trim(which2)//'.dat'
      open (unit=995, file = pathspects, form = 'formatted', access = 'sequential')
         do i=1,ntsrow
            write(995,specform) tstime(i),kex_ek_spec_tb(:,i)/nx/ny
         end do
      close (995)

      pathspects='specdata/ke_eky_'//trim(which2)//'.dat'
      open (unit=994, file = pathspects, form = 'formatted', access = 'sequential')
         do i=1,ntsrow
            write(994,specform) tstime(i),key_ek_spec_tb(:,i)/nx/ny
         end do
      close (994)