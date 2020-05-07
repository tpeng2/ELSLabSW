 ! Initialize forcing
  !forcing
   Fws=20 !sampling extension (beyond f0)
   f_thrhld=0.2*f0
   Omgrange=3000
   n_omega=Fws*Omgrange !20 per day to construct the curve
   allocate(omega(n_omega))
   allocate(A_n(n_omega))
   allocate(phi(n_omega))
   allocate(amp_matrix_rand(nsteps))
   write(*,*) 'omega, A_n, and phi are allocated'
   delta_omega=2*pi/86400./Omgrange
   !construct Amp frequency curve
   amp = 1.0
  
   do i = 1, n_omega !1e*5
      omega(i) = i*delta_omega
      phi(i)= ran2(iseed)*2.0*pi
      if(omega(i).le.f_thrhld) then
            A_n(i) = amp
      else
            A_n(i) = amp*(f_thrhld/omega(i))
      endif
   enddo

   open(unit=21,file='forcing_frequency_dist.dat',access='sequential',form='formatted',action='write')
  
   do i = 1, n_omega !1e*5
      write(21,'(3e15.6)') omega(i)/f0,A_n(i)**2,phi(i)
   end do
   
   close(21)
   
   ! Read ampmatrix
   write(*,*) 'Amp frequency chart is constructed'
   
   if ( ifsteady .eqv. .true. ) then
      write(*,*)'steady forcing'
      amp_matrix(:) = 0.
   elseif ( ifsteady .eqv. .false. ) then
      write(*,*)'unsteady forcing'
      if(iou_method==0) then
         write(*,*)'iou_method=0, read Amp_matrix file'
         open(99, file = 'amp_matrix')
         write(*,*) 'amp_matrix file is opened'
         do it = 1,864000
            read(99,*) amp_matrix(it)
            if (mod(it,2001)==1) write(*,*) 'Reading Amp_matrix',it,'th step'
         enddo
         amp_matrix(:) = 2.0*amp_matrix(:)
         write(*,*) 'amp_matrix file is read'
      else if(iou_method==1) then
         write(*,*)'iou_method=1, generate Amp_matrix array before simulation'
         do itt = 1, nsteps
            time = (itt-1)*dt
            amp_forcing = 0.
            amp_forcing=sum(A_n*sin(omega*time+phi))
            amp_matrix(itt)=amp_forcing
            rms_amp = rms_amp + amp_forcing**2
         enddo
         close(20)
         rms_amp = rms_amp/nsteps
         rms_amp = sqrt(rms_amp)
         ampfactor = c_sigma/rms_amp
         amp_matrix=amp_matrix*ampfactor*2.0;
         write(*,*) 'amp_matrix is devided by c_sigma/rms_amp',c_sigma/rms_amp,'and with a scaling factor',2.0
         open(unit=20,file='amp_matrix.dat',access='sequential',form='formatted',action='write')
         do itt = 1, nsteps
            write(20,'(2e15.6)') time/86400, amp_matrix(itt)
         enddo
         write(*,*) 'Amp_matrix after normalization is stored, scale factor',ampfactor*2.0
         ! call get_tau_amp_AR(time,amp_matrix(its))
         write(*,*) 'time, read first forcing step', time, amp_matrix(its)
      else if(iou_method==2) then
         write(*,*)'iou_method=2, Euler-method for stochastic Lagevin equation'
         ! the first step
         call OU_Euler(amp_load,amp_matrix(its),c_theta,c_mu,c_sigma)
      endif
   endif
