! Initialize forcing
  !forcing
   Fws=Fs_high ! highest frequency: f=30
   allocate(omega(n_omega)) !n_omega defined in parameters.f90
   allocate(A_n(n_omega))
   allocate(phi(n_omega))
   allocate(amp_matrix_rand(nsteps))
   write(*,*) 'omega, A_n, and phi are allocated'
   delta_omega=2*pi/86400.*(Fs_high-Fs_low)/n_omega !incremental frequency step (linearly sampled from Fs_low to Fws)
   !construct Amp frequency curve
   amp = 1.0
   ! Transient forcing is a time series with multiple frequencies and random phases
   ! I.e., \tau(x,y,t)=\tau_0(x,y)+\tau_t, where \tau_t=Sum(A_n(omega(n))*sin(omega(n)*t+phase))
   ! amplitude A_n is a function of frequencies
   ! irand_proc=0, Ornstein-Uhlenbeck; 1, Gaussian;
   ! For example: 
   !  a) Ornstein-Uhlenbeck amplitude: 
   !     ==> A_n(omega)= const. before cut-off; or, prop to omega^-1 after cut-off.
   !  b) Gaussian:
   !     ==> A_n(omega)= sqrt(N_n) exp((omega-f_0)**2/(2*sigma)**2)
   do i = 1, n_omega !1e*5
      omega(i) = 2*pi/86400.*(Fs_low)+i*delta_omega !doesn't need to have a very low frequency
      phi(i)= ran2(iseed)*2.0*pi
      if(irand_proc==0) then
         if(omega(i).le.f_thrhld) then
               A_n(i) = amp
         else
               A_n(i) = amp*(f_thrhld/omega(i))
         endif
      elseif (irand_proc==1) then
         A_n(i)=1/sigma_An/sqrt(2*pi)*exp(-((omega(i)-f_thrhld)/(sigma_An))**2.0/2.0)
      elseif (irand_proc==2) then !Fs_high here is a dummy variable
         if(omega(i).le.Fs_low) then
            A_n(i) = 1e-15
         elseif (omega(i).gt.Fs_low .and.omega(i).le.f_thrhld ) then
            A_n(i) = amp
         elseif (omega(i).gt.f_thrhld ) then
            A_n(i) = amp*(f_thrhld/omega(i))
         endif
      end if !irand_proc
   enddo

   open(unit=21,file='forcing_frequency_dist.dat',access='sequential',form='formatted',action='write')
  
   do i = 1, n_omega !1e*5
      write(21,'(3e15.6)') omega(i)/f0,A_n(i),phi(i)
   end do
   
   close(21)
   
   ! Read ampmatrix
   write(*,*) 'Amp frequency chart is constructed'
   
   if ( ifsteady .eqv. .true. ) then
    !! STEADY  
      write(*,*)'steady forcing'
      amp_matrix(:) = 0.
   elseif ( ifsteady .eqv. .false. ) then
    !! UNSTEADY 
   !<<< istoc_method==0 <<< read amp_matrix (must copy to the case direcotry)
   !<<< istoc_method==1 <<< generate random number and time series for transient forcing 
      write(*,*)'unsteady forcing'
      if(istoc_method==0) then 
         write(*,*)'istoc_method=0, read Amp_matrix file'
         open(99, file = 'amp_matrix')
         write(*,*) 'amp_matrix file is opened'
         do it = 1,nsteps+1
            read(99,*) amp_matrix(it)
            if (mod(it,2001)==1) write(*,*) 'Reading Amp_matrix',it,'th step'
         enddo
         amp_matrix(:) = 2.0*amp_matrix(:)
         write(*,*) 'amp_matrix file is read for', nsteps+1, 'steps'
         close(99)
      else if(istoc_method==1) then
         write(*,*)'istoc_method=1, generate Amp_matrix array before simulation'
         do itt = 1, nsteps
            time_tmp(itt) = (itt-1)*dt
            amp_forcing = 0.
            do i = 1, n_omega !1e*5
            amp_forcing=amp_forcing+(A_n(i)*sin(omega(i)*time_tmp(itt)+phi(i))) !over omega (frequencies)
            end do
            amp_matrix(itt)=amp_forcing
         enddo
         ampfactor=sum(amp_matrix)/nsteps !temporary stored as a mean
         rms_amp=0.0 !initialize
         do itt = 1,nsteps
            rms_amp = rms_amp + (amp_matrix(itt)-ampfactor)**2
         enddo
         rms_amp = rms_amp/nsteps !normalized by timesteps
         rms_amp = sqrt(rms_amp) !sine wave RMS
         ampfactor = c_sigma/rms_amp
         amp_matrix=amp_matrix*ampfactor*2.0;
         write(*,*) 'amp_matrix is scaled by c_sigma/rms_amp',c_sigma/rms_amp,'and with a scaling factor',2.0
         write(*,*) 'RMS of amp_matrix is',4.0*c_sigma**2.0
         open(unit=20,file='amp_matrix.dat',access='sequential',form='formatted',action='write')
         do itt = 1, nsteps
            write(20,'(2e15.6)') time_tmp(itt)/86400, amp_matrix(itt)
         enddo
         close (20)
         write(*,*) 'Amp_matrix after normalization is stored, scale factor',ampfactor*2.0
      endif
   endif
