
       program forcing
     
       parameter (idays  = 3000, idt  = 300)
       parameter (n_length = 86400*idays/idt)! want 86400/idt = integer
       parameter (n_omega = 10*idays)
       real omega(n_omega), phase(n_omega)
       real A_n(n_omega)
       real amp_matrix(n_length) 
       real ran2
       integer iseed,values(8)
       call date_and_time(VALUES=values)
       iseed = -(values(8)+values(7)+values(6))

      !use ifport
      print*, n_length, n_omega
       f = 7.e-5
       dt = float(idt)
       pi = 2.*asin(1.)
       delta_omega =  2.*pi/86400./idays
       amp = 1.
       print*,'max f-day', delta_omega*n_omega/f

       do i = 1, n_omega !10*idays
          omega(i) = i*delta_omega
          phase(i) = 2.*pi*ran2(iseed)
!         if(omega(i).le.0.1*f) then
          if(omega(i).le.0.2*f) then
              A_n(i) = amp
          else
!             A_n(i) = amp*(0.1*f/omega(i))
              A_n(i) = amp*(0.2*f/omega(i))
          endif
          write(21,*) omega(i)/f, A_n(i)**2
       enddo
       print*,'A_n and phase are constructed'

       ! below is for truncating the high frequencies
       ! I don't think we have to redo this stuff
       ! the factor = 1.02309827E-02 was because I wanted to
       ! truncate the series without renormalizing it.
!      do i = 1, n_omega
!         if(omega(i).ge.0.5*f) then  
!             A_n(i) = 0.
!         endif
!      enddo

!      aa = 0. This is just to compare with another way
       ! to calc the time series using random numbers
       ! i,e, an Orhstein-Ulenbeck or AR1 process
       rms_amp = 0.
       do it = 1, n_length
          time = (it-1)*dt
          amp_forcing = 0.
       !    do i = 1, n_omega
       !       phi = omega(i)*time + phase(i)
       !       amp_forcing = amp_forcing + A_n(i)*sin(phi)
       !    enddo
          amp_forcing=sum(A_n*sin(omega*time+phase))
          amp_matrix(it)=amp_forcing
          write(20,*) time/86400, amp_forcing
          rms_amp = rms_amp + amp_forcing**2
!         aa = aa*(1.-0.2*f*dt) + 0.1*sqrt(dt)*(rand()-0.5)
!         write(22,*) t/86400, aa
       enddo

       rms_amp = rms_amp/n_length
       rms_amp = sqrt(rms_amp)
       factor = 0.1/rms_amp
       print*, 'factor' , factor
!      factor = 1.02309827E-02
!      open(99, file = 'amp_matrix_truncated')
       open(99, file = 'amp_matrix_testsum')
       do it = 1, n_length 
       write(99,*) amp_matrix(it)*factor
       enddo

       print*, 'n_length = ', n_length

       end






       !SCRAP

!         if(omega(i).le.0.02*f) then
!             A_n(i) = omega(i)/(amp*0.02*f)
!         elseif(omega(i).le.0.1*f) then
!             A_n(i) = amp
!         else
!             A_n(i) = amp*(0.1*f/omega(i))
!         endif
       !  write(21,*) omega(i)/f, omega(i)*A_n(i)**2

!      do i = 1, n_omega
!         if(omega(i).gt.4.*f) A_n(i) = 0.
!      enddo

function ran2(idum)
       use data_initial
       integer :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
       real :: ran2,AM,EPS,RNMX
       PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,&
              IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,&
              IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
       INTEGER :: idum2,j,k,iv(NTAB),iy
       SAVE iv,iy,idum2
       DATA idum2/123456789/,iv/NTAB*0/,iy/0/

       if (idum .le. 0) then
              idum=max(-idum,1)
              idum2 = idum
              do j = NTAB+8,1,-1
                     k=idum/IQ1
                     idum=IA1*(idum-k*IQ1)-k*IR1
                     if (idum .lt. 0) idum=idum+IM1
                     if (j .le. NTAB) iv(j) = idum
              end do
              iy=iv(1)
       end if
       k=idum/IQ1
       idum=IA1*(idum-k*IQ1)-k*IR1
       if (idum .lt. 0) idum=idum+IM1
       k=idum2/IQ2
       idum2=IA2*(idum2-k*IQ2)-k*IR2
       if (idum2 .lt. 0) idum2=idum2+IM2
       j = 1+iy/NDIV
       iy = iv(j) - idum2
       iv(j) = idum
       if (iy .lt. 1) iy = iy+IMM1
       ran2=min(AM*iy,RNMX)
       return
end function ran2