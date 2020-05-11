

       u(:,:,:,:) = 0.
       v(:,:,:,:) = 0.
       eta(:,:,:,:) = 0.
       u_ag(:,:,:) = 0.
       v_ag(:,:,:) = 0.
       eta_ag(:,:) = 0.
       u_qg(:,:,:) = 0.
       v_qg(:,:,:) = 0.
       eta_qg(:,:) = 0.
       taux(:,:) = 0.
       tauy(:,:) = 0.
       zeta(:,:) = 0.
       div(:,:) = 0.
       B(:,:) = 0.
       B_nl(:,:) = 0.
       grad2u(:,:) = 0.
       grad2v(:,:) = 0.
       grad4u(:,:) = 0.
       grad4v(:,:) = 0.
       array(:,:) = 0.
       thickness(:,:) = 0.
       rhs_eta(:,:,:) = 0.
       Uek(:,:,3) = 0.
       Vek(:,:,3) = 0.

       top(:) = 0.
       bot(:) = 0.
       top(1) = 1.
       bot(nz) = 1.

       count_specs_1 = 0
       count_specs_2 = 0
       count_specs_ek = 0
       count_specs_to = 0
       count_specs_AG = 0
       ke1_spec(:) = 0.
       ke2_spec(:) = 0.
       ke_ek_spec(:) = 0.
       for_to_spec(:) = 0.
       for_ag_spec(:) = 0.

       Htot = 4000.
       H(1) = 1000.
       H(2) = Htot - H(1)

       if (nz.eq.2) then
          gprime(1) = 0.
          gprime(2) = Htot*c_bc**2/H(1)/H(2) 
       else
          print*, 'bug; need to set gprimes'
          stop
       endif

  ! --- Wind Stress 

       f(:) = f0
       taux_steady=0.0;
       do j = 1,ny
          y = -Ly/2. + (j-1)*dy
          if(forcingtype.eq.0) then
            taux_steady(:,j) = tau0
          elseif (forcingtype.eq.1) then
            taux_steady(:,j) = tau0*cos(twopi*y/Ly)
         elseif (forcingtype.eq.2) then
            taux_steady(:,j) = tau0*(1.0/3.0*cos(twopi*y/Ly)&
            +1.0/3.0*cos(twopi*y/Ly*4.0)+1.0/3.0*cos(twopi*y/Ly*8.0))
          endif
       enddo
       open(unit=22,file='taux_steady.dat',access='sequential',form='formatted',action='write')
       
       array = taux_steady
       include 'subs/bndy.f90'
       taux_steady = array

       taux_var(:,:) = tau1
       array = taux_var
       include 'subs/bndy.f90'
       taux_var = array

       ! Output forcing field
       write(fmtstr,"(I6)") nx
       fmtstr='('//trim(fmtstr)//'e15.6)'
       do j = 0,ny+1
          write(22,fmtstr) taux_steady(:,j)
          write(*,*) '(1,j),forcing',j,taux_steady(1,j)
       enddo
       close (22)
       

!   --- Restart
       icount = 0 !for  output file index
       iftcount = 0
       ignucount = 0
       if ( ifrestart .eqv. .false. ) then
         time = 0.  !in second
         restart_from=time
         print*,'New case from',restart_from, 'day','icount,iftcount',icount,iftcount
         ! forcing
         amp_load=0.0
         amp_save=0.0
          ! add seed IC
          tmp(1) = Lx/8.
          do j = 1,ny
             y = (j-1)*dy - Ly/2.
             y = y/tmp(1)
             do i = 1,nx
                x = (i-1)*dx - Lx/2.
                x = x/tmp(1)
                tmp(2) = x**2 + y**2
                eta(i,j,2,1) =  20.0*exp(-tmp(2))
                eta(i,j,1,1) =  0.
             enddo
          enddo
       array = eta(:,:,1,1)
       include 'subs/bndy.f90'
       eta(:,:,1,1) = array

       else !if ifrestart
          open(0,file='restart')
          do j = 0,nny
          do i = 0,nnx
             read(0,*) Uek(i,j,1),Vek(i,j,1),u(i,j,1,1), &
                 &      u(i,j,2,1),v(i,j,1,1),v(i,j,2,1), &
                 &      eta(i,j,2,1)
          enddo
          enddo
          read(0,*) icount_srt,time,nspecfile,iftcount_srt
          close(0)
          if(ifsteady.eqv..false..and.iou_method.eq.1) then
            ! Read forcing data !Only when iou_method == 2
            open(1,file='restart_amp')
            read(1,*) amp_save,amp_load ! inverse order when saving these 
            close(1)
            print*, 'Amp_matrix at the end of the previous simulation is loaded',amp_load
         end if
         restart_from=time/86400
         print*, 'Restart from', restart_from, 'day'
         print*, 'FFT index begins with',iftcount_srt
         print*, 'snapshot index begins with',icount_srt

!         if (restart_from == 999 ) then
!            time = 0
!            icount = 0
!         else 
!            time = (restart_from-100)*iout* dt
!            icount = restart_from - 100
!         end if

!        print*, 'time =', time/86400. , 'days'

!         WRITE(which,'(I3)') restart_from

!         string1 = 'data/u'  // '_' // which(1:3)
!         string2 = 'data/v'  // '_' // which(1:3)
!         string3 = 'data/eta'  // '_' // which(1:3)
! 
!         open(unit=14,file=string1,access='DIRECT',&
!              & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!         read(14,REC=1) u(:,:,:,1)
!         close(14)
!         
!         open(unit=15,file=string2,access='DIRECT',&
!              & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!         read(15,REC=1) v(:,:,:,1)
!         close(15)
! 
!         open(unit=16,file=string3,access='DIRECT',&
!              & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!         read(16,REC=1) eta(:,:,:,1)
!         close(16)

!         eta(:,:,1,:) = 0.

       endif  ! --- restart

!
!    Set bndy conditions
!
       do k = 1,nz
          array(:,:) = u(:,:,k,1)
          include 'subs/bndy.f90'
          u(:,:,k,1) = array(:,:)

          array(:,:) = v(:,:,k,1)
          include 'subs/bndy.f90'
          v(:,:,k,1) = array(:,:)

          array(:,:) = eta(:,:,k,1)
          include 'subs/bndy.f90'
          eta(:,:,k,1) = array(:,:)
       enddo ! k = 1,nz

       array(:,:) = Psurf(:,:)
       include 'subs/bndy.f90'
       Psurf(:,:) = array(:,:)

