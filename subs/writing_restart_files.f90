      write(which2,'(I6)') floor(today) 
      which='_'//trim(adjustl(which2))
      which2='rst/'
      if(its==nsteps) then
        which2="./"
        which=""
      end if
      ! ========= I/O for restart files =========
      which3=trim(which2)//'ke1_spec'//trim(adjustl(which))
      open(99, file =which3 )
      do i = 1,ny
         write(99,'(1i4,1e12.4)')i,ke1_spec(i)/nx/nx/ny/ny/count_specs_1
      enddo
      close(99)

      which3=trim(which2)//'ke_ek_spec'//trim(adjustl(which))
      open(98, file = which3)
      do i = 1,ny
         write(98,'(1i4,1e12.4)') i, ke_ek_spec(i)/nx/nx/ny/ny/count_specs_ek
      enddo
      close(98)
      
      which3=trim(which2)//'ke2_spec'//trim(adjustl(which))
      open(97, file = which3)
      do i = 1,ny
         write(97,'(1i4,1e12.4)') i, ke2_spec(i)/nx/nx/ny/ny/count_specs_2
      enddo
      close(97)

      which3=trim(which2)//'for_to_spec'//trim(adjustl(which))
      open(96, file = which3)
      do i = 1,ny
         write(96,'(1i4,1e12.4)') i, for_to_spec(i)/nx/nx/ny/ny/count_specs_to
      enddo
      close(96)

      which3=trim(which2)//'for_ag_spec'//trim(adjustl(which))
      open(95, file = which3)
      do i = 1,ny
         write(95,'(1i4,1e12.4)') i, for_ag_spec(i)/nx/nx/ny/ny/count_specs_ag
      enddo
      close(95)

      which3=trim(which2)//'ke1_qg_spec'//trim(adjustl(which))
      open(94, file =which3 )
      do i = 1,ny
         write(94,'(1i4,1e12.4)')i,ke1_qg_spec(i)/nx/nx/ny/ny/count_specs_1
      enddo
      close(94)

      which3=trim(which2)//'ke2_qg_spec'//trim(adjustl(which))
      open(93, file =which3 )
      do i = 1,ny
         write(93,'(1i4,1e12.4)')i,ke2_qg_spec(i)/nx/nx/ny/ny/count_specs_1
      enddo
      close(93)

      which3=trim(which2)//'ke1_ag_spec'//trim(adjustl(which))
      open(92, file =which3 )
      do i = 1,ny
         write(92,'(1i4,1e12.4)')i,ke1_ag_spec(i)/nx/nx/ny/ny/count_specs_1
      enddo
      close(92)

      which3=trim(which2)//'ke2_ag_spec'//trim(adjustl(which))
      open(91, file =which3 )
      do i = 1,ny
         write(91,'(1i4,1e12.4)')i,ke2_ag_spec(i)/nx/nx/ny/ny/count_specs_1
      enddo
      close(91)

      which3=trim(which2)//'restart'//trim(adjustl(which))
      open(0,file=which3)
      array = Uek(:,:,3)
      include 'subs/bndy.f90'
      Uek(:,:,3) = array
      array = Vek(:,:,3)
      include 'subs/bndy.f90'
      Vek(:,:,3) = array
      array = u(:,:,1,3)
      include 'subs/bndy.f90'
      u(:,:,1,3) = array
      array = u(:,:,2,3)
      include 'subs/bndy.f90'
      u(:,:,2,3) = array
      array = v(:,:,1,3)
      include 'subs/bndy.f90'
      v(:,:,1,3) = array
      array = v(:,:,2,3)
      include 'subs/bndy.f90'
      v(:,:,2,3) = array
      array = eta(:,:,2,3)
      include 'subs/bndy.f90'
      eta(:,:,2,3) = array
      do j = 0,nny
      do i = 0,nnx
         write(0,*) Uek(i,j,3),Vek(i,j,3),u(i,j,1,3), &
            &      u(i,j,2,3),v(i,j,1,3),v(i,j,2,3), &
            &      eta(i,j,2,3)
      enddo
      enddo
      write(0,*) icount+icount_srt,time,nspecfile,iftcount+iftcount_srt
      close(0)
       
      if(its==nsteps)then
        write(*,*) 'restart files are saved for the final step'
      else 
        write(*,*) 'restart files are saved for day', today
        write(*,*) 'restart files are saved at',which2
      end if
