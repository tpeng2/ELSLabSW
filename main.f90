
       module data_initial
       integer nx, ny, nz, nnx, nny
       integer i_diags
       double precision pi, twopi, Lx, Ly, dx, dy
       real f0, beta, r_drag, Ah, r_invLap, rf
       real tau0, tau1, hek
       integer nsteps,fileperday
       real ndays,totaltime,dt
       real restart_from
       integer subsmpstep,itape,ispechst,iout,itlocal,itsrow,ntsrow,nspecfile
       logical restart, use_ramp, ifsteady
       include 'parameters.f90'
       parameter( ntsrow=itape/ispechst  )! how many lines for a time series table (e.g. spectrum)   
       end module data_initial

       program main

       use data_initial
  
       real u(0:nnx,0:nny,nz,3), v(0:nnx,0:nny,nz,3), eta(0:nnx,0:nny,nz,3)
       real u_ag(0:nnx,0:nny,nz), v_ag(0:nnx,0:nny,nz)
       real u_qg(0:nnx,0:nny,nz), v_qg(0:nnx,0:nny,nz)
       real eta_ag(0:nnx,0:nny), eta_qg(0:nnx,0:nny)
       real Psurf(0:nnx,0:nny), rhs_Psurf(0:nnx,0:nny)
       real div(0:nnx,0:nny), zeta(0:nnx,0:nny)
       real B(0:nnx,0:nny), B_nl(0:nnx,0:nny)
       real div1(0:nnx,0:nny), zeta1(0:nnx,0:nny)
       real div2(0:nnx,0:nny), zeta2(0:nnx,0:nny)
       real div_ek(0:nnx,0:nny), zeta_ek(0:nnx,0:nny) 
       real div_ek_qg(0:nnx,0:nny)
       real grad2u(0:nnx,0:nny), grad4u(0:nnx,0:nny)
       real grad2v(0:nnx,0:nny), grad4v(0:nnx,0:nny)
       real dissi_u(0:nnx,0:nny), dissi_v(0:nnx,0:nny)
       real taux(0:nnx,0:nny), tauy(0:nnx,0:nny)
       real taux_steady(0:nnx,0:nny), taux_var(0:nnx,0:nny)
       real rhs_u(0:nnx,0:nny,nz), rhs_v(0:nnx,0:nny,nz)
       real rhs_eta(0:nnx,0:nny,nz)
       real Uek(0:nnx,0:nny,3), Vek(0:nnx,0:nny,3)
       real Uek_qg(0:nnx,0:nny), Vek_qg(0:nnx,0:nny)
       real rhs_Uek(0:nnx,0:nny), rhs_Vek(0:nnx,0:nny)
       real uu(0:nnx,0:nny), vv(0:nnx,0:nny)
       real uu1(0:nnx,0:nny), vv1(0:nnx,0:nny)
       real uu_old(0:nnx,0:nny), vv_old(0:nnx,0:nny)
       real invLap_u(0:nnx,0:nny), invLap_v(0:nnx,0:nny)
       real uh(0:nnx,0:nny), vh(0:nnx,0:nny)
       real pressure(0:nnx,0:nny), thickness(0:nnx,0:nny)
       real forcing_qg(0:nnx,0:nny), forcing_ag(0:nnx,0:nny)
       real forcing_total(0:nnx,0:nny)
       real tmp(0:10), array(0:nnx,0:nny)
       real q(0:nnx,0:nny,nz), psi(0:nnx,0:nny,nz)
       real qmode(0:nnx,0:nny,nz), psimode(0:nnx,0:nny,nz)
       real u_out(1:(nx/subsmpstep),1:(ny/subsmpstep),nz)
       real v_out(1:(nx/subsmpstep),1:(ny/subsmpstep),nz)
       real eta_out(1:(nx/subsmpstep),1:(ny/subsmpstep),nz)
       real div_ek_out(1:(nx/subsmpstep),1:(ny/subsmpstep))
       real f(0:nny)
       real gprime(nz), Htot, H(nz)
       real top(nz), bot(nz), Fmode(nz)
       real pdf(-100:100)
       real x, y, z, ramp, ramp0, time
       real Lrossby
       real amp_matrix(864000)
       real ke1, ke2, ke1_qg, ke2_qg, pe, pe_qg, etot, etot_qg
       real*4 tmp_out(10)
       !2-D FFT spectra
       real*4 ke1_spec(0:nx), ke2_spec(0:nx), ke_ek_spec(0:nx)
       real*4 for_to_spec(0:nx), for_ag_spec(0:nx)
       !1-D spectra
       real*4 kex1_spec(0:nx/2), kex2_spec(0:nx/2), kex_ek_spec(0:nx/2)
       real*4 key1_spec(0:ny/2), key2_spec(0:ny/2), key_ek_spec(0:ny/2)
       real*4 kex1_spec_tb(0:nx/2,1:ntsrow),key1_spec_tb(0:ny/2,1:ntsrow),tstime(1:ntsrow)
       real*4 kex2_spec_tb(0:nx/2,1:ntsrow),key2_spec_tb(0:ny/2,1:ntsrow)
       real*4 kex_ek_spec_tb(0:nx/2,1:ntsrow),key_ek_spec_tb(0:ny/2,1:ntsrow)

       integer i, j, k, ii, jj, kk, ip, im, jp, jm, kp, km, it, its, ntimes
       integer isubx(1:(nx/subsmpstep)),isuby(1:(ny/subsmpstep))
       integer icount, count_specs_1, count_specs_2, count_specs_ek
       integer count_specs_to, count_specs_ag

       character(20) scrap, which,which2, spectbszx,spectbszy,specform, pathspects
       character(20) string1, string2, string3, string4, string5
       character(20) string6, string7, string8, string9, string10
       character(20) string11, string12, string13, string14, string15
       character(20) string16, string17, string18, string19, string20
        
       include 'fftw_stuff/fft_params.f90'
       include 'fftw_stuff/fft_init.f90'
       ! === Define subsampling range
       isubx=(/(i, i=1,nx, subsmpstep)/)
       isuby=(/(i, i=1,ny, subsmpstep)/)
       ! === 
       !set icount and time based on if restart is true or false
      !  icount = 0 !for  output file index
      !  time = 0.  !in second
       Lrossby = c_bc/f0
       write(*,*) 'Lrossby over Lx/2pi = ', twopi*Lrossby/Lx
       write(*,*) 'Lrossby over dx = ', Lrossby/dx
       write(*,*) 'Lrossby  = ', Lrossby/1000. , 'km'
       write(*,*) 'total step, nstep= ', nsteps
       write(*,*) 'write output data for every',iout, 'steps, for ndays=', ndays
       write(*,*) 'one day = ', 86400/dt, 'time steps'
       write(*,*) 'Spectrum time series has', ntsrow, 'lines'
       write(*,*) 'dx = ', dx

       if ( ifsteady .eqv. .true. ) then
          amp_matrix(:) = 0.
       elseif ( ifsteady .eqv. .false. ) then
          open(99, file = 'amp_matrix')
          do it = 1,864000
             read(99,*) amp_matrix(it)
          enddo
          amp_matrix(:) = 2.0*amp_matrix(:)
       endif

!! --- Initialization !!
       include 'subs/initialize.f90'
       its=1
       itlocal=1
       nspecfile=0
       write(*,*) 'iout,ispechst,ntsrow',iout,ispechst,ntsrow
       !forcing
       taux(:,:) = taux_steady(:,:)*(1+amp_matrix(its))

!==============================================================
!
!      1st time step
!
!==============================================================


       ramp = 1.
       if ( use_ramp .eqv. .true. ) then
          tmp(1) = 10*twopi/f0/dt
          ramp0  =  1./tmp(1)
          ramp = ramp0
       endif

       uu(:,:) = Uek(:,:,1)
       vv(:,:) = Vek(:,:,1)
       uu_old(:,:) = Uek(:,:,1)
       vv_old(:,:) = Vek(:,:,1)
       uu1(:,:) = u(:,:,1,1)
       vv1(:,:) = v(:,:,1,1)
       include 'subs/rhs_ek.f90'

       pressure(:,:) =  0.
       do k = 1,nz
          uu(:,:) = u(:,:,k,1) 
          vv(:,:) = v(:,:,k,1)
          uu_old(:,:) = u(:,:,k,1) 
          vv_old(:,:) = v(:,:,k,1)
  
          if (k.eq.1) then
             thickness(:,:) = H(k) - eta(:,:,k+1,1) 
          elseif (k.eq.nz) then
             pressure(:,:) = pressure(:,:) + gprime(k)*eta(:,:,k,1) 
             thickness(:,:) = H(k) + eta(:,:,k,1)
          else
             pressure(:,:) =  pressure(:,:) + gprime(k)*eta(:,:,k,1) 
             thickness(:,:) =  H(k) + eta(:,:,k,1) - eta(:,:,k+1,1)
          endif
          Uek(:,:,2) = Uek(:,:,1) 
          Vek(:,:,2) = Vek(:,:,1) ! needed for 1st time step in rhs.f90 
          include 'subs/rhs.f90'
       enddo  ! end of the do k = 1,nz loop
       
!
!    bug here if nz \= 2
!    rhs_eta gives nz eqns for nz-1 interfaces
!    set eta(k =1) to be zero (rigid lid)
!    for nz = 2, rhs_eta(k = 2) is just d/dt of eta(k=2)
!    fo  nz > 2, need to solve a system of eqns to get from d/dt of thickness
!    to d/dt of eta.
!
       rhs_eta(:,:,1) = 0.
       u(:,:,:,2) = u(:,:,:,1) + dt*rhs_u(:,:,:)
       v(:,:,:,2) = v(:,:,:,1) + dt*rhs_v(:,:,:)
       eta(:,:,:,2) = eta(:,:,:,1) + dt*rhs_eta(:,:,:)
       Uek(:,:,2) = Uek(:,:,1) + dt*rhs_Uek(:,:)
       Vek(:,:,2) = Vek(:,:,1) + dt*rhs_Vek(:,:)
       time = dt
       its = its + 1
       taux(:,:) = taux_steady(:,:)*(1+amp_matrix(its))

!
!     need to correct u,v for surface pressure 
!
       ilevel = 2
       include 'subs/p_correction.f90'
       include 'subs/p_correction.f90'
          do k = 1,nz ! ? if we really need this
             array = eta(:,:,k,2)
             include 'subs/bndy.f90'
             eta(:,:,k,2) = array
          enddo

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

          
         !  open (998, file = 'ke1y.dat', status = 'unknown', access = 'append')
         !  write(998,*) 'ky',nky1d
         !  open (997, file = 'ke2x.dat', status = 'unknown', access = 'append')
         !  write(997,*) 'kx',nkx1d
         !  open (996, file = 'ke2y.dat', status = 'unknown', access = 'append')
         !  write(996,*) 'ky',nky1d
         !  open (995, file = 'keex.dat', status = 'unknown', access = 'append')
         !  write(995,*) 'kx',nkx1d
         !  open (994, file = 'keey.dat', status = 'unknown', access = 'append')
         !  write(994,*) 'ky',nky1d
!==============================================================      
!
!      subsequent time steps
!
!==============================================================      

       do its = 2, nsteps

         ! determine local time step !tstime
         itlocal=mod(its,itape)
         itsrow=int(itlocal/ispechst)
         if (itsrow==0.and.its.gt.ispechst)then
         itsrow=ntsrow
         end if
         if(itlocal==0) then
            itlocal = itape
            nspecfile=nspecfile+1
         end if
         if(mod(itlocal,ispechst)==0) then
            tstime(itsrow)=time/86400. !time to print
            ! write(*,*) 'its,itsrow,time[day]',its,itsrow,tstime(itsrow)
         end if
!        =================
          if ( mod(its,iout*10).eq.1 ) write(*,*) its, iout
!        =================
          if ( use_ramp .eqv. .true. ) then
             ramp =  min(1.,float(its)*ramp0)
          endif
!        =================

          uu(:,:) = Uek(:,:,2)
          vv(:,:) = Vek(:,:,2)
          uu_old(:,:) = Uek(:,:,1)
          vv_old(:,:) = Vek(:,:,1)
          uu1(:,:) = u(:,:,1,2)
          vv1(:,:) = v(:,:,1,2)
          include 'subs/rhs_ek.f90'

          pressure(:,:) =  0.
          do k = 1,nz
             uu(:,:) = u(:,:,k,2) 
             vv(:,:) = v(:,:,k,2)   
             uu_old(:,:) = u(:,:,k,1)
             vv_old(:,:) = v(:,:,k,1)

             if (k.eq.1) then
                thickness(:,:) = H(k) - eta(:,:,k+1,2) 
             elseif (k.eq.nz) then
                pressure(:,:) = pressure(:,:) + gprime(k)*eta(:,:,k,2) 
                thickness(:,:) = H(k) + eta(:,:,k,2)
             else
                pressure(:,:) = pressure(:,:) + gprime(k)*eta(:,:,k,2) 
                thickness(:,:) = H(k) + eta(:,:,k,1) - eta(:,:,k+1,2)
             endif
             include 'subs/rhs.f90' 
          enddo  ! k

!
!   same bug as in first time step if nz /= 2
!
          rhs_eta(:,:,1) = 0.
          u(:,:,:,3) = u(:,:,:,1) + 2.*dt*rhs_u(:,:,:)
          v(:,:,:,3) = v(:,:,:,1) + 2.*dt*rhs_v(:,:,:)
          eta(:,:,:,3) = eta(:,:,:,1) + 2.*dt*rhs_eta(:,:,:)
          Uek(:,:,3) = Uek(:,:,1) + 2*dt * rhs_Uek(:,:)
          Vek(:,:,3) = Vek(:,:,1) + 2*dt * rhs_Vek(:,:)
          time = its*dt
          taux(:,:) = taux_steady(:,:)*(1+amp_matrix(its+1))

!
!     stuff for surface pressure correction
!
          ilevel = 3
          include 'subs/p_correction.f90'
          include 'subs/p_correction.f90'

          do k = 1,nz
             array = eta(:,:,k,3)
             include 'subs/bndy.f90'
             eta(:,:,k,3) = array
          enddo

          array = Uek(:,:,3)
          include 'subs/bndy.f90'
          Uek(:,:,3) = array
          array = Vek(:,:,3)
          include 'subs/bndy.f90'
          Vek(:,:,3) = array

          u(:,:,:,2) = u(:,:,:,2) + rf*(u(:,:,:,1)+u(:,:,:,3)-2*u(:,:,:,2))
          v(:,:,:,2) = v(:,:,:,2) + rf*(v(:,:,:,1)+v(:,:,:,3)-2*v(:,:,:,2))
          eta(:,:,:,2) = eta(:,:,:,2)   &
              &        + rf*(eta(:,:,:,1)+eta(:,:,:,3)-2*eta(:,:,:,2))
          Uek(:,:,2) = Uek(:,:,2) + rf*(Uek(:,:,1)+Uek(:,:,3)-2*Uek(:,:,2))
          Vek(:,:,2) = Vek(:,:,2) + rf*(Vek(:,:,1)+Vek(:,:,3)-2*Vek(:,:,2))
       
          u(:,:,:,1) = u(:,:,:,2)
          v(:,:,:,1) = v(:,:,:,2)
          eta(:,:,:,1) = eta(:,:,:,2)
          Uek(:,:,1) = Uek(:,:,2)
          Vek(:,:,1) = Vek(:,:,2)
          u(:,:,:,2) = u(:,:,:,3)
          v(:,:,:,2) = v(:,:,:,3)
          eta(:,:,:,2) = eta(:,:,:,3)
          Uek(:,:,2) = Uek(:,:,3)
          Vek(:,:,2) = Vek(:,:,3) 

          ke1 = 0.
          ke2 = 0.
          do j = 1, ny
          do i = 1, nx
             ke1 = ke1+(u(i,j,1,3)**2+v(i,j,1,3)**2)/2.
             ke2 = ke2+(u(i,j,2,3)**2+v(i,j,2,3)**2)/2.
          enddo
          enddo
 
          write(300,*), time/86400., ke1/nx/ny, ke2/nx/ny
          call flush(300)

          if(nsteps.lt.start_movie) then
             if ( mod(its,iout).eq.0 ) then  ! output 
                 include 'subs/div_vort.f90'
                 include 'subs/dump_gnu1a.f90'
             endif
         endif
          if ( its .gt. start_movie ) then
             if ( mod(its,iout).eq.0 ) then  ! output 
                include 'subs/div_vort.f90'
!                include 'subs/tmp_complex.f90'
!                include 'subs/calc_q.f90'
                include 'subs/diags.f90'
                icount = icount + 1
                include 'subs/dump_bin.f90'
                print*, 'writing data No.', icount
             endif
          endif  !output

          if ( time.gt.0*86400. ) then
             if ( mod(its,ispechst).eq.0 ) then
                count_specs_1 = count_specs_1 + 1
                count_specs_2 = count_specs_2 + 1 
                count_specs_ek = count_specs_ek + 1
                count_specs_to = count_specs_to + 1
                count_specs_ag = count_specs_ag + 1
         
                ! 1-Dx
                if(mod(itlocal,ispechst)==0) then
                kex1_spec=0.0
                kex2_spec=0.0
                do j = 1,ny
                  ! u, v at iz=1
                datrx(:)=u(1:nx,j,1,2)
                include 'fftw_stuff/specx.f90'
                kex1_spec =  spectrumx
                datrx(:)=v(1:nx,j,1,2)
                include 'fftw_stuff/specx.f90'
                kex1_spec =  kex1_spec + spectrumx
                kex1_spec_tb(:,itsrow)=kex1_spec
                !u,v at iz=2
                datrx(:)=u(1:nx,j,2,2)
                include 'fftw_stuff/specx.f90'
                kex2_spec =  kex2_spec + spectrumx
                datrx(:)=v(1:nx,j,2,2)
                include 'fftw_stuff/specx.f90'
                kex2_spec =  kex2_spec + spectrumx
                kex2_spec_tb(:,itsrow)=kex2_spec
                !Uek 
                datrx(:) = Uek(1:nx,j,2)       
                include 'fftw_stuff/specx.f90'
                kex_ek_spec =  kex_ek_spec + spectrumx
         
                datrx(:) = Vek(1:nx,j,2)       
                include 'fftw_stuff/specx.f90'
                kex_ek_spec =  kex_ek_spec + spectrumx
                kex_ek_spec_tb(:,itsrow)=kex_ek_spec
                enddo
               
               !  ! 1-Dy
                key1_spec=0.0
                key2_spec=0.0
                do i = 1,nx
                  ! u, v at iz=1
                datry(:)=u(i,1:ny,1,2)
                include 'fftw_stuff/specy.f90'
                key1_spec =  key1_spec + spectrumy
                datry(:)=v(i,1:ny,1,2)
                include 'fftw_stuff/specy.f90'
                key1_spec =  key1_spec + spectrumy
                key1_spec_tb(:,itsrow)=key1_spec
                !u,v at iz=2
                datry(:)=u(i,1:ny,2,2)
                include 'fftw_stuff/specy.f90'
                key2_spec =  key2_spec + spectrumy
                datry(:)=v(i,1:ny,2,2)
                include 'fftw_stuff/specy.f90'
                key2_spec =  key2_spec + spectrumy
                key2_spec_tb(:,itsrow)=key2_spec
                !Uek,Vek
                datry(:) = Uek(i,1:ny,2)  
                include 'fftw_stuff/specy.f90'
                key_ek_spec =  key_ek_spec + spectrumy
                datry(:) = Vek(i,1:ny,2)  
                include 'fftw_stuff/specy.f90'
                key_ek_spec =  key_ek_spec + spectrumy
                key_ek_spec_tb(:,itsrow)=key_ek_spec
               enddo

               ! Write 1-D spectra
               WRITE(which2,'(I5)') 10000 + nspecfile
               if(itsrow==ntsrow) then
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
               
               end if !itsrow==ntsrow
               end if !every ispechst

               !  write(998,*) time/86400.,key1_spec
               !  write(997,*) time/86400.,kex2_spec
               !  write(996,*) time/86400.,key2_spec
               !  write(995,*) time/86400.,kex_ek_spec
               !  write(994,*) time/86400.,key_ek_spec
         
         
                         ! 2D
                datr(:,:) = u(1:nx,1:ny,1,2)       
                include 'fftw_stuff/spec1.f90'
                ke1_spec =  ke1_spec + spectrum

                datr(:,:) = v(1:nx,1:ny,1,2)       
                include 'fftw_stuff/spec1.f90'
                ke1_spec =  ke1_spec + spectrum

                datr(:,:) = u(1:nx,1:ny,2,2)
                include 'fftw_stuff/spec1.f90'
                ke2_spec = ke2_spec + spectrum

                datr(:,:) = v(1:nx,1:ny,2,2)
                include 'fftw_stuff/spec1.f90'
                ke2_spec = ke2_spec + spectrum
                          
                datr(:,:) = forcing_total(1:nx,1:ny)
                include 'fftw_stuff/spec1.f90'
                for_to_spec = for_to_spec + spectrum

                datr(:,:) = forcing_ag(1:nx,1:ny)
                include 'fftw_stuff/spec1.f90'
                for_ag_spec = for_ag_spec + spectrum

                datr(:,:) = Uek(1:nx,1:ny,2)       
                include 'fftw_stuff/spec1.f90'
                ke_ek_spec =  ke_ek_spec + spectrum

                datr(:,:) = Vek(1:nx,1:ny,2)       
                include 'fftw_stuff/spec1.f90'
                ke_ek_spec =  ke_ek_spec + spectrum
             endif
          endif

          if ( time.gt.0*86400. ) then
            if ( mod(its,ispechst).eq.0 ) then
            endif
         endif


!          if ( mod(its,i_diags).eq.0 ) then
!             include 'subs/diags.f90'
!          endif

       enddo ! its

       open(99, file = 'ke1_spec')
       do i = 1,ny
          write(99,*) i, ke1_spec(i)/nx/nx/ny/ny/count_specs_1
       enddo
       close(99)

       open(98, file = 'ke_ek_spec')
       do i = 1,ny
          write(98,*) i, ke_ek_spec(i)/nx/nx/ny/ny/count_specs_ek
       enddo
       close(98)
       
       open(97, file = 'ke2_spec')
       do i = 1,ny
          write(97,*) i, ke2_spec(i)/nx/nx/ny/ny/count_specs_2
       enddo
       close(97)




 !      open(96, file = 'for_to_spec')
 !      do i = 1,ny
 !         write(96,*) i, for_to_spec(i)/nx/nx/ny/ny/count_specs_to
 !      enddo
 !      close(96)
   
 !      open(95, file = 'for_ag_spec')
 !      do i = 1,ny
 !         write(95,*) i, for_ag_spec(i)/nx/nx/ny/ny/count_specs_ag
 !      enddo
 !      close(95)

       include 'fftw_stuff/fft_destroy.f90'

       open(0,file='restart')
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
       write(0,*) icount,time
       close(0)


      end program main
