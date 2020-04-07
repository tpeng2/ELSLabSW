      module data_initial
      integer nx, ny, nz, nnx, nny
      integer i_diags
      double precision pi, twopi, Lx, Ly, dx, dy
      real f0, beta, r_drag, Ah, r_invLap, rf
      real tau0, tau1, hek
      integer nsteps,fileperday,start_movie
      real ndays,totaltime,dt
      real restart_from
      integer subsmprto,itape,ispechst,iout,itlocal,itsrow,ntsrow,nspecfile
      integer ftsubsmprto
      logical restart, use_ramp, ifsteady
      logical calc1Dspec,save_movie,save2dfft
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
      real u_out(1:(nx/subsmprto),1:(ny/subsmprto),nz)
      real v_out(1:(nx/subsmprto),1:(ny/subsmprto),nz)
      real eta_out(1:(nx/subsmprto),1:(ny/subsmprto),nz)
      real div_ek_out(1:(nx/subsmprto),1:(ny/subsmprto))
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
      double complex,dimension(nx/2+1,ny,nz) :: ufft,vfft
      double complex,dimension(nx/2+1,ny) :: etafft,ftotalfft,fagfft

      integer i, j, k, ii, jj, kk, ip, im, jp, jm, kp, km, it, its, ntimes, inkrow
      !subsampling arrays
      integer,allocatable:: isubx(:),isuby(:),iftsubkl(:,:)
      integer rdsubk,rdsubl !temporary variables for reading (k,l) pair 
      !subsampling size
      integer szsubx,szsuby,szftrdrow,szftrdcol
      !I/O info
      integer icount,iftcount, count_specs_1, count_specs_2, count_specs_ek
      integer count_specs_to, count_specs_ag

      character(88) scrap, which,which2,which3, spectbszx,spectbszy,specform, pathspects
      character(88) string1, string2, string3, string4, string5
      character(88) string1i, string2i, string3i, string4i, string5i
      character(88) string6, string7, string8, string9, string10
      character(88) string11, string12, string13, string14, string15
      character(88) string16, string17, string18, string19, string20

      include 'fftw_stuff/fft_params.f90'
      include 'fftw_stuff/fft_init.f90'
      ! === Allocate variables
      szsubx=ceiling(nx/(subsmprto+1e-15))
      szsuby=ceiling(ny/(subsmprto+1e-15))
      ! szftrdrow=ceiling((nx/2+1)/(ftsubsmprto+1e-15))
      ! szftrdcol=ceiling(ny/(ftsubsmprto+1e-15))
      ! === Define subsampling range in spatial space
      isubx=(/(i, i=1,nx, subsmprto)/)
      isuby=(/(i, i=1,ny, subsmprto)/)
      ! === FFT subsmpling, fed with (k,l) pair
      include 'subs/read_kxky_subsmp.f90'

      ! === 
      !set icount and time based on if restart is true or false !moved to initialize
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
      if(restart==.false.)   nspecfile=0
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

         include 'subs/IOheader.f90' 
         
      !==============================================================      
      !
      !      subsequent time steps
      !
      !==============================================================      

      do its = 2, nsteps

         ! determine local time step !tstime
         itlocal=mod(its,itape)       ! for 1-D time series, relative its for itape
         itsrow=int(itlocal/ispechst) ! xxth row in a 1-D time series file
         ! forced to set 0-mod to the last row
         if (itsrow==0.and.its.gt.ispechst) itsrow=ntsrow

         if(itlocal==0) then
            itlocal = itape       ! for 0-mod
            nspecfile=nspecfile+1 ! increase the file index by 1
         end if

         if(mod(itlocal,ispechst)==0) then
            tstime(itsrow)=time/86400. !record each time when print out 
            ! write(*,*) 'its,itsrow,time[day]',its,itsrow,tstime(itsrow)
         end if

         !        =================
         if ( mod(its,iout*10).eq.1 ) write(*,*) 'current time step, iout', its, iout
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

         if(nsteps.lt.start_movie.and.save_movie) then
            if ( mod(its,iout).eq.0 ) then  ! output 
               include 'subs/div_vort.f90'
               include 'subs/dump_gnu1a.f90'
            end if
         end if

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
               if(mod(itlocal,ispechst)==0.and.calc1Dspec) then
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

                  if(itsrow==ntsrow) then
                     include 'subs/IOoutput.f90'
                  endif !itsrow==ntsrow
               endif !calc1Dspectrum and every ispechst

               ! 2D !ufft,vfft,etafft,ftotalfft,fagfft
               datr(:,:) = u(1:nx,1:ny,1,2)       
               include 'fftw_stuff/spec1.f90'
               ke1_spec =  ke1_spec + spectrum
               ufft(:,:,1)=datc

               datr(:,:) = v(1:nx,1:ny,1,2)       
               include 'fftw_stuff/spec1.f90'
               ke1_spec =  ke1_spec + spectrum
               vfft(:,:,1)=datc

               datr(:,:) = u(1:nx,1:ny,2,2)
               include 'fftw_stuff/spec1.f90'
               ke2_spec = ke2_spec + spectrum
               ufft(:,:,2)=datc

               datr(:,:) = v(1:nx,1:ny,2,2)
               include 'fftw_stuff/spec1.f90'
               ke2_spec = ke2_spec + spectrum
               vfft(:,:,2)=datc

               datr(:,:) = forcing_total(1:nx,1:ny)
               include 'fftw_stuff/spec1.f90'
               for_to_spec = for_to_spec + spectrum
               ftotalfft=datc

               datr(:,:) = forcing_ag(1:nx,1:ny)
               include 'fftw_stuff/spec1.f90'
               for_ag_spec = for_ag_spec + spectrum
               fagfft=datc

               datr(:,:) = Uek(1:nx,1:ny,2)       
               include 'fftw_stuff/spec1.f90'
               ke_ek_spec =  ke_ek_spec + spectrum

               datr(:,:) = Vek(1:nx,1:ny,2)       
               include 'fftw_stuff/spec1.f90'
               ke_ek_spec =  ke_ek_spec + spectrum

               iftcount=iftcount+1
               write(*,*) '2D FFT/spec done, its,time,iftcount',its,time/86400,iftcount
               
               ! Write 2-D FFT fields
               if(save2dfft) then
                  include 'subs/dump_bin_spec2d.f90'
               endif !itsrow==ntsrow
            endif
         endif


         ! ========= From the original code
         !          if ( mod(its,i_diags).eq.0 ) then
         !             include 'subs/diags.f90'
         !          endif

      enddo ! its

      ! ========= I/O after time loop =========
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

      open(96, file = 'for_to_spec')
      do i = 1,ny
         write(96,*) i, for_to_spec(i)/nx/nx/ny/ny/count_specs_to
      enddo
      close(96)

      open(95, file = 'for_ag_spec')
      do i = 1,ny
         write(95,*) i, for_ag_spec(i)/nx/nx/ny/ny/count_specs_ag
      enddo
      close(95)

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
      write(0,*) icount,time,nspecfile,iftcount
      close(0)

      end program main
