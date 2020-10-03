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

               datr(:,:) = u_ag(1:nx,1:ny,1)
               include 'fftw_stuff/spec1.f90'
               ke1_ag_spec = ke1_ag_spec + spectrum
               u_agfft(:,:,1)=datc        
   
               datr(:,:) = u_ag(1:nx,1:ny,2)
               include 'fftw_stuff/spec1.f90'
               ke2_ag_spec = ke2_ag_spec + spectrum
               u_agfft(:,:,2)=datc
   
               datr(:,:) = v_ag(1:nx,1:ny,1)
               include 'fftw_stuff/spec1.f90'
               ke1_ag_spec = ke1_ag_spec + spectrum
               v_agfft(:,:,1)=datc       
   
               datr(:,:) = v_ag(1:nx,1:ny,2)
               include 'fftw_stuff/spec1.f90'
               ke2_ag_spec = ke2_ag_spec + spectrum
               v_agfft(:,:,2)=datc
   
               datr(:,:) = u_qg(1:nx,1:ny,1)
               include 'fftw_stuff/spec1.f90'
               ke1_qg_spec = ke1_qg_spec + spectrum
               u_qgfft(:,:,1)=datc        
   
               datr(:,:) = u_qg(1:nx,1:ny,2)
               include 'fftw_stuff/spec1.f90'
               ke2_qg_spec = ke2_qg_spec + spectrum
               u_qgfft(:,:,2)=datc
   
               datr(:,:) = v_qg(1:nx,1:ny,1)
               include 'fftw_stuff/spec1.f90'
               ke1_qg_spec = ke1_qg_spec + spectrum
               v_qgfft(:,:,1)=datc       
   
               datr(:,:) = v_qg(1:nx,1:ny,2)
               include 'fftw_stuff/spec1.f90'
               ke2_qg_spec = ke2_qg_spec + spectrum
               v_qgfft(:,:,2)=datc
               
               datr(:,:) = Psurf_csum(1:nx,1:ny)
               include 'fftw_stuff/spec1.f90'
               etafft(:,:,1)=datc

               datr(:,:) = eta(1:nx,1:ny,2,3)
               include 'fftw_stuff/spec1.f90'
               etafft(:,:,2)=datc

               datr(:,:) = eta_qg(1:nx,1:ny)
               include 'fftw_stuff/spec1.f90'
               eta_qgfft(:,:)=datc

               u_agfft_bc=u_agfft(:,:,2)-u_agfft(:,:,1)
               v_agfft_bc=v_agfft(:,:,2)-v_agfft(:,:,1)
               ! write(*,*) '2D FFT/spec done, its,time,iftcount',its,time/86400,iftcount