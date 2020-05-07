
            ! 2D !ufft,vfft,etafft,ftotalfft,f_agfft
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
            f_agfft=datc

            datr(:,:) = Uek(1:nx,1:ny,2)       
            include 'fftw_stuff/spec1.f90'
            ke_ek_spec =  ke_ek_spec + spectrum

            datr(:,:) = Vek(1:nx,1:ny,2)       
            include 'fftw_stuff/spec1.f90'
            ke_ek_spec =  ke_ek_spec + spectrum

            datr(:,:) = u_ag(1:nx,1:ny,1)
            include 'fftw_stuff/spec1.f90'
            u_agfft(:,:,1)=datc        

            datr(:,:) = u_ag(1:nx,1:ny,2)
            include 'fftw_stuff/spec1.f90'
            u_agfft(:,:,2)=datc

            datr(:,:) = v_ag(1:nx,1:ny,1)
            include 'fftw_stuff/spec1.f90'
            v_agfft(:,:,1)=datc       

            datr(:,:) = v_ag(1:nx,1:ny,2)
            include 'fftw_stuff/spec1.f90'
            v_agfft(:,:,2)=datc
