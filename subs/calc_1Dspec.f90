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