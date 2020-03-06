!
! top layer
!
       uu = u1_filtered(:,:)
       vv = v1_filtered(:,:)

       array = uu
       include 'subs/bndy.f90'
       uu = array

       array = vv
       include 'subs/bndy.f90'
       vv = array

       datr(:,:) = uu(1:nx,1:ny)
       include 'fftw_stuff/invLaplacian.f90'
       invLap_u(1:nx,1:ny) = datr(:,:)

       datr(:,:) = vv(1:nx,1:ny)
       include 'fftw_stuff/invLaplacian.f90'
       invLap_v(1:nx,1:ny) = datr(:,:)

       u1diss(:,:) =  r_invLap*invLap_u(:,:) 
       v1diss(:,:) =  r_invLap*invLap_v(:,:) 

!
! bottom layer
!

       uu = u2_filtered(:,:)
       vv = v2_filtered(:,:)

       array = uu
       include 'subs/bndy.f90'
       uu = array

       array = vv
       include 'subs/bndy.f90'
       vv = array

       datr(:,:) = uu(1:nx,1:ny)
       include 'fftw_stuff/invLaplacian.f90'
       invLap_u(1:nx,1:ny) = datr(:,:)

       datr(:,:) = vv(1:nx,1:ny)
       include 'fftw_stuff/invLaplacian.f90'
       invLap_v(1:nx,1:ny) = datr(:,:)

       u2diss(:,:) =  r_invLap*invLap_u(:,:)!& 
!      &            - r_drag*uu(:,:)
       v2diss(:,:) =  r_invLap*invLap_v(:,:)!&
!      &            - r_drag*vv(:,:)

