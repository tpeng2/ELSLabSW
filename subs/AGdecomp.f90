!
!      if called after time steps switched, so 2 and 3 are both at "n+1"
!      time level
!
       do k = 1,nz

!      uu(:,:) = u(:,:,k,2)  
!      vv(:,:) = v(:,:,k,2)
       if(k.eq.1) then
       uu = u1_snap - u1_filtered
       vv = v1_snap - v1_filtered
       thickness =  -(eta_snap - eta_filtered)/H(1)
       elseif(k.eq.nz) then
       thickness(:,:) =  eta(:,:,nz,2)/H(nz)
       uu = u2_snap - u2_filtered
       vv = v2_snap - v2_filtered
       thickness =  (eta_snap - eta_filtered)/H(nz)
       else
           print*, 'whoops'
           stop
       thickness(:,:) =  (eta(:,:,k,2)-eta(:,:,k+1,2))/H(k)
       endif

       array = uu
       include 'subs/bndy.f90'
       uu = array

       array = vv
       include 'subs/bndy.f90'
       vv = array

       array = thickness
       include 'subs/bndy.f90'
       thickness = array

       do j = 1,ny
       do i = 1,nx

       zeta(i,j) =  (vv(i,j)-vv(i-1,j))/dx    &
       &         -  (uu(i,j)-uu(i,j-1))/dy 
       q(i,j,k) = zeta(i,j) -0.25*f0*(thickness(i,j)+    &
                &                     thickness(i-1,j)+  &
                &                     thickness(i,j-1)+  &
                &                     thickness(i-1,j-1) )

       enddo
       enddo

       array(:,:) = q(:,:,k)
       include 'subs/bndy.f90'
       q(:,:,k) = array(:,:)

       enddo ! k loop

!      modal decomposition is only set up for nz = 2
       if(nz.eq.2) then
       qmode(:,:,1) = (H(1)*q(:,:,1) + H(2)*q(:,:,2))/Htot
       qmode(:,:,2) =  q(:,:,2) - q(:,:,1)
       Fmode(1) = 0.
       Fmode(2) =  (f0/c_bc)**2

       !  invert to get psimode

       do k = 1,nz

       datr(:,:) = qmode(1:nx,1:ny,k)
       include 'fftw_stuff/invertQ.f90'
       psimode(1:nx,1:ny,k) = datr(:,:)

       array(:,:) = psimode(:,:,k)
       include 'subs/bndy.f90'
       psimode(:,:,k) = array(:,:)

       enddo

       psi(:,:,1) = psimode(:,:,1) - H(2)*psimode(:,:,2)/Htot
       psi(:,:,2) = psimode(:,:,1) + H(1)*psimode(:,:,2)/Htot

       endif ! if nz = 2


       do j = 1,ny
       do i = 1,nx
       uu_G(i,j) =  (psimode(i,j,2)-psimode(i,j+1,2))/dy
       vv_G(i,j) =  (psimode(i+1,j,2)-psimode(i,j,2))/dx
       enddo
       enddo
       uu_A(:,:) =  u(:,:,2,2)-u(:,:,1,2) - uu_G(:,:)
       vv_A(:,:) =  v(:,:,2,2)-v(:,:,1,2) - vv_G(:,:)

       ! below gives eta at zeta points, need to interpolate
       ! and iterate to get at eta points
       do j = 1,ny
       do i = 1,nx
       eta_G(i,j) =  psimode(i,j,2)*f0/gprime(2)
       enddo
       enddo

       array = uu_G(:,:)
       include 'subs/bndy.f90'
       uu_G(:,:) = array

       array = vv_G(:,:)
       include 'subs/bndy.f90'
       vv_G(:,:) = array

       array = eta_G(:,:)
       include 'subs/bndy.f90'
       eta_G(:,:) = array

       !  iterate to get eta_G at eta points

       array1 = eta_G
       array2 = 0.
       do ii = 1,10

       do j = 1,ny
       do i = 1,nx
       array(i,j) = 0.25*(array1(i,j)+array1(i+1,j)+ &
           &              array1(i,j+1)+array1(i+1,j+1))
       enddo
       enddo
       include 'subs/bndy.f90'
       array2 = array2 + array ! estimate of eta_G at eta points

       do j = 1,ny ! calc eta_G at zeta points as implied by array2
       do i = 1,nx
       array(i,j) = 0.25*(array2(i,j)+array2(i,j-1)+ &
           &              array2(i-1,j)+array2(i-1,j-1))
       enddo
       enddo
       include 'subs/bndy.f90'
       array1 = eta_G - array

       enddo
       eta_G = array2  ! this is eta for the G modes at eta points
       eta_A =  eta_snap - eta_filtered - eta_G
