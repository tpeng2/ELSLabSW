!
!      calc pv
!
       
       do k = 1,nz

          uu(:,:) = u(:,:,k,2)  
          vv(:,:) = v(:,:,k,2)
          if (k.eq.1) then
             thickness(:,:) =  H(1)-eta(:,:,k+1,2)
          elseif (k.eq.nz) then
             thickness(:,:) =  H(2)+eta(:,:,k,2)
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
             q(i,j,k) = (f0+zeta(i,j))/(0.25*(thickness(i,j)+    &
                      &                       thickness(i-1,j)+  &
                      &                       thickness(i,j-1)+  &
                      &                       thickness(i-1,j-1)))

          enddo
          enddo

          array(:,:) = q(:,:,k)
          include 'subs/bndy.f90'
          q(:,:,k) = array(:,:)
  
       enddo ! k loop


