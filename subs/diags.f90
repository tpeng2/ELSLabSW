!
!      calc pv
!      called after time steps switched, so 2 and 3 are both at "n+1"
!      time level
!
       do k = 1,nz
          uu(:,:) = u(:,:,k,3)  
          vv(:,:) = v(:,:,k,3)
          if (k.eq.1) then
             thickness(:,:) =  -eta(:,:,k+1,3)/H(1)
          elseif (k.eq.nz) then
             thickness(:,:) =  eta(:,:,k,3)/H(k)
          else
             thickness(:,:) =  (eta(:,:,k,3)-eta(:,:,k+1,3))/H(k)
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
            enddo
          enddo

          array = zeta
          include 'subs/bndy.f90'
          zeta = array
 
          do j = 1, ny
            do i = 1, nx
               q(i,j,k) = zeta(i,j) -0.25*f(j)*(thickness(i,j)+    &
                        &                       thickness(i-1,j)+  &
                        &                       thickness(i,j-1)+  &
                        &                       thickness(i-1,j-1) )
            enddo
          enddo

          array(:,:) = q(:,:,k)
          include 'subs/bndy.f90'
          q(:,:,k) = array(:,:)
       enddo ! k-loop

       if (nz.eq.2) then
          qmode(:,:,1) = (H(1)*q(:,:,1) + H(2)*q(:,:,2))/Htot
          qmode(:,:,2) = q(:,:,2) - q(:,:,1)
          Fmode(1) = 0.
          Fmode(2) = (f0/c_bc)**2
       else
          print*,  'error, nx /= 2'
       endif

       !  invert to get psimode
 
       do imode = 1,2
          datr(:,:) = qmode(1:nx,1:ny,k)
          include 'fftw_stuff/invertQ.f90'
          psimode(1:nx,1:ny,k) = datr(:,:)

          array(:,:) = psimode(:,:,k)
          include 'subs/bndy.f90'
          psimode(:,:,k) = array(:,:)
       enddo

       psi(:,:,1) = psimode(:,:,1) - H(2)*psimode(:,:,2)/Htot
       psi(:,:,2) = psimode(:,:,1) + H(1)*psimode(:,:,2)/Htot

       do k = 1,nz ! for k=1,2 two layers
          array(:,:) = psi(:,:,k)
          include 'subs/bndy.f90'
          psi(:,:,k) = array(:,:)
          do j = 1,ny
            do i = 1,nx
               v_qg(i,j,k) =  (psi(i+1,j,k)-psi(i,j,k))/dx
               u_qg(i,j,k) = -(psi(i,j+1,k)-psi(i,j,k))/dy
            enddo
          enddo

          array(:,:) = u_qg(:,:,k)
          include 'subs/bndy.f90'
          u_qg(:,:,k) = array(:,:)

          array(:,:) = v_qg(:,:,k)
          include 'subs/bndy.f90'
          v_qg(:,:,k) = array(:,:)
      
          !  bndy on u and v is not necessary
          !  had been done in p_correction for each time level
          u_ag(:,:,k) = u(:,:,k,2)-u_qg(:,:,k)  
          v_ag(:,:,k) = v(:,:,k,2)-v_qg(:,:,k)

          array(:,:) = u_ag(:,:,k)
          include 'subs/bndy.f90'
          u_ag(:,:,k) = array(:,:)

          array(:,:) = v_ag(:,:,k)
          include 'subs/bndy.f90'
          v_ag(:,:,k) = array(:,:)
       end do

       
       !essentially psi(,,2)-psi(,,1) is psimode(,,2)
       do j = 1,ny
         do i = 1,nx    
            eta_qg(i,j) = (f0/gprime(2))*0.25*(psi(i,j,2)-psi(i,j,1)   &
                        &               + psi(i+1,j,2)-psi(i+1,j,1)   &
                        &               + psi(i,j+1,2)-psi(i,j+1,1)   &
                        &               + psi(i+1,j+1,2)-psi(i+1,j+1,1))
         enddo
       enddo  

       !  iterate to get eta_G at eta points
       array1 = eta_qg
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
         array1 = eta_qg - array

       enddo
       eta_qg = array2  ! this is eta for the G modes at eta points


       array(:,:) = eta_qg(:,:)
       include 'subs/bndy.f90'
       eta_qg(:,:) = array(:,:)
  
       array(:,:) = eta(:,:,2,3)
       include 'subs/bndy.f90'
       eta(:,:,2,3) = array(:,:) 
   
       eta_ag(:,:) = eta(:,:,2,3)-eta_qg(:,:)
      
       array(:,:) = eta_ag(:,:)
       include 'subs/bndy.f90'
       eta_ag(:,:) = array(:,:)

       do k=1,nz
         uu=u_qg(:,:,k)
         vv=v_qg(:,:,k)
         ! zeta_G
         do j = 1,ny
            do i = 1,nx
               zeta_G(i,j,k) =  (vv(i,j)-vv(i-1,j))/dx    &
               &         -  (uu(i,j)-uu(i,j-1))/dy 
            enddo ! i-loop
         enddo ! j-koop
         array = zeta_G(:,:,k)
         include 'subs/bndy.f90'
         zeta_G(:,:,k)= array
       enddo ! k-loop



!     ke1 = 0.
!     ke2 = 0.
!     do j = 1, ny
!     do i = 1, nx
!        ke1 = ke1 + (u(i,j,1,2)**2 + v(i,j,1,2)**2)/2.
!        ke2 = ke2 + (u(i,j,2,2)**2 + v(i,j,2,2)**2)/2.
!     enddo
!     enddo

!     ke1_qg = 0.
!     ke2_qg = 0.
!     do j = 1, ny
!     do i = 1, nx
!        ke1_qg = ke1_qg + ((psi(i+1,j,1)-psi(i,j,1))/dx)**2   &
!        &      + ((psi(i,j+1,1)-psi(i,j,1))/dy)**2 
!        ke2_qg = ke2_qg + ((psi(i+1,j,2)-psi(i,j,2))/dx)**2   &
!        &      + ((psi(i,j+1,2)-psi(i,j,2))/dy)**2 
!     enddo
!     enddo
!     ke1_qg = ke1_qg/2.
!     ke2_qg = ke2_qg/2.

!     pe_qg = 0.
!     pe = 0.
!     do j = 1,ny
!     do i = 1,nx
!        pe_qg = pe_qg + (f0/c_bc)**2 * (psi(i,j,2)-psi(i,j,1))**2 
!        pe = pe + gprime(2)*eta(i,j,2,2)**2
!     enddo
!     enddo
!     pe_qg = pe_qg/2.
!     pe = pe/2.
!     pe = pe*Htot/H(1)/H(2)
!     etot = (H(1)*ke1 + H(2)*ke2)/Htot + pe
!     etot_qg = (H(1)*ke1_qg + H(1)*ke2_qg)/Htot + pe_qg

!      write(300,*), time/86400., ke1/nx/ny, ke2/nx/ny, pe/nx/ny   &
!          &                   , etot/nx/ny
!      write(301,*), time/86400., ke1_qg/nx/ny, ke2_qg/nx/ny   &
!          &                   , pe_qg/nx/ny, etot_qg/nx/ny 
!      write(302,*), time/86400., eta(1,ny/2,2,1), eta(1,ny/4,2,1)   &
!          &                   , eta(1,1,2,1)
!      call flush(300)
!      call flush(301)
!      call flush(302)
!
!
