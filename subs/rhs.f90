
       tmp(1) = minval(thickness)/H(k)
       if(tmp(1).le.0.02) then
       print*, 'thickness too small'
       stop
       endif

       array = uu
       include 'subs/bndy.f90'
       uu = array

       array = vv
       include 'subs/bndy.f90'
       vv = array

       array = pressure
       include 'subs/bndy.f90'
       pressure = array

       array = thickness
       include 'subs/bndy.f90'
       thickness = array

       array = uu_old
       include 'subs/bndy.f90'
       uu_old = array

       array = vv_old
       include 'subs/bndy.f90'
       vv_old = array

       datr(:,:) = uu_old(1:nx,1:ny)
       include 'fftw_stuff/invLaplacian.f90'
       invLap_u(1:nx,1:ny) = datr(:,:)

       datr(:,:) = vv_old(1:nx,1:ny)
       include 'fftw_stuff/invLaplacian.f90'
       invLap_v(1:nx,1:ny) = datr(:,:)

!
!   note: div is just for diagnostics
!
       do j = 1,ny
       do i = 1,nx

       div(i,j) = (uu(i+1,j)-uu(i,j))/dx   &
       &        + (vv(i,j+1)-vv(i,j))/dy 

       B(i,j) = 0.25*(uu(i,j)**2+uu(i+1,j)**2+vv(i,j)**2+vv(i,j+1)**2)
       B(i,j) = B(i,j) + pressure(i,j)

       zeta(i,j) = (vv(i,j)-vv(i-1,j))/dx   &
       &         - (uu(i,j)-uu(i,j-1))/dy 

       grad2u(i,j) = (uu_old(i+1,j)+uu_old(i-1,j)-2.*uu_old(i,j))/dx/dx   &
       &           + (uu_old(i,j+1)+uu_old(i,j-1)-2.*uu_old(i,j))/dy/dy   
          
       grad2v(i,j) = (vv_old(i+1,j)+vv_old(i-1,j)-2.*vv_old(i,j))/dx/dx   &
       &           + (vv_old(i,j+1)+vv_old(i,j-1)-2.*vv_old(i,j))/dy/dy   

       uh(i,j) = 0.5*(thickness(i,j)+thickness(i-1,j))*uu(i,j)
       vh(i,j) = 0.5*(thickness(i,j)+thickness(i,j-1))*vv(i,j)
     
       enddo
       enddo


       array = uh
       include 'subs/bndy.f90'
       uh = array

       array = vh
       include 'subs/bndy.f90'
       vh = array

       array = B
       include 'subs/bndy.f90'
       B = array

       array = zeta
       include 'subs/bndy.f90'
       zeta = array

       array = grad2u
       include 'subs/bndy.f90'
       grad2u = array

       array = grad2v
       include 'subs/bndy.f90'
       grad2v = array


       do j = 1, ny
          jp = j+1
          jm = j-1
       do i = 1, nx
          ip = i+1
          im = i-1
       grad4u(i,j) = (grad2u(ip,j)+grad2u(im,j)-2.*grad2u(i,j))/dx/dx   &
       &           + (grad2u(i,jp)+grad2u(i,jm)-2.*grad2u(i,j))/dy/dy   
       grad4v(i,j) = (grad2v(ip,j)+grad2v(im,j)-2.*grad2v(i,j))/dx/dx   &
       &           + (grad2v(i,jp)+grad2v(i,jm)-2.*grad2v(i,j))/dy/dy   
       enddo
       enddo

       do j = 1, ny
          jp = j+1
          jm = j-1
       do i = 1, nx
          ip = i+1
          im = i-1

       rhs_u(i,j,k) = -(B(i,j)-B(im,j))/dx   &
       &            + 0.25*(f(j)+zeta(i,j))*(vv(i,j)+vv(im,j))   &
       &            + 0.25*(f(jp)+zeta(i,jp))*(vv(i,jp)+vv(im,jp))   &
       &            - Ah*grad4u(i,j)   &
       &            + r_invLap*invLap_u(i,j)*sig_rinvL(k)  &
       &            - r_drag*uu_old(i,j)*bot(k)
 
       rhs_v(i,j,k) = -(B(i,j)-B(i,jm))/dy   &
       &            - 0.25*(f(j)+zeta(i,j))*(uu(i,j)+uu(i,jm))   &
       &            - 0.25*(f(jp)+zeta(ip,j))*(uu(ip,j)+uu(ip,jm))   &
       &            - Ah*grad4v(i,j)   &
       &            + r_invLap*invLap_v(i,j)*sig_rinvL(k)   &
       &            - r_drag*vv_old(i,j)*bot(k)

       rhs_eta(i,j,k) = -(uh(ip,j)-uh(i,j))/dx   &
       &                -(vh(i,jp)-vh(i,j))/dy   &
       &                - top(k)*(Uek(ip,j,2)-Uek(i,j,2))/dx   & 
       &                - top(k)*(Vek(i,jp,2)-Uek(i,j,2))/dy 
 
       enddo
       enddo

