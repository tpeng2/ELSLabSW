
!
!    bndy condition: periodic in x,y
!
       
       array = u(:,:,1,3)
       include 'subs/bndy.f90'
       u(:,:,1,3) = array
       array = v(:,:,1,3)
       include 'subs/bndy.f90'
       v(:,:,1,3) = array
       array = u(:,:,2,3)
       include 'subs/bndy.f90'
       u(:,:,2,3) = array
       array = v(:,:,2,3)
       include 'subs/bndy.f90'
       v(:,:,2,3) = array
       array = Uek(:,:,3)
       include 'subs/bndy.f90'
       Uek(:,:,3) = array
       array = Vek(:,:,3)
       include 'subs/bndy.f90'
       Vek(:,:,3) = array


!
!   note: div is just for diagnostics
!
       do j = 1, ny
          jp = j+1
          jm = j-1
       do i = 1, nx
          ip = i+1
          im = i-1

       div1(i,j) = (u(ip,j,1,3)-u(i,j,1,3))/dx + (v(i,jp,1,3)-v(i,j,1,3))/dy

       div2(i,j) = (u(ip,j,2,3)-u(i,j,2,3))/dx + (v(i,jp,2,3)-v(i,j,2,3))/dy 
    
       div_ek(i,j) = (Uek(ip,j,3)-Uek(i,j,3))/dx + (Vek(i,jp,3)-Vek(i,j,3))/dy
   
       zeta1(i,j) = (v(i,j,1,3)-v(im,j,1,3))/dx - (u(i,j,1,3)-u(i,jm,1,3))/dy 

       zeta2(i,j) = (v(i,j,2,3)-v(im,j,2,3))/dx - (u(i,j,2,3)-u(i,jm,2,3))/dy

       zeta_ek(i,j) = (Vek(i,j,3)-Vek(im,j,3))/dx - (Uek(i,j,3)-Uek(i,jm,3))/dy

       enddo
       enddo

       array = div1(:,:)
       include 'subs/bndy.f90'
       div1(:,:) = array
       array = div2(:,:)
       include 'subs/bndy.f90'
       div2(:,:) = array
       array = div_ek(:,:)
       include 'subs/bndy.f90'
       div_ek(:,:) = array
       array = zeta1(:,:)
       include 'subs/bndy.f90'
       zeta1(:,:) = array
       array = zeta2(:,:)
       include 'subs/bndy.f90'
       zeta2(:,:) = array
       array = zeta_ek(:,:)
       include 'subs/bndy.f90'
       zeta_ek(:,:) = array


