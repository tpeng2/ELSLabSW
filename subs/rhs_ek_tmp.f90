!
!    bndy condition: periodic in x,y
!
       array = uu
       include 'subs/bndy.f90'
       uu = array

       array = vv
       include 'subs/bndy.f90'
       vv = array

       array = uu1
       include 'subs/bndy.f90'
       uu1 = array

       array = vv1
       include 'subs/bndy.f90'
       vv1 = array

       array = uu_old
       include 'subs/bndy.f90'
       uu_old = array

       array = vv_old
       include 'subs/bndy.f90'
       vv_old = array

       print*, 'in rhs_ek_tmp'

       do j = 1,ny
       do i = 1,nx

       B(i,j) =  0.5*(uu(i,j)*uu1(i,j) + uu(i+1,j)*uu1(i+1,j) &
           &         +vv(i,j)*vv1(i,j) + vv(i,j+1)*vv1(i,j+1))
         
       B_nl(i,j) = 0.25*(uu(i,j)*uu(i,j)+uu(i+1,j)*uu(i+1,j) &
          &            + vv(i,j)*vv(i,j)+vv(i,j+1)*vv(i,j+1))

       zeta(i,j) =  (vv(i,j)-vv(i-1,j))/dx    &
       &         -  (uu(i,j)-uu(i,j-1))/dy 
       zeta1(i,j) =  (vv1(i,j)-vv1(i-1,j))/dx    &
       &          -  (uu1(i,j)-uu1(i,j-1))/dy 

       div(i,j) = (uu1(i+1,j)-uu1(i,j))/dx &
              & + (vv1(i,j+1)-vv1(i,j))/dy ! dns

       grad2u(i,j) = (uu_old(i+1,j)+uu_old(i-1,j)-2.*uu_old(i,j))/dx/dx   &
       &           + (uu_old(i,j+1)+uu_old(i,j-1)-2.*uu_old(i,j))/dy/dy   
          
       grad2v(i,j) = (vv_old(i+1,j)+vv_old(i-1,j)-2.*vv_old(i,j))/dx/dx   &
       &           + (vv_old(i,j+1)+vv_old(i,j-1)-2.*vv_old(i,j))/dy/dy   

       enddo
       enddo
       print*, 'hello'

       array = div !dns
       include 'subs/bndy.f90'
       div = array

       array = B
       include 'subs/bndy.f90'
       B = array
       
       array = B_nl
       include 'subs/bndy.f90'
       B_nl = array

       array = zeta
       include 'subs/bndy.f90'
       zeta = array

       array = zeta1
       include 'subs/bndy.f90'
       zeta1 = array

       array = grad2u
       include 'subs/bndy.f90'
       grad2u = array

       array = grad2v
       include 'subs/bndy.f90'
       grad2v = array


       print*, 'hello1'
       do j = 1,ny
          jp = j+1
          jm = j-1
       do i = 1,nx
          ip = i + 1
          im = i - 1
       grad4u(i,j) = (grad2u(ip,j)+grad2u(im,j)-2.*grad2u(i,j))/dx/dx   &
       &           + (grad2u(i,jp)+grad2u(i,jm)-2.*grad2u(i,j))/dy/dy   
       grad4v(i,j) = (grad2v(ip,j)+grad2v(im,j)-2.*grad2v(i,j))/dx/dx   &
       &           + (grad2v(i,jp)+grad2v(i,jm)-2.*grad2v(i,j))/dy/dy   
       enddo
       enddo
       print*, 'hello2'

       do j = 1,ny
       print*, j
          jp = j+1
          jm = j-1
       do i = 1,nx
          ip = i+1
          im = i-1

       rhs_Uek(i,j) =  -(B(i,j)-B(im,j))/dx  &
       &               -(1/hek)*(B_nl(i,j)-B_nl(im,j))/dx  &   
       &            +  0.25*(f(j)+zeta1(i,j))*(vv(i,j)+vv(im,j))  &
       &            +  0.25*(f(jp)+zeta1(i,jp))*(vv(i,jp)+vv(im,jp))  &
       &            +  0.25*zeta(i,j)*(vv1(i,j)+vv1(im,j))  &
       &            +  0.25*zeta(i,jp)*(vv1(i,jp)+vv1(im,jp))  &
       &            +  0.25*(1/hek)*zeta(i,j)*(vv(i,j)+vv(im,j))  &
       &            +  0.25*(1/hek)*zeta(i,jp)*(vv(i,jp)+vv(im,jp))  &
       &            -  0.5*(div(i,j)+div(im,j))* uu(i,j) &   !dns
       &            +  taux(i,j)*ramp  & 
       &            -  Ah*grad4u(i,j)        

       rhs_Vek(i,j) =  -(B(i,j)-B(i,jm))/dy  &
       &               -(1/hek)*(B_nl(i,j)-B_nl(i,jm))/dy  &
       &            -  0.25*(f(j)+zeta1(i,j))*(uu(i,j)+uu(i,jm))  &
       &            -  0.25*(f(jp)+zeta1(ip,j))*(uu(ip,j)+uu(ip,jm))  &
       &            -  0.25*zeta(i,j)*(uu1(i,j)+uu1(i,jm))  &
       &            -  0.25*zeta(ip,j)*(uu1(ip,j)+uu1(ip,jm))  &
       &            -  0.25*(1/hek)*zeta(i,j)*(uu(i,j)+uu(i,jm))  &
       &            -  0.25*(1/hek)*zeta(ip,j)*(uu(ip,j)+uu(ip,jm))  &
       &            -  0.5*(div(i,j)+div(i,jm))* vv(i,j) &   !dns
       &            -  Ah*grad4v(i,j)       


       enddo
       enddo

