       
       thickness(:,:) = H(1) - eta(:,:,2,3)
       array = thickness
       include 'subs/bndy.f90'
       thickness = array

       do j = 1,ny
       do i = 1,nx
          Uek_qg(i,j) = 0.
          Vek_qg(i,j) = -taux(i,j)/f0
       enddo
       enddo

       array = Uek_qg
       include 'subs/bndy.f90'
       Uek_qg = array
       array = Vek_qg
       include 'subs/bndy.f90'
       Vek_qg = array

       do j = 1,ny
       do i = 1,nx

       div_ek(i,j) = (Uek(i+1,j,3)-Uek(i,j,3))/dx   &
       &           + (Vek(i,j+1,3)-Vek(i,j,3))/dy

       div_ek_qg(i,j) = (Uek_qg(i+1,j)-Uek_qg(i,j))/dx   &
       &              + (Vek_qg(i,j+1)-Vek_qg(i,j))/dy
 
       forcing_total(i,j) = (f0+zeta1(i,j))/thickness(i,j)/thickness(i,j)*div_ek(i,j)

       forcing_qg(i,j) = (f0+zeta1(i,j))/thickness(i,j)/thickness(i,j)*div_ek_qg(i,j)
 
       forcing_ag(i,j) = forcing_total(i,j)-forcing_qg(i,j)
 
       enddo 
       enddo

       array = forcing_total
       include 'subs/bndy.f90'
       forcing_total = array       
       array = forcing_qg
       include 'subs/bndy.f90'
       forcing_qg = array
       array = forcing_ag
       include 'subs/bndy.f90'
       forcing_ag = array


