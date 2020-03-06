
!
!    bndy condition: periodic in x,y
!

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

       do j = 1, ny
          jp = j+1
          jm = j-1
       do i = 1, nx
          ip = i+1
          im = i-1

       array(i,j) = (uh(ip,j)-uh(i,j))/dx   &
       &          + (vh(i,jp)-vh(i,j))/dy 

       enddo
       enddo


