
       WRITE(which,'(I3)') 100 + icount
       string1 = 'data/zeta1'  // '_' // which(1:3)
       string2 = 'data/div1'  // '_' // which(1:3)
       string3 = 'data/zeta2'  // '_' // which(1:3)
       string4 = 'data/div2'  // '_' // which(1:3)
       open(14, file = string1, status = 'unknown')
       open(15, file = string2, status = 'unknown')
       open(16, file = string3, status = 'unknown')
       open(17, file = string4, status = 'unknown')


       do j = 14,17  !write header on files
       write(j,*) 'set pm3d map'
       write(j,*) 'splot "-"'
       enddo

       uu(:,:) = u(:,:,1,3)
       vv(:,:) = v(:,:,1,3)

       do j = 1,ny
       do i = 1,nx
       div(i,j) = (uu(i+1,j)-uu(i,j))/dx + &
           &      (vv(i,j+1)-vv(i,j))/dy
       zeta(i,j) = (vv(i,j)-vv(i-1,j))/dx - &
           &       (uu(i,j)-uu(i,j-1))/dy
       enddo
       enddo

       array = zeta
       include 'subs/bndy.f90'
       zeta = array

       array = div
       include 'subs/bndy.f90'
       div = array

       pdf(:) = 1.e-8 ! ie., zero, but for plotting on log axis
       do j = 1,ny
       do i = 1,nx
       tmp(1) = zeta(i,j)/f0 
       tmp(1) = tmp(1)*100.
       k = int(tmp(1)+0.5)
       pdf(k) = pdf(k) + 1.
       enddo
       enddo
       pdf(:) = pdf(:)/nx/ny
       open (70, file = 'pdf_zeta1')
       do i = -100,100
       write(70,*) float(i)/100., pdf(i)
       enddo
       close(70)
       
       do j = 1,ny
       y =  -Ly/2. + (j-1)*dx

       do i = 1,nx
       x =  -Lx/2. + (i-1)*dx
       write(14,*) x, y,  zeta(i,j)/f0
       write(15,*) x, y,  div(i,j)/f0
       enddo  ! end of row

       write(14,*) ! write blank line between rows 
       write(15,*) 
       enddo  ! end of j loop

       uu(:,:) = u(:,:,2,3)
       vv(:,:) = v(:,:,2,3)

       do j = 1,ny
       do i = 1,nx
       div(i,j) = (uu(i+1,j)-uu(i,j))/dx + &
           &      (vv(i,j+1)-vv(i,j))/dy
       zeta(i,j) = (vv(i,j)-vv(i-1,j))/dx - &
           &       (uu(i,j)-uu(i,j-1))/dy
       enddo
       enddo


       array = zeta
       include 'subs/bndy.f90'
       zeta = array

       array = div
       include 'subs/bndy.f90'
       div = array

       do j = 1,ny
       y =  -Ly/2. + (j-1)*dx

       do i = 1,nx
       x =  -Lx/2. + (i-1)*dx
       write(16,*) x, y,  zeta(i,j)/f0
       write(17,*) x, y,  div(i,j)/f0
       enddo  ! end of row

       write(16,*) ! write blank line between rows 
       write(17,*) 
       enddo  ! end of j loop

       do j = 14,17  !forces write from buffer
       flush(j) 
       enddo
       do j = 14,17
       close(j)
       enddo

