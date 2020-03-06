!      Compute wavenumbers
       do j = 1,ny
        do i = 1,nx/2+1
         if (i < (nx/2+2)) then
          kx(i,j) = twopi*(i-1)/Lx
         else
          kx(i,j) = twopi*(-nx+i-1)/Lx
         endif
         if (j < ny/2+2) then
          ky(i,j) = twopi*(j-1)/Ly 
         else
          ky(i,j) = twopi*(-ny+j-1)/Ly
         endif
        enddo
       enddo
       k2 = kx**2 + ky**2
       nkx = kx*Lx/twopi
       nky = ky*Ly/twopi
       nk2 = nkx**2 + nky**2
       nk = ifix(sqrt(float(nk2)))

!      do j = 1,ny
!      do i = 1,nx/2+1
!      if(nk2(i,j).le.1) then
!          print*, nkx(i,j), nky(i,j), i, j
!      endif
!      enddo
!      enddo
!      stop

