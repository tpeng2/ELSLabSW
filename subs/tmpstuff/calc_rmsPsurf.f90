
       tmp(1) = 0.
       xavg(:) = 0.

       do  j = 1,ny
       do  i = 1,nx
       xavg(j) = xavg(j) + Psurf(i,j)
       enddo
       xavg(j) = xavg(j)/nx
       array(:,j) = Psurf(:,j) - xavg(j)
       enddo

       xavg(:) = 0.
       do  j = 1,ny
       do  i = 1,nx
       xavg(j) = xavg(j) + Psurf(i,j)**2
       enddo
       xavg(j) = xavg(j)/nx
       write(30,*) j, xavg(j)
       enddo

       xavg(:) = 0.
       do  j = 1,ny
       do  i = 1,nx
       xavg(j) = xavg(j) + array(i,j)**2
       enddo
       xavg(j) = xavg(j)/nx
       write(31,*) j, xavg(j)
       enddo

       

