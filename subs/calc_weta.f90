
       array(:,:) = w_f1(:,:)*eta_f1(:,:)
       weta_avg1(:,:) = weta_avg1(:,:) + array(:,:)
       w_avg1(:,:) = w_avg1(:,:) + w_f1(:,:) 
       eta_avg1(:,:) = eta_avg1(:,:) + eta_f1(:,:) 
       icount_wp1 = icount_wp1 + 1

       if(mod(icount_wp1,5).eq.0) then
           open(33, file = 'Weta',  status = 'unknown') 
           rewind(33)
           do j = 1,ny
           tmp(1) = 0.
           tmp(2) = 0.
           do i = 1,nx
           tmp(1) = tmp(1) + weta_avg1(i,j)
           tmp(2) = tmp(2) + w_avg1(i,j)*eta_avg1(i,j)
           enddo
           tmp(1) = tmp(1)/icount_wp/nx
           tmp(2) = tmp(2)/icount_wp/icount_wp/nx
           tmp(3) = tmp(1) - tmp(2)

           write(33,*) j, tmp(2), tmp(3), tmp(1)
           enddo
           close(33)
           call flush(33)
       endif
       if(mod(icount_wp,20).eq.0) then
        include 'subs/dump1_gnu.f90'
       endif



