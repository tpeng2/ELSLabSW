
       array(:,:) = w_ek(:,:)*Psurf(:,:)
       wp_avg1(:,:) = wp_avg1(:,:) + array(:,:)
       w_avg1(:,:) = w_avg1(:,:) + w_ek(:,:) 
       p_avg1(:,:) = p_avg1(:,:) + Psurf(:,:) 
       icount_wp1 = icount_wp1 + 1

       if(mod(icount_wp1,1000).eq.0) then
           open(33, file = 'WP1',  status = 'unknown') 
           rewind(33)
           do j = 1,ny
           tmp(1) = 0.
           tmp(2) = 0.
           do i = 1,nx
           tmp(1) = tmp(1) + wp_avg1(i,j)
           tmp(2) = tmp(2) + w_avg1(i,j)*p_avg1(i,j)
           enddo
           tmp(1) = tmp(1)/icount_wp1/nx
           tmp(2) = tmp(2)/icount_wp1/icount_wp1/nx
           tmp(3) = tmp(1) - tmp(2)

           write(33,*) j, tmp(2), tmp(3), tmp(1)
           enddo
           close(33)
           call flush(33)
       endif



