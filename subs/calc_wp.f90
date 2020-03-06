
       array(:,:) = w_filtered(:,:)*p_filtered(:,:)
       wp_avg(:,:) = wp_avg(:,:) + array(:,:)
       w_avg(:,:) = w_avg(:,:) + w_filtered(:,:) 
       p_avg(:,:) = p_avg(:,:) + p_filtered(:,:) 
       icount_wp = icount_wp + 1
!       include 'subs/dump_gnu.f90'
!       stop

       if(mod(icount_wp,10).eq.0) then
           open(33, file = 'WP',  status = 'unknown') 
           rewind(33)
           do j = 1,ny
           tmp(1) = 0.
           tmp(2) = 0.
           do i = 1,nx
           tmp(1) = tmp(1) + wp_avg(i,j)
           tmp(2) = tmp(2) + w_avg(i,j)*p_avg(i,j)
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
        include 'subs/dump_gnu.f90'
       endif



