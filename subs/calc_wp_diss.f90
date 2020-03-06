
       array(:,:) = w_filtered(:,:)*p_filtered(:,:)
       wp_avg(:,:) = wp_avg(:,:) + array(:,:)
       w_avg(:,:) = w_avg(:,:) + w_filtered(:,:) 
       p_avg(:,:) = p_avg(:,:) + p_filtered(:,:) 

       include 'subs/calc_diss.f90'
       u1diss(:,:) = u1diss(:,:)*H(1)
       u2diss(:,:) = u2diss(:,:)*H(2)
       v1diss(:,:) = v1diss(:,:)*H(1)
       v2diss(:,:) = v2diss(:,:)*H(2)

       array(:,:) = u1_filtered(:,:)*u1diss(:,:) &
           &      + u2_filtered(:,:)*u2diss(:,:) &
           &      + v1_filtered(:,:)*v1diss(:,:) &
           &      + v2_filtered(:,:)*v2diss(:,:) 

       diss_avg(:,:) = diss_avg(:,:) + array(:,:)
       u1_avg(:,:) = u1_avg(:,:) + u1_filtered(:,:)
       u2_avg(:,:) = u2_avg(:,:) + u2_filtered(:,:)
       v1_avg(:,:) = v1_avg(:,:) + v1_filtered(:,:)
       v2_avg(:,:) = v2_avg(:,:) + v2_filtered(:,:)

       icount_wp = icount_wp + 1
!      if(mod(icount_wp,10).eq.0) then
           open(33, file = 'WP',  status = 'unknown') 
           rewind(33)
           do j = 1,ny
           tmp(:) = 0.
           do i = 1,nx
           tmp(1) = tmp(1) + wp_avg(i,j)
           tmp(2) = tmp(2) + w_avg(i,j)*p_avg(i,j)
           tmp(4) = tmp(4) + diss_avg(i,j)
           enddo
           tmp(1) = tmp(1)/icount_wp/nx
           tmp(2) = tmp(2)/icount_wp/icount_wp/nx
           tmp(3) = tmp(1) - tmp(2)
           tmp(4) = tmp(4)/icount_wp/nx

           write(33,*) j, tmp(2), tmp(3), tmp(1), tmp(4)
           enddo
           close(33)
           call flush(33)
!      endif

       if(mod(icount_wp,20).eq.0) then
        include 'subs/dump_gnu.f90'
       endif

