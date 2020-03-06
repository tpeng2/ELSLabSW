!
!   use rhs_Psurf for w_ek (calc between ilevel 2,3)
       do j = 1, ny
       do i = 1, nx
          w_ek(i,j) =      (Uek(i+1,j,2)-Uek(i,j,2))/dx &
              &          + (Vek(i,j+1,2)-Vek(i,j,2))/dy 
       enddo
       enddo
       array(:,:) = w_ek(:,:)
       include 'subs/bndy.f90'
       w_ek(:,:) = array(:,:)
       Psurf(:,:) = Psurf(:,:)/2/dt
       

        param = (time-tcenter)/1.1/86400.
        param = param**2

      !normalization 2304 for dt = 300, window = 8 days width = 1.1

       w_filtered(:,:) = w_filtered(:,:) &
           &   + w_ek(:,:)*exp(-param)/561.512451
       p_filtered(:,:) = p_filtered(:,:) &
           &   + Psurf(:,:)*exp(-param)/561.5122451

       u1_filtered(:,:) = u1_filtered(:,:) &
           &   + u(:,:,1,2)*exp(-param)/561.512451
       u2_filtered(:,:) = u2_filtered(:,:) &
           &   + u(:,:,2,2)*exp(-param)/561.512451
       v1_filtered(:,:) = v1_filtered(:,:) &
           &   + v(:,:,1,2)*exp(-param)/561.512451
       v2_filtered(:,:) = v2_filtered(:,:) &
           &   + v(:,:,2,2)*exp(-param)/561.512451

       Uek_filtered(:,:) = Uek_filtered(:,:) &
           &   + Uek(:,:,2)*exp(-param)/561.512451
       Vek_filtered(:,:) = Vek_filtered(:,:) &
           &   + Vek(:,:,2)*exp(-param)/561.512451
       eta_filtered(:,:) = eta_filtered(:,:) &
           &   + eta(:,:,2,2)*exp(-param)/561.512451

        include 'subs/calc_wp1.f90'

        if(time.ge.tcenter.and.time.lt.tcenter+dt) then
            u1_snap(:,:) = u(:,:,1,2)
            v1_snap(:,:) = v(:,:,1,2)
            u2_snap(:,:) = u(:,:,2,2)
            v2_snap(:,:) = v(:,:,2,2)
            Uek_snap(:,:) = Uek(:,:,2)
            Vek_snap(:,:) = Vek(:,:,2)
            Psurf_snap(:,:) = Psurf(:,:)
            wek_snap(:,:) = w_ek(:,:)
            eta_snap(:,:) = eta(:,:,2,2)
        endif

        if(time.ge.tstop) then
        include 'subs/calc_wp_diss.f90'
        include 'subs/AGdecomp_hp.f90'
        Psurf_hp = Psurf_snap - p_filtered 
        Psurf_bt = Psurf_hp   &
            &    + H(2)*gprime(2)*(eta_G + eta_A)/Htot
        Psurf_A = gprime(2)*eta_A
        Psurf_G = gprime(2)*eta_G

        include 'subs/gnu/dump_w.f90'
        include 'subs/gnu/dump_p.f90'
        include 'subs/gnu/dump_u.f90'
        include 'subs/gnu/dump_v.f90'
        stop
!       include 'subs/calc_wp.f90'
!       include 'subs/calc_diss.f90'
        tstart = time + dt
        tstop = tstart + 8.*86400.
        tcenter = (tstart + tstop)/2.
        w_filtered(:,:) = 0.
        p_filtered(:,:) = 0.
        u1_filtered(:,:) = 0.
        u2_filtered(:,:) = 0.
        v1_filtered(:,:) = 0.
        v2_filtered(:,:) = 0.
        Uek_filtered(:,:) = 0.
        endif
