!! make a restart file with doubled resolution
!! for two-layer models

program double_res
    integer nx,ny
    integer nnx,nny,nz,nnx2,nny2,nz2
    parameter(nx = 512, ny = 512)
    parameter(nnx = nx+1, nny = ny+1, nnx2 = 2*nx+1, nny2 = 2*ny+1)
    integer ix,iy,iz,i,j,k,ixedge(4),iyedge(4)
    integer ix2,iy2,iz2,i2,j2,k2,ix2edge(4),iy2edge(4)
    real    Uek(0:nnx,0:nny),Vek(0:nnx,0:nny),u(0:nnx,0:nny,2),v(0:nnx,0:nny,2),eta_low(0:nnx,0:nny)
    real    Uek2(0:nnx2,0:nny2),Vek2(0:nnx2,0:nny2),u2(0:nnx2,0:nny2,2),v2(0:nnx2,0:nny2,2),eta_low2(0:nnx2,0:nny2)
    integer icount_srt,nspecfile,iftcount_srt
    real time
    character :: fname_in(25),fname_out(25)    !  Read
    open(0,file='restart')
          do j = 0,nny
          do i = 0,nnx
             read(0,*) Uek(i,j),Vek(i,j),u(i,j,1), &
                 &      u(i,j,2),v(i,j,1),v(i,j,2), &
                 &      eta_low(i,j)
          enddo
          enddo
          read(0,*) icount_srt,time,nspecfile,iftcount_srt
    close(0)
    restart_from=time/86400
    print*, 'Restart from', restart_from, 'day'
    !copy the left-lower corner if it's not on the original grids 
    do j2 = 0,nny2
        j=floor(j2/2.0)
        do i2 = 0,nnx2
            i=floor(i2/2.0)
            !exclude
            u2(i2,j2,:)=u(i,j,:)
            v2(i2,j2,:)=v(i,j,:)
            Uek2(i2,j2)=Uek(i,j)
            Vek2(i2,j2)=Uek(i,j)
            eta_low2(i2,j2)=Uek(i,j)
        end do
    end do
    !Reset for initialization
    time=0.0
    icount_srt=0
    iftcount_srt=0
    !Write output
    open(1,file='init_2xres')
        do j = 0,nny2
        do i = 0,nnx2
           write(1,*) Uek2(i,j),Vek2(i,j),u2(i,j,1), &
               &      u2(i,j,2),v2(i,j,1),v2(i,j,2), &
               &      eta_low2(i,j)
        enddo
        enddo
        write(1,*) icount_srt,time,nspecfile,iftcount_srt
    close(1)
    write(*,*) 'Initial file is set for nx=',nnx2-1
end program
