
    integer*8        :: pr2c, pc2r
    real  spectrum(0:nx/2)
    double complex,   dimension(nx/2+1,ny) :: datc, c_array
    double precision, dimension(nx/2+1,ny) :: kx, ky, k2
    double precision, dimension(nx,ny)     :: datr, r_array
    integer, dimension(nx/2+1,ny) :: nkx,nky,nk2,nk
    INTEGER ikx,iky
