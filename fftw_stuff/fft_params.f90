    integer*8 :: pr2c, pc2r,vr2cx,vr2cy,vc2rx,vc2ry
    real  spectrum(0:nx),nkoutput(0:nx)
    real  spectrumx(0:nx/2),spectrumy(0:ny/2)
    double complex,   dimension(nx/2+1,ny) :: datc, c_array
    double complex,   dimension(nx) :: datcx, cx_array
    double complex,   dimension(ny) :: datcy, cy_array
    double precision, dimension(nx/2+1,ny) :: kx, ky, k2
    double precision, dimension(nx) :: datrx,rx_array
    double precision, dimension(nx/2+1) :: kx1d

    double precision, dimension(ny) :: datry,ry_array
    double precision, dimension(ny/2+1) :: ky1d

    double precision, dimension(nx,ny)     :: datr, r_array
    integer, dimension(nx/2+1,ny) :: nkx,nky
    integer, dimension(nx/2+1) :: nkx1d ! for 1-D, nk1 for x, nk2 for y
    integer, dimension(ny/2+1) :: nky1d ! for 1-D, nk1 for x, nk2 for y
    integer, dimension(nx/2+1,ny) :: nk2,nk
    real, dimension(nx/2+1,ny) :: nkf
    INTEGER ikx,iky
    double complex M