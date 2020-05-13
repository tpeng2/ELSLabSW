

 ! --- Misc ---

  parameter ( pi = 2.*asin(1.), twopi = 2.*pi )

 ! ---  Grid ---

  parameter ( Lx = 2.0e6, Ly = Lx )

  parameter ( nx = 2**8, ny = 2**8, nz = 2 )

  parameter ( dx = Lx/nx, dy = Ly/ny )

  parameter ( nnx = nx+1, nny = ny+1 )


 ! --- Paraterers ---

  parameter ( tau0 = 1.e-4, tau1 = 1.e-5 )

  parameter ( f0 = 7.e-5, beta = 0. )

  parameter ( r_drag = 1.e-7 )

  parameter ( r_invLap = 1.e-6*twopi**2/Ly**2 )

  parameter ( Ah = 0.05e-5*dx**4 )

  parameter ( rf = 0.001 )

  parameter ( c_bc = 2. )

  parameter ( hek = 50.)   

 ! ---  Time ---

  parameter ( dt = 300. )
 
  parameter ( ndays= 1200., totaltime = 86400 * ndays )

  parameter ( nsteps = totaltime/dt,fileperday=12 )

! parameter ( iout = 9 , i_diags = ifix(86400./16/dt) )
  parameter ( iout = nsteps/ndays/fileperday, i_diags = ifix(86400./16/dt))
 
  parameter (itape=86400*10/dt,ispechst=12) !spectrum output file, output one spectra per ispechst

  parameter(save2dfft=.true.,calc1Dspec=.false. )

  parameter ( start_movie = 1. , subsmprto=4, ftsubsmprto=2, save_movie=.true. )

  parameter ( ifsteady = .true., forcingtype=1) ! forcingtype =1, sine wave, =2 O-U process (AR1)

  parameter ( restart = .true. , use_ramp = .false. )
!  parameter ( restart = .false. , use_ramp = .false. )

