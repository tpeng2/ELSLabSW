

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

  parameter ( Ah = 1.e-5*dx**4 )

  parameter ( rf = 0.001 )

  parameter ( c_bc = 2. )

  parameter ( hek = 50.)   

 ! ---  Time ---

  parameter ( dt = 300. )
 
  parameter ( totaltime = 86400 * 100. )

  parameter ( nsteps = totaltime/dt )

! parameter ( iout = 9 , i_diags = ifix(86400./16/dt) )
  parameter ( iout = nsteps/1000/16 , i_diags = ifix(86400./16/dt) )

  parameter ( start_movie = 1. , subsmpstep=4 )

  parameter ( ifsteady = .true. )

  parameter ( restart = .false. , use_ramp = .false. )
!  parameter ( restart = .false. , use_ramp = .false. )

