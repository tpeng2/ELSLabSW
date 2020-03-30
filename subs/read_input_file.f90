

      subroutine read_input_file
        use data_initial    
        character(48)  label
        character(180)  params_dir
        namelist /Grid/ Lx,Ly,nx,ny,nz
        namelist /Parameters/ tau0,tau1,f0,beta,r_drag,co_rinvLap, &
        co_Ah, rf,c_bc,hek
        namelist /Time/ dt,totalday,ioutperday
        namelist /Switches/ ifsteady,restart,use_ramp


        call get_command_argument(1,params_dir)
        open(711,file=params_dir,status="old")
  
        read(711,nml=Grid)
        print Grid
  
        read(711,nml=Parameters)
        print Parameters
  
        read(711,nml=Time)
        print Time
  
        read(711,nml=Switches)
        print Switches

        dx=Lx/nx
        dy=Ly/ny
        nnx=nx+1
        nny=ny+1
        r_invLap = co_rinvLap*twopi**2/Ly**2
        Ah = co_Ah*dx**4 
        totaltime=totalday*86400
        nsteps=totaltime/dt
        iout=nsteps/totalday/ioutperday
        i_diags= ifix(86400./16/dt) 
        start_movie = 1.*nsteps/5.

      end 