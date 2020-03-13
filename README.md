# ELSLabSW
Slab model + Two-layer shallow water equations

| Model name | Coupled? (alpha) | self-advective  |
| --- | --- | --- |
| T | No (alpha=0) | N/A  |
| C1 | Yes (alpha=1) |  No (`hek` set to a very large number)  |
| C2 | Yes (alpha=1) | Yes  |

## Dependent files
- Parameters and input

| File  name | Description | 
| --- | --- |
| `./parameters.f90` | parameter files |
| `amp_matrix` | Transient forcing (psuedo-stochastic)  |

- Subroutines

| File  name | Description | 
| --- | --- |
| `./subs/initialize.f90` | Initialization files | 
| `./subs/rhs_ek.f90` | Computing RHS for Ekman Layer | 
| `./subs/rhs.f90` | Computing RHS for SWE | 
| `./subs/p_correction.f90` | Correct pressure, using `fftw_stuff/invLaplacian.f90` | 
| `./subs/bndy.f90` | Periodic boundary conditions, temp. var `array`, call it everytime if necessary| 
| `./subs/div_vort.f90` | ... files | 
| `./subs/rhs.dump_gnu1a` | ... files | 
| `./subs/rhs.diags` | ... files | 
| `./subs/rhs.dump_bin` | ... files | 

- FFT library

| File  name | Description | 
| --- | --- |
| `./fftw_stuff/fft_params.f90` | fft params files | 
| `./fftw_stuff/fft_init.f90` | fft init files | 
| `./fftw_stuff/spec1.f90` | fft init files | 
| `./fftw_stuff/fft_destroy.f90` | fft init files | 

- Output file

| File name  | Description | 
| --- | --- | 
| `./ke1_spec` | Output | 
| `./ke_ek_spec` | Output | 
| `./ke2_spec` | Output | 

## Module data_initial
Table of variables
| Var. name | Type | Description |
| --- | --- | --- |
| (`nx`,`ny`,`nz`) | Integer |  grid number in (x,y,z) direction |
| `nnx`,`nny` | Integer |  ...? |
| `pi` | double |  const pi |
| `twopi` | double |  const 2*`pi` |
| (`Lx`,`Ly`) | double |  (horiz.) Domain length|
| (`dx`,`dy`) | double |  (horiz.) Grid size|
| `f0` | real |  Coriolis f_0|
| `beta` | real |  Coriolis f=f_0+beta*y|
| `r_drag` | real |  ...? |
| `Ah` | real |  Eddy diffusivity ? |
| `r_invLap` | real |   ...? |
| `rf` | real |   ...? |
| `tau0` | real |   surface stress ? |
| `tau1` | real |    ? |
| `hek` | real |    Ekman layer depth? |

Logical parameters
| Var. name | Type | Description |
| --- | --- | --- |
|`restart`| logical |  ... |
|`use-ramp`| logical |  ... |
|`ifsteady`| logical |  ... |

The parameter file is named as `parameters.f90`.
Misc.
## Parameter file setting (with example values)
| Parameter |  Value |  Note  | 
| --- | --- | --- | 
| `pi` | 2.*asin(1.) | pi | 
| `twopi` | 2.*pi | 2*pi |

Grid
| Parameter |  Value |  Note/Unit  | 
| --- | --- | --- | 
| `Lx` | 2.0e6 | [m] | 
| `Ly` | Lx | [m] | 
| `nx` | 2**9 |  | 
| `ny` | 2**9  |  | 
| `nz` | 2 |  | 
| `dx` | Lx/nx |  | 
| `dy` | Ly/ny |  | 
| `nnx` | nx+1 |  | 
| `nny` | ny+1 |  | 

Parameters
| Parameter |  Value |  Note/Unit  | 
| --- | --- | --- | 
| `tau0` |  1.e-4 |  stress  | 
| `tau1` |  1.e-5 |  stress  | 
| `f0` |  7.e-5 |  Coriolis: f_0  | 
| `beta` |  0 |  Coriolis: beta  | 
| `r_drag` |  1.e-7 |  Linear Drag from bottom  | 
| `r_invLap` |  1.e-6*twopi**2/Ly**2 |  Inverse-Laplacian damping coeff.  | 
| `Ah` |  1.e-5*dx**4 |  eddy viscosity  | 
| `rf` |  0.001 |    | 
| `c_bc` |  2. |    | 
| `h_ek` |  50. |  Ekman layer depth   | 

Time
|Parameter |  Value |  Note/Unit  | 
| --- | --- | --- | 
| `dt` |  300 |  time step size [s]  | 
| `total time` |  86400 * 9 |  total time (1 day in sec * days) [s]  | 
| ` nsteps` |  totaltime/dt |  how many steps  | 

I/O
|Parameter |  Value |  Note/Unit  | 
| --- | --- | --- | 
| `iout` |  nsteps/5 |  when to generate output files | 
| `i_diags` |  ifix(86400./16/dt) |  time step size [s]  | 
| `start_movie` |  7.*nsteps/6. |  time step when generating snapshots. no movie when > nsteps  | 
| `ifsteady` |  .true. |  .true. for steady forcing, .false. requires output  | 
| `restart` |  .false. |  if restart  | 
| `use_ramp` |  .false. |  if use ramp  | 


## Discretizations
### Spatial: Finite difference with staggered grids
See `./subs/rhs.f90` and `./main.f90`

* Relative location given indices (i,j)

|   | iu | iv | 
| --- | --- |  --- | 
| ju | `u`| `eta` (div)|
| jv |  `zeta` | `v` |

- variables in the same row share the same `j` index
- variables in the same line share the same `i` index


### Temporal: Leap frog
Central time central space
| time (n) / location (m) | *m-1* | *m* | *m+1* |
| --- | --- | --- |--- |
| *n-1* | o |  `+` |  o |
|  *n*  | `+` |  `*` |  `+` |
| *n+1* | o |  `X` |  o |

- Current variable: `*`
- Dependent variable: `+`
- Predicting variable: `X`

## Program structure
### 1. Configuration
- Calculate Rossiby length
- Determine time step
- Determine if the forcing is steady or not, if transient read `amp_matrix`.

### 2. Initialization  
Initialize following variables:
* pressure  
* thickness for two layers
* horizontal u, `uu(:.:)`
* meridional v, `vv(:.:)`
* old horizontal u, `uu_old(:.:)`
* old meridional v, `vv_old(:.:)`
* Slab current U_ek,
* Ekman transport V_ek

Then, calculate RHS for the first time step.

Finally, correct u,v for surface pressure: call    `subs/p_correction.f90`.

### 3. Subsequent time steps
#### Computing RHS of *u*, *v*, and *eta*
* *u* and *v* in Bernourlli forms: e.g., `rhs_u`= `dB/dx`+`(f+zeta)*v` + `Biharmonic hyperviscosity` + `Inverse Laplacian for damping low-frequency (barotropic mode)` + `linear drag from bottom (when k=2)` 
* *eta* equation: `rhs_eta` = `-(d(u_h)/dx+-d(v_h)/dy)`+ `stress (body-force, when alpha=1)`

#### pressure correction
- Define: Interface height for two layers *H(k)* and *k=1,2*
- Then *thickness = H(k) - eta* for the first layer, and *H(k)+eta* for the second layer.

## I/O files
### Writing Gnuplot files *./subs/dump_gnu.f90*
