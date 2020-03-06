# ELSLabSW

## Dependent files
- Parameters and input
| File name name | Description | 
| --- | --- |
| `./parameters.f90` | parameter files |

- Subroutines
| File name name | Description | 
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
| File name name | Description | 
| --- | --- |
| `./fftw_stuff/fft_params.f90` | fft params files | 
| `./fftw_stuff/fft_init.f90` | fft init files | 
| `./fftw_stuff/spec1.f90` | fft init files | 
| `./fftw_stuff/fft_destroy.f90` | fft init files | 

- Output file
| File name name | Description | 
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

## Program structure
