# ELSLabSW
## Module data_initial
Table of variables
| Var. name | Type | Description |
| --- | --- | --- |
| ('nx','ny','nz') | Integer |  grid number in (x,y,z) direction |
| 'nnx','nny' | Integer |  ...? |
| 'pi' | double |  const pi |
| 'twopi' | double |  const 2*'pi' |
| ('Lx','Ly') | double |  (horiz.) Domain length|
| ('dx','dy') | double |  (horiz.) Grid size|
| 'f0' | real |  Coriolis f_0|
| 'beta' | real |  Coriolis f=f_0+beta*y|
| 'r_drag' | real |  ...? |
| 'Ah' | real |  Eddy diffusivity ? |
| 'r_invLap' | real |   ...? |
| 'rf' | real |   ...? |
| 'tau0' | real |   surface stress ? |
| 'tau1' | real |    ? |
| 'hek' | real |    Ekman layer depth? |

Logical parameters
| Var. name | Type | Description |
| --- | --- | --- |
|'restart'| logical |  ... |
|'use-ramp'| logical |  ... |
|'ifsteady'| logical |  ... |

The parameter file is named as 'parameters.f90'.

## Program structure
