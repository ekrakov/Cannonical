#  cannonical band calculation



example input

```
Structure Constants for BCC Lattice
2.5 3       17.8       1          0
3   1       3445555555 1          0
3   1.0     7.
0.5 0.5     0.5        -0.5       0.5
-0.5 0.5    -0.5       0.5        0.5
0   0       0          0          0.
0   0       17         0          0 ```







#Input files (will add)

## Structure Constant calculation

### construction of k-mesh in k space
    

### Transformation from symmetry to rectangular cords
    unsure what this is doing 

### alpha calculation???

### convergence functions

### factorial calculation
    Currently having an issue with this!
    
### cross product
    working test script 

### generate vectors of direct and real space

### final calculation of structure constant matrix 

There are currently three files:


Input
RMAX, GMAX (FIRST LINE) radii for k and q respectively (Rest zeros). 
Third line
Nl = 1-l max  (so NL=3 includes s,p, and d electrons)
nq is the number of atoms per primitive cell
NGLN?
ALAMDA? (DETERMINED WITH TEST RUNS WITH ICNVRG=1)
FOURTH LINE

LAT (SELECTS BZ TO BE USED IN CALCULATION ACCORDING TO TABLE IN SUBROUTINE MESH)
BOA B/A 
COA C/A

.....
TRANSLATIONAL VECTORS SPANNING UNIT CELL




### spherical:

describe spherical harmonic function

### factorial:

At the moment, this code is in fortran 77. To run,

```bash
gfortran -lgfortran -o factorial factorial.F
```

```bash
./factorial
```

## LMTO

## Cannonical bands 
