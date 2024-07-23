# Adaptive Dynamic Relaxation

## About
Dynamic Relaxation (DR) is a niche finite element analysis (FEA) technique which allos analysis of flexible structures with very large structural deflections. This archive contains ANSI C and FORTRAN programs for carrying out adaptive DR analysis. This process can be broken down into following steps.

* Specifying a coarse 2D mesh with boundary nodes and out-of-plane loaded nodes.
* Computing the deformed coarse 3D mesh, from the loading, and elemental stresses by *DR*.
* Conducting adaptive analysis.
* Generating a refined mesh using the meshing parameters computed during adaptive analysis.
* Mapping the refined 2D adaptive mesh to coarse 3D mesh to form the final 3D refined mesh.
## Initial Mesh
![coarse 2D input mesh](/assets/PNG/walldr.png)
## Dynamic Relaxation
Copy over a *dot-dat* file and call it `input.dat`
```
cp ../dot-dat/walldr.dat 
cp walldr.dat input.dat
```
Run the DR program `flat2`
```
../flat2 
 input number of required iterations
1000
 INPUT MASS ADDITION
1000
 INPUT DENSITY
1
 OUTPUT /SLACK LENGTHS T/F
f


 INITIALISE LOAD,MASS,FOR FIRST TIME
  ENERGY AT ITERATION    1    EQUALS       0.13556818E+09
...
...
ENERGY AT PEAK   ITERATION  228    EQUALS       0.46047986E-09
  ENERGY AT ITERATION  229    EQUALS       0.42586760E-10
```
 ## View DR results
```
cp data.dat input.dat
../plotp
output device is postscript file pltf*.plt
 Brightness controll factor [0.]...[1.]
0.7
 plot membrane elements ? (t/f)
t
 plot element errors or Domain dec? (t/f)
f
 plot stress ? (t/f)
f
 plot cable elements ?    (t/f)
f
 plot g strings ?   (t/f)
f
 number nodes ?           (t/f)
f
 number elements ? (t/f)
f
 element reduction factor ? (0.0 - 1.0)
0.
 auto scale ?             (t/f)
t
 plot single view ? t/f
t
 view no   1/2/3/4 ?
4
  isoyes =   T
 input vertical angle
30
 input horizontal angle
30
    pltf4.plt     opened 
 enter the title
walldr.dat after DR
   4.6354998053636374E-310   0.0000000000000000        0.0000000000000000        0.0000000000000000
```
## Display DR Output
```
ps2pdf pltf4.plt
okular pltf4.plt.pdf
```
![DR displaced mesh walldr-o.dat](/assets/PNG/walldr-o.png)
## Adaptive Analysis
```
cp data.dat walldr-o.dat
cp stress.dat walldr-o.str
../adpgs1 
Enter the input file (.dat assumed)
walldr-o
Example rect                    
0.007000 0.000000 
opening walldr-o.str for reading the stresses
 1.11530e+05  1.22930e+04  3.18390e+04
...
...
9.17860e+03  1.14650e+05  2.63560e+04

 nita (percent) from the current mesh = 72.263493
 he_min = 1.238584  he_max = 16.746790
Enter your own hmin 
1.25
Magnification =   1.0 (hmin =   1.2)
*** Element errors are stored in walldr-o.err ***
*** Mesh parameters are in walldr-o.me 
*** Adaptivity results are in walldr-o.adp
```
## Mesh Refinement
```
../faopgs2 
Max of 50000 nodes and 70000 elements can be meshed

***Node para are going to be used*** 
Enter the project file name (without extension)
walldr-o
Example rect                    
0.007000 0.000000 
opening walldr-o.men for reading 
cannot open walldr-o.men  opening walldr-o.me for read
Enter [1] if coarse background mesh is to be refined
1
Enter [1] if mesh post-processing required
1
The program is to be run without debug

 The No. of p_stack elems = 0
Coarse BG mesh elem no.1 resulted in --> nn = 63 and ne = 94
...
...
Coarse BG mesh elem no.28 resulted in --> nn = 63 and ne = 94

 over all diagonal exchange started 

 overall smoothing started
iln = 10

 Time taken for mesh generation =     0.0000 min



  *** Mesh Compiled with 676 Nodes and 1198 Element *** 
*** File walldr-oo.dat with 676 nodes & 1198 elements has been generated ***
*** Neural net training data saved in walldr-o.inf ***
```
## Map the 2D refined mesh to its 3D geometry
```
../surfgs1 
Enter the (3D) coarse mesh file name (without extension)
walldr-o
Example rect                    
0.000000 0.000000 
Enter the (2D) refined mesh file name (without extension)
walldr-oo
Example rect                    
0.000000 0.000000 
Dmax = 141.421356
Nq = 77.781746
Higher Nq value can improve the quality of mesh
Enter Nq or RETURN to accept computed value
500
Nq is set to: 500.000000
***input.dat file has been generated***
```
## View the 3D adaptive mesh
```
../plotp 
 output device is postscript file pltf*.plt
 Brightness controll factor [0.]...[1.]
0.7
 plot membrane elements ? (t/f)
t
 plot element errors or Domain dec? (t/f)
f
 plot stress ? (t/f)
f
 plot cable elements ?    (t/f)
f
 plot g strings ?   (t/f)
f
 number nodes ?           (t/f)
f
 number elements ? (t/f)
f
 element reduction factor ? (0.0 - 1.0)
0.
 auto scale ?             (t/f)
t
 plot single view ? t/f
t
 view no   1/2/3/4 ?
4
  isoyes =   T
 input vertical angle
30
 input horizontal angle
30
    pltf4.plt     opened 
 enter the title
Adaptive Refined Mesh walldr-oo.dat     
   4.6357120049427339E-310   0.0000000000000000        0.0000000000000000        0.0000000000000000
```
## Display adaptively refined mesh in 3D
```
ps2pdf pltf4.plt
okular pltf4.plt.pdf 
```
![adapative refined mesh walldr-oo.png](/assets/PNG/walldr-oo.png)
## Remarks
The mapping of any 2D unstructured mess on to a 3D surface is explained in Khan and Topping [1]. Full details of adaptive DR are available in Topping and Khan [2].

## References
[[1](/assets/PDF/surf.pdf)] A.I. Khan and B.H.V. Topping Three Dimensional Adaptive Surface Re-Meshing For Large Displacement Finite Element Analysis.
[[2](https://www.saxe-coburg.co.uk/pubs/descrip/btak.htm)]: PARALLEL FINITE ELEMENT COMPUTATIONS B.H.V. Topping and A.I. Khan