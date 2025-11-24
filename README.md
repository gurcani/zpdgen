zpdgen
======

A generalized plasma dispersion function formulation of Inm's, which are defined for Im(za)>0 as:

```
                              /inf /inf
Inm(za,zb,b) = (2/sqrt(pi)) * | dy | dx x**n*y**m*exp(-x**2-y**2)*J0(sqrt(2*b)*x)**2/(y**2+x**2/2+za-zb*y)
                          -inf/   0/
```

The idea is to reduce the above integral to a 1D integral:

```
                 /inf
Inm(za,zb,b) = 2*| dx x**n*exp(-x**2)*J0(sqrt(2*b)*x)**2*Gm(z1(x,za,zb),z2(x,za,zb))
                0/
```

where 

```
z1=1/2*(zb+sqrt(zb**2-2*(x**2+2*za)))
z2=1/2*(zb-sqrt(zb**2-2*(x**2+2*za)))
```

and 

```
                           /inf
Gm(z1,z2) = (1/sqrt(pi)) * | dx x**m*exp(-x**2)/((x-z1)*(x-z2))
                       -inf/
```

which is a simple generalization of the usual plasma dispersion function.

How to Compile:
============

you can compile the main library and the test binaries as follows:

```
./configure
make
```

You probably need gfortran. (i.e. `sudo pacman -Sy gcc-fortran`, or `sudo apt install gfortran`

configure script has some special flags in particular

```
--enable-python : generates the python wrapper interface.
--enable-debug : enables the debug symbols to be built.
--enable-profile : enables profiling.
```

In other words, for example, you can enable the python wrapper interface using:

./configure --enable python
make

Note that you need to install meson if you enable the python interface (e.g. `sudo pacman -Sy meson` or `sudo apt install meson ninja-build`)

How to Use:
==========
If you want to use th Inm functions in your project, you can either include zpdgen by hand or install it system wide.

In order to include it in your C/C++ project you need to include the header file "gpdf.h" and link against the generated libraries libgpdf.a, libamos.a an libquadp.a. 
You also need to add -L /path/to/zpdgen -L/path/to/zpdgen/amos -L/path/to/zpdgen/quadp as linker flags to your project.

If you want to use the python wrapper, you will need numpy (http://www.numpy.org/). After compiling and generating the python wrappers (see How to Compile), you can use it as follows. 
Here is the output from ipython:

```
In [1]: import gpdf as gp

In [2]: import numpy as nm

In [3]: za=nm.arange(-1,1,0.1)

In [4]: gp.Inm(za,0.0,0.01,1,0)
Out[4]: 
array([-0.74052221-1.41161561j, -0.65900773-1.56713712j,
       -0.54465944-1.73075998j, -0.38901067-1.89832962j,
       -0.18160146-2.06277514j,  0.09048294-2.2122345j ,
        0.44318584-2.326756j  ,  0.89620376-2.37168455j,
        1.47394180-2.28122997j,  2.20698953-1.90182412j,
        3.14186311-0.03544759j,  1.64492834+0.j        ,
        1.30212152+0.j        ,  1.10111201+0.j        ,
        0.96285427+0.j        ,  0.85984540+0.j        ,
        0.77921516+0.j        ,  0.71392125+0.j        ,
        0.65970808+0.j        ,  0.61381924+0.j        ])
```

Notice how the wrapper is able to treat arrays. While this is not a feature tested extensively (So I suggest verifying the result at least once with the result from a hand made loop in the case of a complex array shape), it should work for arrays of higer rank as well (tested with arrays with two and three dimensions).

How to Install Systemwide
=========================

```
./configure --prefix dest_path
make
make install
```

will install libgpdf.a, ligquadp.a, libamos.a in dest_path/lib and gpdf.h in dest_path/include
