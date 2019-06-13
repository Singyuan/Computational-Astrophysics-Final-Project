Computational Astrophysics Final Project
===

## One-dimension hydrodynamics
This is simple simulation about 1D hydrodynamics.
1. Visit file "maincode".
2. The constant parameter can be tuned in "ConstPara.cpp".
3. Data reconstruction PLM or PPM can be chosen in "main.cpp".
4. Compile and excutable "main.cpp" and output the text "anidata.txt".
5. Visit file "pyani" and excutable the "anicode.py" and it will  output the animation of one-dimension Hydrodynamics.


## One-dimension magnetohydrodynamics
1. Visit file "mainMHD".
2. Compile and excutable "main_MHD.cpp"

## Two-dimension hydrodynamics
This is use unsplit mathod so it might be a little ustable
1. Visit file "main3d"
2. The constant parameter can be tuned in "ConstPara.cpp". Note that this is 3D algorithm but the z component is commennted out.
3. Choose "density", "velocity" or "pressure" data you want to make animation.
4. Compile and excutable "main3D2.cpp" and output the text "anidata.txt".
5. Run matlab code "shock3d.m" and it will  output the animation.

## Two-dimension magnetohydrodynamics
1. Visit file "maincode".
2. Run the matlab code "Hydro2D.m"

## One-dimension hydrodynamics with OpenMP
It is similar to "One-dimension hydrodynamics"

## Two-dimension hydrodynamics with MPI
It is similar to "Two-dimension hydrodynamics"
