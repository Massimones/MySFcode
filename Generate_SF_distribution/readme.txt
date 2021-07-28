%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The scripts generate a file 'EDCMP_sf_input.txt' containing the distribution of forces
in a compatible format for the "edcmp_sf.inp" file.

Just cut and paste the content of 'EDCMP_sf_input.txt' in the "edcmp_sf.inp" files, below 

#===============================================================================
NUMBER OF SF
#         coord. origin: (40.739N, 30.05E)
#-------------------------------------------------------------------------------
# no  Force (N)        xs     ys     zs   azimuth(to east)  inclination(depth>0)
#-------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The script "prepare_general_shape_4_0.m" requires the following matlab addons:


- stlwrite (provided) !! Please use the function provided !!

- inpolyhedron (provided)

- Partial Differential Equation Toolbox (https://it.mathworks.com/products/pde.html)


