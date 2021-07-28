

              EDGRN_SF/EDCMP_SF
              Massimo Nespoli
           massimo.nespoli2@unibo.it



%%%%HOW TO COMPILE EDGRN_SF/EDCMP_SF %%%%%%%%%%%%%
From the shell,
go to edcmpf77_2.0 and/or edgrnf77_2.0 folders and:

Use makefiles with "make"
---------------------------------------------------------
add -no-pie flag in makefile in case of problems with PIE:

$(PROGRAM): 	$(OBJECTS)
		$(FC) -no-pie $(FFLAGS) $(OBJECTS) -o $@
----------------------------------------------------------

		
or

alternatives for LINUX:

gfortran *.f 

or

f77 -c *.f 
f77 -o output *.o 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Run "edgrn_sf" with inputfile "edgrn_sf.inp" to generate 
the numerical Green functions. They will be saved 
in the edgrnfcts directory.


Run "edcmp_sf" with inputfile "edcmp_sf.inp" to compute
Displacement, Stress and Strain due to single forces.
The output will be saved in three different files:

izmhs.disp    displacement
izmhs.strn    strain
izmhs.strss   stress

Each inputfile contains detailed instructions
