

                EFGRN/EFCMP
              Massimo Nespoli
           massimo.nespoli2@unibo.it



%%%%HOW TO COMPILE EFGRN/EFCMP %%%%%%%%%%%%%
From the shell,
go to efcmpf77 and/or efgrnf77 folders and:

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

Run "efgrn" with inputfile "efgrn.inp" to generate 
the numerical Green functions. They will be saved 
in the efgrnfcts directory.


Run "efcmp" with inputfile "efcmp.inp" to compute
Displacement, Stress and Strain due to single forces.
The output will be saved in three different files:

izmhs.disp    displacement
izmhs.strn    strain
izmhs.strss   stress

Each inputfile contains detailed instructions
