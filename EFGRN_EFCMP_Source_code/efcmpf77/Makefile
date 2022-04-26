FORFILES =  edcdisc.f    edcgrn.f     edcmain.f    edcoutput.f  getdata.f  \
            mind.f      mindlin.f       

OBJECTS =	$(FORFILES:%.f=%.o)

#FFLAGS =	-O3

PROGRAM =	efcmp

$(PROGRAM): 	$(OBJECTS)
		$(FC) $(FFLAGS) $(OBJECTS) -o $@

edcdisc.o:edcglobal.h edcdisc.f
edcgrn.o:edcglobal.h edcgrn.f
edcmain.o:edcglobal.h edcmain.f
edcoutput.o:edcglobal.h edcoutput.f
getdata.o:getdata.f
mind.o:mind.f
mindlin.o:mindlin.f






