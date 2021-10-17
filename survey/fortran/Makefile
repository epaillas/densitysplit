MLIBS = -lm

#################

EXES = tophat_filter.exe

all: $(EXES)

tophat_filter.exe: tophat_filter.f90
ifeq ($(FC_COMPILER), ifort)
	$(FC_COMPILER) -qopenmp procedures.f90 tophat_filter.f90 -o tophat_filter.exe $(MLIBS)
else
	$(FC_COMPILER) -fopenmp procedures.f90 tophat_filter.f90 -o tophat_filter.exe $(MLIBS)
endif


install: $(EXES)
	mv $(EXES) ../bin

clean:
	rm -f $(EXES) *.o *~ *.mod
