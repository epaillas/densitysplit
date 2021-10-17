export FC_COMPILER=gfortran

all: survey

.PHONY: survey

survey:
	make -C survey/fortran -f Makefile install

clean:
	make -C survey/fortran -f Makefile clean
	rm -f survey/bin/*.exe
