F95 = gfortran
FFLAGS = -O -C
GNUPLOT = gnuplot

%.o: %.f95
	$(F95) $(FFLAGS) -c $*.f95

objfiles = constants.o random.o vmc_input.o n_drude_wf.o vmc.o

all: vmc

vmc: $(objfiles)
	$(F95) $(FFLAGS) $(objfiles) -o vmc.x

clean:
	rm -f $(objfiles) core* *.mod *~ *.o
