F95 = gfortran
FFLAGS = -O -C
GNUPLOT = gnuplot

%.o: %.f95
	$(F95) $(FFLAGS) -c $*.f95

objfiles = n_drude_wf.o n_drude_test.o

all: out.txt

out.txt: main makefile
	./main.x > out.txt
	cat out.txt

main: $(objfiles)
	$(F95) $(FFLAGS) $(objfiles) -o main.x

clean:
	rm -f $(objfiles) core* main.x *.mod *~
