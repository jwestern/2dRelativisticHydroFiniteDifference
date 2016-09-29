PROG =	wave

OBJS =	boundary.o derivs.o evolve.o initial.o main.o params.o rhs.o fourier_repackage.o fourier_tools.o

MAKEFLAGS = -r

LIBS = -L/usr/local/rnpletal/lib -lbbhutil -lrnpl -lutilio -lvutil -lxvs -L$(HOME)/lib/ -lfftpack5.1

F90 = gfortran 
F77 = gfortran 
F90FLAGS = -fdefault-real-8 -O3 -ffree-line-length-165 # -g -pg -p
# -fpstkchk -fpe0 -ftz
#LDFLAGS = -L/opt/intel/fce/9.0/lib
LIBPATH = -L/usr/local/lib/

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod *.o

.SUFFIXES: .o .mod .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

.f90.mod:
	$(F90) $(F90FLAGS) -c $<

.f.o:
	$(F90) $(F90FLAGS) -c $<

.f.mod:
	$(F90) $(F90FLAGS) -c $<

boundary.o: params.o
derivs.o: params.o
evolve.o: boundary.o derivs.o params.o rhs.o
initial.o: derivs.o
main.o: derivs.o evolve.o initial.o params.o
rhs.o: params.o fourier_repackage.o fourier_tools.o
fourier_tools.o: fourier_repackage.o
