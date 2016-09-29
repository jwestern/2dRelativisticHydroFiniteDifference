PROG =	wave

OBJS = get_mkl_vsl.o boundary.o derivs.o evolve.o initial.o main.o params.o rhs.o fourier_tools.o poisson_solver.o scaling_gen.o

MAKEFLAGS = -r

LIBS = -mkl -L/usr/local/rnpletal/lib -lrnpl -lutilio -lvutil -lxvs -lbbhutil
#-L/lib/x86_64-linux-gnu -lfftpack5.1d

F90 = ifort 
F77 = ifort
#F90 = gfortran
#F77 = gfortran
F90FLAGS = -r8 -O3 -132 -heap-arrays -traceback -g # -check all
#-gen-interfaces -warn all -r8 -O3 -132 
#F90FLAGS = -fbounds-check -fdefault-real-8 -O3 -ffree-line-length-165 # 
#-g3 -pg -p
# -fbacktrace -ffp3-trap=zero,overflow,underflow # -fpstkchk -fpe0 -ftz
#LDFLAGS = -L/opt/intel/fce/9.0/lib
LIBPATH = -L/usr/local/lib/

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o  $@ $(OBJS) $(LIBS)

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
evolve.o: boundary.o derivs.o params.o rhs.o poisson_solver.o fourier_tools.o
initial.o: derivs.o
main.o: scaling_gen.o derivs.o evolve.o initial.o params.o poisson_solver.o fourier_tools.o
rhs.o: params.o fourier_tools.o #fourier_repackage.o
scaling_gen.o:
fourier_tools.o: #fourier_repackage.o
