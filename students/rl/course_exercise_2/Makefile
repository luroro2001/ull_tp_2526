## Compilers
FC = gfortran
MPIFC   = mpif90

## Flags
FCFLAGS = -O2
OMPFLAGS = -fopenmp 

## Source files
#OBJECTS = geometry.o particle.o bh.o main.o
SRC = openmp/geometry.f90 openmp/particle.f90 openmp/bh.f90 openmp/main.f90
SRC_MPI    = mpi/geometry.f90 mpi/particle.f90 mpi/bh.f90 mpi/main_mpi.f90

## Targets (default: serial)
all: nbody_bh

# In serial, I enabled openmp just to avoid issues with 'omp_get_wtime':
nbody_bh: $(SRC)
	$(FC) $(FCFLAGS) $(OMPFLAGS) -o $@ $(SRC) 

# OpenMP version
nbody_bh_omp: $(SRC)
	$(FC) $(FCFLAGS) $(OMPFLAGS) -o $@ $(SRC)

# MPI version
mpi: nbody_bh_mpi

nbody_bh_mpi: $(SRC_MPI)
	$(MPIFC) $(FCFLAGS) -o $@ $(SRC_MPI)


# Module dependencies
#%.o: %.f90
#	$(FC) $(FCFLAGS) -c $<

## Utilities

# Remove compiled files
clean:
	rm -f *.o *.mod nbody_bh nbody_bh_omp nbody_bh_mpi output.dat

# Run tests
test: nbody_bh
	./nbody_bh < input.dat

test_omp: nbody_bh_omp
	OMP_NUM_THREADS=4 ./nbody_bh_omp < input.dat

test_mpi: nbody_bh_mpi
	mpirun -np 4 ./nbody_bh_mpi < input.dat

# Note: use 'make mpi' when using that. INCKUDE THIS LATER IN TUTORIAL