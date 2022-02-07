#Set the compilers and linker
# FC=gfortran
# CC=gcc
# LD=gfortran

FC=gfortran
CC=gcc
LD=gfortran

#Set the objects
OBJS= dummy_grad.o sundials_solve.o sundials_test.o

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fortran Compiler options
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MY_OPTIONS    =           -ffree-line-length-800 -o3
FC_FLAGS    =           -fcheck=bounds -ffree-line-length-0 -fimplicit-none -fbacktrace -O3 -g -DDEBUG -fopenmp
FC_FLAGS_DEBUG = -fimplicit-none -fcheck=bounds -ffree-line-length-0 -fcheck=all -fbacktrace -g -DDEBUG
FC_FLAGS_HOME = -fbounds-check -ffree-line-length-0 -DHOME -g -O3

debug:	FC_FLAGS = $(FC_FLAGS_DEBUG)
home:		FC_FLAGS = $(FC_FLAGS_HOME)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# LINKER - FORTRAN LINKER (MPI w/ GNUCOMPILER gfortran)
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

LD_FLAGS = -O3 -g

LD_FLAGS_DEBUG = -g -Og

debug:	LD_FLAGS = $(LD_FLAGS_DEBUG)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# PETSC (STATIC & DYNAMIC) LIBRARY
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# debug: PETSC_DIR = $(PETSC_DIR_DEBUG)

SUNDIALS_INCLUDE_PATH = -I/mnt/d/code/sundials/instdir/include/

SUNDIALS_LINK_PATH = -L/mnt/d/code/sundials/instdir/lib/
SUNDIALS_MOD_FILES = -I/mnt/d/code/sundials/instdir/fortran/

#--dev
SUNDIALS_LOAD_PATH = -Wl,-rpath=/mnt/d/code/sundials/instdir/lib/

# home:	PETSC_LOAD_PATH = $(PETSC_LOAD_PATH_HOME)

PETSC_LIB = -lsundials_fcvode_mod -lsundials_fsunlinsolklu_mod

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# PATHS
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

INCLUDE_PATH = $(SUNDIALS_INCLUDE_PATH)
LINK_PATH = $(SUNDIALS_LINK_PATH)
LOAD_PATH = $(SUNDIALS_LOAD_PATH)
MOD_PATH = $(SUNDIALS_MOD_FILES)
LIBS = $(PETSC_LIB)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEFINE MODS: OBJS .o REPLACED BY .mod
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
MODS= $(OBJS:.o=.mod)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEFINE EXECUTABLE
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
EXEC=test

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAKEFILE VARIABLE
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
DEFAULT=makefile

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAKE EXECUTABLE
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
all : $(EXEC)
debug : $(EXEC)
home: $(EXEC)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# LINK THE SOURCE (.o) INTO THE TARGET (EXEC) - explicit rule
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
$(EXEC): $(OBJS)
	$(LD) $(FC_FLAGS) $(LINK_PATH)  $(MOD_PATH) -o $@ $^ $(LOAD_PATH) $(LIBS)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# COMPILE INFERRED SOURCE (.f90) INTO TARGET (.o) - inference rule
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%.o : %.f90
	$(FC) $(FC_FLAGS) $(MOD_PATH) -cpp -dM -c $< -o $@ $(INCLUDE_PATH) $(LINK_PATH) $(LIBS)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# UPDATE OBJS
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%.f90 : $(DEFAULT)

#Clean the directory and any directories searched
.PHONY : clean
clean :
	rm -rf $(OBJS) $(EXEC) $(MODS) *.o *.mod
