.DEFAULT_GOAL := all

DIR = .

TGT = main.exe

FC = mpif90

OBJDIR = obj
SRCDIR = src
INCLUDEDIR = inc

MOD = -J$(DIR)/$(INCLUDEDIR)
WALL = -Wall
OPT = -O3
OG = -Og
DEBUG = -g -fcheck=all -fimplicit-none -fbacktrace -pedantic -Wall

SRC =  dg_param.f90 \
       dg_mpi.f90 \
       dg_basis.f90 \
       dg_nodal_2d_storage.f90 \
       dg_constructor.f90 \
       dg_spatial_derivative.f90 \
       dg_riemann_solver.f90 \
       dg_flux_vector.f90 \
       dg_external_state.f90 \
       dg_time_derivative.f90 \
       dg_step_by_RK3.f90 \
       dg_user_defined.f90 \
       dg_advection_diffusion_driver.f90 \
       dg_verification.f90 \
       dg_end_games.f90 \
       dg_main_loop.f90

#SOURCE = $(wildcard $(DIR)/$(makefile formatSRCDIR)/$(SRC))
#SOURCE  = $(DIR)/$(SRCDIR)
SOURCE = $(patsubst %, $(DIR)/$(SRCDIR)/%, $(SRC))

OBJ = $(addprefix $(DIR)/$(OBJDIR)/, $(notdir $(SRC:.f90=.o)))


$(DIR)/$(OBJDIR)/$(TGT) : $(OBJ)
#	$(FC) $(MOD) -o $@ $^
	$(FC) $(OPT) $(WALL) $(MOD) -o $(TGT) $^
 
$(DIR)/$(OBJDIR)/%.o : $(DIR)/$(SRCDIR)/%.f90
	$(FC) $(OPT) $(WALL) $(MOD) -c $< -o $@



.PHONY : help run clean all 

all : $(DIR)/$(OBJDIR)/$(TGT)
#	$(DIR)/$(OBJDIR)/$(TGT)
	@echo "------------------------------"
	@echo "Makefile succeed"
	@echo "------------------------------"

run : $(TGT)
	mpirun -np 1 $(TGT)

drun : $(TGT)
	mpirun -np 1 xterm -e gdb $(TGT)

help : 
	@echo "source : $(SOURCE)"
	@echo "src : $(SRC)"
	@echo "obj : $(OBJ)"


debug : 
	make "OPT = -g -fcheck=all -fimplicit-none -fbacktrace -pedantic -Wall"

clean :
	rm -rf $(OBJ) 
	rm -rf $(DIR)/$(INCLUDEDIR)/*.mod
	rm -rf *.dat *.txt $(TGT)
