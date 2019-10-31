.DEFAULT_GOAL := all

DIR = .

TGT = main.exe

FC = mpif90

OBJDIR = obj
SRCDIR = src
INCLUDEDIR = inc
OUTPUT = output

MOD = -J$(DIR)/$(INCLUDEDIR)
WALL = -Wall
OPT = -O3
OG = -Og
DEBUG = -g -fcheck=all -fimplicit-none -fbacktrace -pedantic -Wall

SRC =  dg_param.f90 \
       dg_mpi.f90 \
       dg_local_storage.f90 \
       dg_affine_map.f90 \
       dg_basis.f90 \
       dg_interpolate_to_new_point.f90 \
       dg_ghost.f90 \
       dg_user_defined.f90 \
       dg_external_state.f90 \
       dg_poly_level_and_order.f90 \
       dg_hilbert_curve.f90 \
       dg_nodal_2d_storage.f90 \
       dg_construct_mpi_boundary.f90 \
       dg_get_dual_coord.f90 \
       dg_gen_dual_graph.f90 \
       dg_write_data.f90 \
       dg_read_mesh_2d.f90 \
       dg_hilbert_sort.f90 \
       dg_distribute_elements.f90 \
       dg_start_parallel.f90 \
       dg_prepare_hilbert_scheme.f90 \
       dg_basis_storage.f90 \
       dg_constructor.f90 \
       dg_riemann_solver.f90 \
       dg_interfaces_construct.f90 \
       dg_flux_vector.f90 \
       dg_numerical_flux.f90 \
       dg_spatial_derivative.f90 \
       dg_a_times_spatial_der.f90 \
       dg_time_derivative_global.f90 \
       dg_step_by_RK3.f90 \
       dg_output.f90 \
       dg_io.f90 	\
       dg_advection_diffusion_driver.f90 \
       dg_verification.f90 \
       dg_end_games.f90 \
       dg_main_loop.f90
#       dg_central_flux.f90 \
#       dg_Lax_Friedrichs_flux.f90 \
#       dg_time_derivative.f90 \

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
	mpirun -np 2 $(TGT)

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
	rm -rf $(OUTPUT)/*.dat
	rm -rf *.dat *.txt $(TGT)
