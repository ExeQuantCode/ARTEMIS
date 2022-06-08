##########################################
# CODE DIRECTORIES AND FILES
##########################################
mkfile_path := $(abspath $(firstword $(MAKEFILE_LIST)))
mkfile_dir := $(dir $(mkfile_path))
BIN_DIR := ./bin
SRC_DIR := ./src
LIB_DIR := ./lib
BUILD_DIR = ./obj
LIBS := mod_constants.f90 \
	mod_misc.f90 \
	mod_misc_maths.f90 \
	mod_misc_linalg.f90 \
	mod_tools_infile.f90 \
	mod_rw_geom.f90 \
	mod_edit_geom.f90 \
	mod_sym.f90
OBJS := $(addprefix $(LIB_DIR)/,$(LIBS))
#$(info VAR is $(OBJS))
SRCS := io.F90 \
	aspect.f90 \
	mod_help.f90 \
	mod_intf_identifier.f90 \
	mod_plane_matching.f90 \
	mod_lat_compare.f90 \
	mod_swapping.f90 \
	mod_shifting.f90 \
	default_infile.f90 \
	inputs.f90 \
	interfaces.f90 \
        main.f90
SRCS := $(OBJS) $(SRCS)
OBJS := $(addprefix $(SRC_DIR)/,$(SRCS))


##########################################
# COMPILER CHOICE SECTION
##########################################
FFLAGS = -O2
#PPFLAGS = -cpp
FC?=gfortran
ifeq ($(FC),ifort)
	MPIFLAG = -qopenmp
	MODULEFLAG = -module
	DEVFLAGS = -check all -warn #all
	DEBUGFLAGS = -check all -fpe0 -warn -tracekback -debug extended # -check bounds
else
	MPIFLAG = -fopenmp
	MODULEFLAG = -J
	DEVFLAGS = -g -fbacktrace -fcheck=all #-g -static -ffpe-trap=invalid
	DEBUGFLAGS = -fbounds-check -Wall -Wno-maybe-uninitialized
endif


##########################################
# LAPACK SECTION
##########################################
MKLROOT?="/usr/local/intel/parallel_studio_xe_2017/compilers_and_libraries_2017/linux/mkl/lib/intel64_lin"
LLAPACK = $(MKLROOT)/libmkl_lapack95_lp64.a \
	-Wl,--start-group \
	$(MKLROOT)/libmkl_intel_lp64.a \
	$(MKLROOT)/libmkl_sequential.a \
	$(MKLROOT)/libmkl_core.a \
	-Wl,--end-group \
	-lpthread

#$(MKLROOT)/libmkl_scalapack_lp64.a \
#$(MKLROOT)/libmkl_solver_lp64_sequential.a \


##########################################
# COMPILATION SECTION
##########################################
INSTALL_DIR?=$(HOME)/bin
ARTEMIS = artemis
programs = $(BIN_DIR)/$(ARTEMIS)

.PHONY: all debug install uninstall dev mpi clean

all: $(programs)

$(BIN_DIR):
	mkdir -p $@

$(BUILD_DIR):
	mkdir -p $@

$(BIN_DIR)/$(ARTEMIS): $(OBJS) | $(BIN_DIR) $(BUILD_DIR)
	$(FC) $(MODULEFLAG) $(BUILD_DIR) $(OBJS) -o $@

install: $(OBJS) | $(INSTALL_DIR) $(BUILD_DIR)
	$(FC) $(MODULEFLAG) $(BUILD_DIR) $(OBJS) -o $(INSTALL_DIR)/$(ARTEMIS)

debug: $(OBJS) | $(BIN_DIR) $(BUILD_DIR)
	$(FC) $(DEBUGFLAGS) $(MODULEFLAG) $(BUILD_DIR) $(OBJS) -o $(programs)

dev: $(OBJS) | $(BIN_DIR) $(BUILD_DIR)
	$(FC) $(DEVFLAGS) $(MODULEFLAG) $(BUILD_DIR) $(OBJS) -o $(programs)

mpi: $(OBJS) | $(BIN_DIR) $(BUILD_DIR)
	$(FC) $(MPIFLAG) $(MODULEFLAG) $(BUILD_DIR) $(OBJS) -o $(programs)

clean: $(BUILD_DIR) $(BIN_DIR)
	rm -r $(BUILD_DIR)/ $(BIN_DIR)/

uninstall: $(INSTALL_DIR)/$(ARTEMIS)
	rm $(INSTALL_DIR)/$(ARTEMIS)
