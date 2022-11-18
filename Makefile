# Replace mpiicpc (intel mpi compiler) with the compiler used in your system.
# Put appropriate compiler options in CXXFLAGS.
# Please be sure that the following environment variables are defined:
# 	METIS_LIB		- METIS library directory
#	PARMETIS_LIB	- ParMETIS library directory
#	PETSC_LIB		- PETSc library directory
#	IDEA_LIB		- IDEA library directory
#
#	METIS_INC		- METIS include (metis.h) directory
#	PARMETIS_INC	- ParMETIS include (parmetis.h) directory
#	PETSC_INC		- PETSc include (petsc.h) directory
#	PETSC_CONF_INC	- PETSc build-specific include (e.g., petscconf.h) directory
#	IDEA_INC		- IDEA include (idea.h) directory

CXX=mpiicpc
CXXFLAGS=-std=c++14 -w -Ofast -mkl=sequential

OBJ_DIR=obj
BIN_DIR=bin
SRC_DIR=source
AVOCADO=Avocado
DENEB=Deneb

A_OBJ_DIR=$(OBJ_DIR)/$(AVOCADO)
A_SRC_DIR=$(SRC_DIR)/$(AVOCADO)/src
A_INC_DIR=$(SRC_DIR)/$(AVOCADO)/inc
D_OBJ_DIR=$(OBJ_DIR)/$(DENEB)
D_SRC_DIR=$(SRC_DIR)/$(DENEB)/src
D_INC_DIR=$(SRC_DIR)/$(DENEB)/inc

AVOCADO_OBJ=$(A_OBJ_DIR)/avocado.o\
            $(A_OBJ_DIR)/avocado_argument.o\
            $(A_OBJ_DIR)/avocado_blas.o\
            $(A_OBJ_DIR)/avocado_config.o\
            $(A_OBJ_DIR)/avocado_dual_number.o\
            $(A_OBJ_DIR)/avocado_file.o\
            $(A_OBJ_DIR)/avocado_memory_profiler.o\
            $(A_OBJ_DIR)/avocado_mpi.o\
            $(A_OBJ_DIR)/avocado_string.o\
            $(A_OBJ_DIR)/avocado_timer.o\

DENEB_OBJ=$(D_OBJ_DIR)/deneb.o\
          $(D_OBJ_DIR)/deneb_artificial_viscosity.o\
          $(D_OBJ_DIR)/deneb_basis.o\
          $(D_OBJ_DIR)/deneb_contour.o\
          $(D_OBJ_DIR)/deneb_data.o\
          $(D_OBJ_DIR)/deneb_DRM.o\
          $(D_OBJ_DIR)/deneb_element.o\
          $(D_OBJ_DIR)/deneb_equation.o\
          $(D_OBJ_DIR)/deneb_equation_equilibriumns2d.o\
          $(D_OBJ_DIR)/deneb_equation_euler2d.o\
          $(D_OBJ_DIR)/deneb_equation_glmmhd2d.o\
          $(D_OBJ_DIR)/deneb_equation_ns2d.o\
          $(D_OBJ_DIR)/deneb_equation_ns3d.o\
          $(D_OBJ_DIR)/deneb_grid_builder.o\
          $(D_OBJ_DIR)/deneb_grid_reader.o\
          $(D_OBJ_DIR)/deneb_jacobian.o\
          $(D_OBJ_DIR)/deneb_limiter.o\
          $(D_OBJ_DIR)/deneb_point.o\
          $(D_OBJ_DIR)/deneb_polynomial.o\
          $(D_OBJ_DIR)/deneb_pressurefix.o\
          $(D_OBJ_DIR)/deneb_quadrature.o\
          $(D_OBJ_DIR)/deneb_saveload.o\
          $(D_OBJ_DIR)/deneb_system_matrix.o\
          $(D_OBJ_DIR)/deneb_timescheme.o\
          $(D_OBJ_DIR)/deneb_timescheme_impeuler.o\
          $(D_OBJ_DIR)/deneb_timescheme_rosenbrock.o\
          $(D_OBJ_DIR)/deneb_timescheme_ssprk.o\
          $(D_OBJ_DIR)/deneb_timescheme_steady.o\

OBJ=$(AVOCADO_OBJ)\
    $(DENEB_OBJ)\
    $(OBJ_DIR)/main.o\

TARGET=Deneb
LINKFLAGS=-L${PARMETIS_LIB} -L${METIS_LIB} -L${IDEA_LIB} -L${PETSC_LIB} -lparmetis -lmetis -lpthread -lm -ldl -lidea -mkl=sequential -lpetsc
IFLAGS=-I./$(A_INC_DIR) -I./$(D_INC_DIR) -I${PARMETIS_INC} -I${METIS_INC} -I${PETSC_INC} -I${PETSC_CONF_INC} -I${IDEA_INC}

.PHONY : all clean

all : folder $(TARGET)

folder :
	-mkdir $(OBJ_DIR)
	-mkdir $(A_OBJ_DIR)
	-mkdir $(D_OBJ_DIR)
	-mkdir $(BIN_DIR)

$(TARGET) : $(OBJ)
	$(CXX) -o $(BIN_DIR)/$(TARGET) $(OBJ) $(LINKFLAGS)

$(A_OBJ_DIR)/%.o : $(A_SRC_DIR)/%.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS) $(IFLAGS)

$(D_OBJ_DIR)/%.o : $(D_SRC_DIR)/%.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS) $(IFLAGS)

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS) $(IFLAGS)

clean :
	-rm $(A_OBJ_DIR)/*.o
	-rm $(D_OBJ_DIR)/*.o
	-rm $(OBJ_DIR)/*.o
	-rmdir $(A_OBJ_DIR)
	-rmdir $(D_OBJ_DIR)
	-rmdir $(OBJ_DIR)
	-rm $(BIN_DIR)/$(TARGET)
	-rmdir $(BIN_DIR)
