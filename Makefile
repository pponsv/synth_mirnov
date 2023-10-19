.PHONY : clean compile_lib compile_py rebuild

LIBNAME = synth_mirnov_lib
# OPT_FLAG = -Wall -Og -fbacktrace 
OPT_FLAG = -O2
SRC_LIB = ./src_lib
SRC_F2PY = ./src
OBJ_DIR = ./bld
LIB_DIR = ./lib
F90_LIB = $(wildcard $(SRC_LIB)/*.f90)
OFILES_LIB = $(patsubst $(SRC_LIB)/%.f90,$(OBJ_DIR)/%.o,$(F90_LIB))

OFILE_FLAGS = -fPIC $(OPT_FLAG) -J$(OBJ_DIR)/ -I$(OBJ_DIR)/ -ffree-form -fopenmp -I/usr/include -std=f2008 -ffree-line-length-none -Wall
F2PY_FLAGS = -fPIC $(OPT_FLAG) -ffree-form -fopenmp -L./lib/ -I$(OBJ_DIR)/  -L/usr/lib -flto -std=f2008 -ffree-line-length-none

rebuild : clean compile_lib compile_py 

compile_py : $(LIB_DIR)/lib$(LIBNAME).a $(SRC_F2PY)/synthetic_mirnov.f90
	f2py -c -m synth_mirnov --debug-capi --f90flags='$(F2PY_FLAGS)' -lgomp -lfftw3 -lm $(SRC_F2PY)/synthetic_mirnov.f90 $(LIB_DIR)/libsynth_mirnov_lib.a

compile_lib : $(LIB_DIR)/lib$(LIBNAME).a

$(LIB_DIR)/lib$(LIBNAME).a :$(OFILES_LIB) $(LIB_DIR) 
	ar rcv $(LIB_DIR)/lib$(LIBNAME).a $(OFILES_LIB)

$(OBJ_DIR)/%.o: $(SRC_LIB)/%.f90 $(OBJ_DIR)
	gfortran $(OFILE_FLAGS) -c -o $@ $<

$(OBJ_DIR)/types.o : $(OBJ_DIR)/constants.o
$(OBJ_DIR)/global.o : $(OBJ_DIR)/types.o $(OBJ_DIR)/constants.o
$(OBJ_DIR)/helper.o : $(OBJ_DIR)/constants.o
$(OBJ_DIR)/potential.o : $(OBJ_DIR)/constants.o $(OBJ_DIR)/helper.o $(OBJ_DIR)/global.o
$(OBJ_DIR)/fft_mod.o : $(OBJ_DIR)/global.o
# gfortran $(OFILE_FLAGS) -c -o $@ $(SRC_LIB)/fft_mod.f90 
$(OBJ_DIR)/derivatives.o : $(OBJ_DIR)/fft_mod.o $(OBJ_DIR)/helper.o $(OBJ_DIR)/types.o 
$(OBJ_DIR)/main.o : $(OBJ_DIR)/global.o $(OBJ_DIR)/potential.o $(OBJ_DIR)/derivatives.o $(OBJ_DIR)/types.o 

$(LIB_DIR) : 
	mkdir -p $(LIB_DIR)

$(OBJ_DIR) : 
	mkdir -p $(OBJ_DIR)

clean:
	rm -rf *.mod *.o $(OBJ_DIR) $(LIB_DIR) ./bin/* ./*.cpython* ./__*__/