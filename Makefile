LIBNAME = synth_mirnov_lib
# OPT_FLAG = -fcheck=all -Og -fbacktrace 
OPT_FLAG = -O2
SRC_LIB = ./src_lib
OBJ_DIR = ./bld
F90_LIB = $(wildcard $(SRC_LIB)/*.f90)
OFILES_LIB = $(patsubst $(SRC_LIB)/%.f90,$(OBJ_DIR)/%.o,$(F90_LIB))

OFILE_FLAGS = -fPIC $(OPT_FLAG) -J./bld/ -I./bld/ -ffree-form -fopenmp -I/usr/include
F2PY_FLAGS = -fPIC $(OPT_FLAG) -ffree-form -fopenmp -L./lib/ -I./bld/  -L/usr/lib -flto

run_py : compile_py
	python3 synth_mirnov.py

compile_py : clean compile_lib ./src/synthetic_mirnov.f90
	f2py -c -m synth_mirnov  --f90flags='$(F2PY_FLAGS)' -lgomp -lfftw3 ./src/synthetic_mirnov.f90 ./lib/libsynth_mirnov_lib.a

compile_lib : ./lib/lib$(LIBNAME).a

test : ./lib/lib$(LIBNAME).a ./src/test.f90
	gfortran $(OPT_FLAG) -J./bld/ -flto -fPIC -I/usr/include -o ./bin/xtest ./src/test.f90 -fopenmp
	./bin/xtest

libtest : ./lib/lib$(LIBNAME).a ./src/test.f90
	gfortran $(OFILE_FLAGS) -o ./bin/xtest ./src/test.f90  -L./lib/ -l$(LIBNAME)
	./bin/xtest

./lib/lib$(LIBNAME).a : $(OFILES_LIB)
	ar rcv ./lib/lib$(LIBNAME).a $(OFILES_LIB)

$(OBJ_DIR)/%.o: $(SRC_LIB)/%.f90
	gfortran $(OFILE_FLAGS) -c -o $@ $<

$(OBJ_DIR)/types.o : $(OBJ_DIR)/constants.o
$(OBJ_DIR)/global.o : $(OBJ_DIR)/types.o $(OBJ_DIR)/constants.o
$(OBJ_DIR)/helper.o : $(OBJ_DIR)/constants.o
$(OBJ_DIR)/potential.o : $(OBJ_DIR)/constants.o $(OBJ_DIR)/helper.o $(OBJ_DIR)/global.o
$(OBJ_DIR)/fft_mod.o : $(OBJ_DIR)/global.o
	gfortran $(OFILE_FLAGS) -c -o $@ $(SRC_LIB)/fft_mod.f90 
$(OBJ_DIR)/derivatives.o : $(OBJ_DIR)/fft_mod.o $(OBJ_DIR)/helper.o $(OBJ_DIR)/types.o 
$(OBJ_DIR)/main.o : $(OBJ_DIR)/global.o $(OBJ_DIR)/potential.o $(OBJ_DIR)/derivatives.o $(OBJ_DIR)/types.o 

clean:
	rm -rf *.mod *.o ./bld/* ./lib/* ./bin/* *.cpython*