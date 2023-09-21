LIBNAME = synth_mirnov_lib

SRC_LIB = ./src_lib
OBJ_DIR = ./bld
F90_LIB = $(wildcard $(SRC_LIB)/*.f90)
OFILES_LIB = $(patsubst $(SRC_LIB)/%.f90,$(OBJ_DIR)/%.o,$(F90_LIB))

OFILE_FLAGS = -fPIC -O2 -J./bld/ -I./bld/ -ffree-form -fopenmp -I/usr/include -flto
F2PY_FLAGS = -fPIC -O2 -ffree-form -fopenmp -L./lib/ -I./bld/  -L/usr/lib -flto

run_py : compile_py
	python3 synth_mirnov.py

compile_py : clean compile_lib ./src/synthetic_mirnov.f90
	f2py -c -m synth_mirnov  --f90flags='$(F2PY_FLAGS)' -lgomp -lfftw3 ./src/synthetic_mirnov.f90 ./lib/libsynth_mirnov_lib.a

compile_lib : ./lib/lib$(LIBNAME).a

test : ./lib/lib$(LIBNAME).a ./src/test.f90
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
$(OBJ_DIR)/fft.o : $(OBJ_DIR)/global.o
	gfortran $(OFILE_FLAGS) -c -o $@ $(SRC_LIB)/fft.f90 
$(OBJ_DIR)/derivatives.o : $(OBJ_DIR)/fft.o $(OBJ_DIR)/helper.o $(OBJ_DIR)/types.o 


clean:
	rm -rf *.mod *.o ./bld/* ./lib/* ./bin/* *.cpython*