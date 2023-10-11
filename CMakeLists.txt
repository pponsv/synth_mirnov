cmake_minimum_required(VERSION 3.18)

enable_language(Fortran)

project(synth_mirnov)

find_package(OpenMP)

INCLUDE_DIRECTORIES(/usr/include)
LINK_DIRECTORIES(/usr/lib/x86_64-linux-gnu)

add_library(fftw3 STATIC IMPORTED)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

# synth_mirnov_lib
file(GLOB lib_src "./src_lib/*.f90")
add_library(synth_mirnov_lib STATIC ${lib_src})
target_link_libraries(synth_mirnov_lib PUBLIC fftw3)

if(OpenMP_Fortran_FOUND)
    target_link_libraries(synth_mirnov_lib PUBLIC OpenMP::OpenMP_Fortran)
endif()

# F2PY
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)
execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
    COMMAND "${Python_EXECUTABLE}"
    -c "import numpy.f2py; print(numpy.f2py.get_include())"
    OUTPUT_VARIABLE F2PY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

set(module_name "synth_mirnov")
set(module_src "${CMAKE_SOURCE_DIR}/src/synthetic_mirnov.f90")

# # Add the f2py module target
# add_library(<your_f2py_module_target_name> MODULE ${f2py_module_source_file})

# # Link the f2py module target with NumPy
# target_link_libraries(<your_f2py_module_target_name> PRIVATE Python::NumPy)

# # Set the include directories for the f2py module target
# target_include_directories(<your_f2py_module_target_name> PRIVATE ${NUMPY_INCLUDE_DIR})