project('synth_mirnov', 'c',
  version : '0.1',
  license: 'BSD-3',
  meson_version: '>=1.4.0',
  default_options : [
    'warning_level=2', 
    'default_library=static', 
    'optimization=2', 
    'fortran_std=f2008',
  ],
)
MODULE_NAME = 'synth_mirnov'

add_languages('fortran')
fortran_compiler = meson.get_compiler('fortran')

if fortran_compiler.get_id() == 'gcc'
  add_project_arguments('-ffree-line-length-none', language: 'fortran')
endif
 
py_mod = import('python')
py = py_mod.find_installation(pure: true)
py_dep = py.dependency()

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

fftw_dep = dependency('fftw3')
openmp = dependency('openmp')

lib_synth_mirnov_sources = files(['src_lib/main.f90',
'src_lib/constants.f90', 
'src_lib/derivatives.f90', 
'src_lib/fft_mod.f90', 
'src_lib/global.f90', 
'src_lib/helper.f90', 
'src_lib/initialization.f90',
'src_lib/potential.f90'])

lib_synth_mirnov = library(MODULE_NAME,
  lib_synth_mirnov_sources,
  include_directories: '/usr/include/',
  dependencies : [fftw_dep, openmp],
  install : false,
  pic : true
)

synth_mirnov_custom_target = custom_target('synth_mirnov.c',
  input : ['src/synthetic_mirnov.f90'], 
  output : [MODULE_NAME+'module.c', 
#   MODULE_NAME+'-f2pywrappers.f', 
  MODULE_NAME+'-f2pywrappers2.f90'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', MODULE_NAME, '--lower']
)

inc_np = include_directories(incdir_numpy, incdir_f2py)

py.extension_module(MODULE_NAME,
  ['src/synthetic_mirnov.f90', synth_mirnov_custom_target],
  incdir_f2py / 'fortranobject.c', # f2py file
  include_directories: inc_np,
  dependencies : [py_dep, openmp],
  link_with : lib_synth_mirnov,
  install : true,
  install_dir : meson.source_root() + '/bin/',
)
