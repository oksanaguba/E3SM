if (compile_threaded)
  #string(APPEND CFLAGS " -fopenmp")
  string(APPEND FFLAGS " -fopenmp")
  string(APPEND CXXFLAGS " -fopenmp")
  string(APPEND LDFLAGS " -fopenmp")
endif()

string(APPEND SLIBS " -L$ENV{PNETCDF_PATH}/lib -lpnetcdf")
set(NETCDF_PATH "$ENV{NETCDF_DIR}")
set(PNETCDF_PATH "$ENV{PNETCDF_DIR}")
set(PIO_FILESYSTEM_HINTS "gpfs")
string(APPEND CXX_LIBS " -lstdc++")

SET(CMAKE_C_COMPILER "cc" CACHE STRING "")
SET(CMAKE_Fortran_COMPILER "ftn" CACHE STRING "")
SET(CMAKE_CXX_COMPILER "hipcc" CACHE STRING "")

string(APPEND LDFLAGS " -L$ENV{ROCM_PATH}/lib -lamdhip64")
string(APPEND CXXFLAGS " -I$ENV{ROCM_PATH}/include")

# Crusher: this resolves a crash in mct in docn init
if (NOT DEBUG)
  string(APPEND CFLAGS " -O2 -hnoacc -hfp0 -hipa0")
  string(APPEND FFLAGS " -O2 -hnoacc -hfp0 -hipa0")
endif()

string(APPEND CPPDEFS " -DCPRCRAY")
