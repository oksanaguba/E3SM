string(APPEND CONFIG_ARGS " --enable-filesystem-hints=lustre")
string(APPEND CPPDEFS " -DLINUX")
if (DEBUG)
  string(APPEND FFLAGS " -check all -ftrapuv")
endif()
set(NETCDF_PATH "$ENV{NETCDF_HOME}")
set(PIO_FILESYSTEM_HINTS "lustre")
set(PNETCDF_PATH "$ENV{PNETCDFROOT}")
string(APPEND SLIBS " -L${NETCDF_PATH}/lib -lnetcdf -lnetcdff -L$ENV{MKL_PATH}/lib/intel64 -lmkl_rt")
if (MPILIB STREQUAL mpich2)
  string(APPEND SLIBS " -mkl=cluster")
endif()
if (MPILIB STREQUAL mpi-serial)
  string(APPEND SLIBS " -mkl")
endif()
