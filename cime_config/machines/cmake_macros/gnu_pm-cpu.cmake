string(APPEND CONFIG_ARGS " --host=cray")
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DHAVE_GETTIMEOFDAY")
endif()
string(APPEND CMAKE_C_FLAGS_RELEASE " -O2 -g")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -O2 -g")
set(MPICC "cc")
set(MPICXX "CC")
set(MPIFC "ftn")
set(SCC "gcc")
set(SCXX "g++")
set(SFC "gfortran")
set(ALBANY_PATH "/global/common/software/fanssie/albany-e3sm-serial-release-gcc")
