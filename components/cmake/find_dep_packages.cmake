# This file is for finding pacakges needed by E3SM. It should be included
# from the main CMakeLists.txt file.
#
# Finding the correct packages will likely depend on the ${Package}_ROOT
# environment variable being set by config_machines.xml for your machine. Note
# that is environment var is case sensitive.

# Machine env vars should follow the following pattern:
#   <env name="Package_ROOT">$SHELL{if [ -z "$Package_ROOT" ]; then echo /default/install/location; else echo "$Package_ROOT"; fi}</env>
#
# This will allow users to easily specify a different location for all their cases by
# simply setting ${Package}_ROOT in their shell.

# Kokkos' find_package needs to come before other find_packages
# that may define Kokkos targets so we can avoid duplicate target
# errors.
if (USE_KOKKOS)

  # Kokkos will be built in the sharedlibs if Kokkos_ROOT is
  # unset.
  if (NOT DEFINED ENV{Kokkos_ROOT})
    set(ENV{Kokkos_ROOT} ${INSTALL_SHAREDPATH})
  endif()

  find_package(Kokkos REQUIRED)
endif()

# Albany depends on Trilinos
if (USE_ALBANY OR USE_TRILINOS)
  find_package(Trilinos REQUIRED)
endif()

if (USE_ALBANY)
  find_package(Albany REQUIRED)
endif()

if (USE_PETSC)
  find_package(PETSc REQUIRED)
endif()

# Placeholder code until spio gets a proper config.cmake
if (NOT PIO_LIBDIR)
  set(PIO_LIBDIR "${INSTALL_SHAREDPATH}/lib")
endif()
if (PIO_VERSION STREQUAL 2)
  # This is a pio2 library
  set(PIOLIBS "${PIO_LIBDIR}/libpiof.a;${PIO_LIBDIR}/libpioc.a")
else()
  # This is a pio1 library
  set(PIOLIBS "${PIO_LIBDIR}/libpio.a")
endif()

find_package(MCT REQUIRED)
find_package(CsmShare REQUIRED)
