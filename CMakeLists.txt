cmake_minimum_required(VERSION 3.17.5)

# define build environments
set( CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local/"
  CACHE STRING "Select where to install the library." )
set(CMAKE_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}
  CACHE STRING "Select where to build the library." )
set(MODULE_DIR ${CMAKE_BUILD_DIR}/mod)

if (DEFINED SKBUILD_PROJECT_NAME)
	set(SKBUILD_PROJECT_NAME ${SKBUILD_PROJECT_NAME})
else()
	set(SKBUILD_PROJECT_NAME "artemis")
endif()


# set compiler
set(CMAKE_Fortran_COMPILER gfortran
  CACHE STRING "Select Fortran compiler." )  # Change this to your desired compiler
set(CMAKE_C_COMPILER gcc
  CACHE STRING "Select C compiler." )  # Change this to your desired compiler
set(CMAKE_Fortran_STANDARD 2018)

# set the project version
file(READ "fpm.toml" ver)
string(REGEX MATCH "version = \"([0-9]+.[0-9]+.[0-9]+)\"" _ ${ver})
set(PROJECT_VERSION ${CMAKE_MATCH_1})
message(STATUS "Project version: ${PROJECT_VERSION}")

# set the project name
project(artemis
    VERSION ${PROJECT_VERSION}
    LANGUAGES C Fortran
)

# set the library name
set( LIB_NAME ${PROJECT_NAME} )
set( PROJECT_DESCRIPTION
  "Fortran materials lattice matcher" )
set( PROJECT_URL "https://github.com/ExeQuantCode/artemis" )
set( CMAKE_CONFIGURATION_TYPES "Release" "Dev" "Debug"
  CACHE STRING "List of configurations types." )
set( CMAKE_BUILD_TYPE "Release"
  CACHE STRING "Select which configuration to build." )

# set options
option(BUILD_EXECUTABLE "Build the Fortran executable" On)

# Define the sources
set(SRC_DIR src)
set(FORTRAN_SRC_DIR ${SRC_DIR}/fortran)
set(LIB_DIR ${FORTRAN_SRC_DIR}/lib)

# Library source files
set(LIB_FILES
    mod_constants.f90
    mod_misc.f90
    mod_io_utils.F90
    mod_help.f90
    mod_misc_maths.f90
    mod_misc_linalg.f90
    mod_geom_utils.f90
    mod_io_utils_extd.F90
    mod_sym.f90
    mod_terminations.f90
    mod_misc_types.f90
    mod_intf_identifier.f90
    mod_plane_matching.f90
    mod_lat_compare.f90
    mod_swapping.f90
    mod_shifting.f90
)

# Main source files
set(SPECIAL_LIB_FILES
    mod_geom_rw.f90
    mod_term_generator.f90
    mod_intf_generator.f90
)


foreach(lib ${LIB_FILES})
    list(APPEND PREPENDED_LIB_FILES ${LIB_DIR}/${lib})
endforeach()
foreach(lib ${SPECIAL_LIB_FILES})
    list(APPEND PREPENDED_LIB_FILES ${LIB_DIR}/${lib})
endforeach()
message(STATUS "Modified LIB_FILES: ${PREPENDED_LIB_FILES}")



set(SRC_FILES
    artemis.f90
)
foreach(lib ${SPECIAL_LIB_FILES})
    list(APPEND F90WRAP_FORTRAN_SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/${LIB_DIR}/${lib})
endforeach()
foreach(src ${SRC_FILES})
    list(APPEND F90WRAP_FORTRAN_SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/${FORTRAN_SRC_DIR}/${src})
    list(APPEND PREPENDED_SRC_FILES ${FORTRAN_SRC_DIR}/${src})
endforeach()


set(EXECUTABLE_FILES
    mod_tools_infile.f90
    default_infile.f90
    inputs.f90
    aspect.f90
    main.f90
)
set(APP_DIR app)
foreach(src ${EXECUTABLE_FILES})
    list(APPEND PREPENDED_EXECUTABLE_FILES ${APP_DIR}/${src})
endforeach()


# initialise flags
set(CPPFLAGS "")
set(CFLAGS "")
set(MODULEFLAGS "")
set(MPFLAGS "")
set(WARNFLAGS "")
set(DEVFLAGS "")
set(DEBUGFLAGS "")
set(MEMFLAGS "")
set(OPTIMFLAGS "")
set(FASTFLAGS "")

# set flags based on compiler
if (CMAKE_Fortran_COMPILER MATCHES ".*gfortran.*" OR CMAKE_Fortran_COMPILER MATCHES ".*gcc.*")
    message(STATUS "Using gfortran compiler")
    set(PPFLAGS -cpp)
    set(MPFLAGS -fopenmp -lgomp -floop-parallelize-all -ftree-parallelize-loops=32)
    set(WARNFLAGS -Wall)
    set(DEVFLAGS -g -fbacktrace -fcheck=all -fbounds-check -Og)
    set(DEBUGFLAGS -fbounds-check)
    set(MEMFLAGS -mcmodel=large)
    set(OPTIMFLAGS -O3 -march=native)
    set(FASTFLAGS -Ofast -march=native)
    set(PYTHONFLAGS -c -O3 -fPIC)
elseif (CMAKE_Fortran_COMPILER MATCHES ".*nag.*")
    message(STATUS "Using nag compiler")
    set(PPFLAGS -f2018 -fpp)
    set(MPFLAGS -openmp)
    set(WARNFLAGS -Wall)
    set(DEVFLAGS -g -mtrace -C=all -colour -O0)
    set(DEBUGFLAGS -C=array)
    set(MEMFLAGS -mcmodel=large)
    set(OPTIMFLAGS -O3)
    set(FASTFLAGS -Ofast)
elseif (CMAKE_Fortran_COMPILER MATCHES ".*ifort.*" OR CMAKE_Fortran_COMPILER MATCHES ".*ifx.*")
    message(STATUS "Using intel compiler")
    set(PPFLAGS -fpp)
    set(MPFLAGS -qopenmp)
    set(WARNFLAGS -warn all)
    set(DEVFLAGS -check all -warn)
    set(DEBUGFLAGS -check all -fpe0 -warn -tracekback -debug extended)
    set(MEMFLAGS -mcmodel=large)
    set(OPTIMFLAGS -O3)
    set(FASTFLAGS -Ofast)
else()
    # Code for other Fortran compilers
    message(STATUS "Using a different Fortran compiler")
endif()



set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${PPFLAGS}")


# create the library
add_library(${PROJECT_NAME} STATIC ${PREPENDED_LIB_FILES} ${PREPENDED_SRC_FILES})
set_target_properties(${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${MODULE_DIR})
target_link_libraries(${PROJECT_NAME} PUBLIC)

# replace ".f90" with ".mod"
string(REGEX REPLACE "\\.[^.]*$" ".mod" MODULE_FILES "${SRC_FILES}")

set(ETC_MODULE_FILES "")
# Loop through each Fortran file
foreach(FILE ${PREPENDED_LIB_FILES})
    # Read the content of the Fortran file
    file(READ "${FILE}" FILE_CONTENTS)

    # Use a regular expression to extract the module name
    string(REGEX MATCH "^module[ \t]+([a-zA-Z0-9_]+)" MATCH "${FILE_CONTENTS}")

    # If a match is found, extract the module name (the first capture group)
    if(MATCH)
        string(REGEX REPLACE "module[ \t]+([a-zA-Z0-9_]+)" "\\1" MODULE_NAME "${MATCH}")

        # Append the module name with .mod to the list
        list(APPEND ETC_MODULE_FILES "${MODULE_DIR}/${MODULE_NAME}.mod")
    endif()
endforeach()

# installation
install(FILES ${MODULE_DIR}/${MODULE_FILES} DESTINATION ${SKBUILD_PROJECT_NAME}/include)
install(FILES ${ETC_MODULE_FILES} DESTINATION ${SKBUILD_PROJECT_NAME}/etc)
install(TARGETS ${PROJECT_NAME} DESTINATION ${SKBUILD_PROJECT_NAME}/lib)
set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${PROJECT_VERSION})

# set compile options based on different build configurations
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Release>:${OPTIMFLAGS}>")
# target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Release>:${MPFLAGS}>")
# target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Release>:${PYTHONFLAGS}>")

target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Dev>:${DEVFLAGS}>")

target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Debug>:${DEBUGFLAGS}>")
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Debug>:${WARNFLAGS}>")
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Debug>:${MPFLAGS}>")
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Debug>:${PYTHONFLAGS}>")

if (BUILD_EXECUTABLE)
    add_executable(${PROJECT_NAME}_executable ${PREPENDED_EXECUTABLE_FILES})
    target_link_libraries(${PROJECT_NAME}_executable PRIVATE ${PROJECT_NAME})
    install(TARGETS ${PROJECT_NAME}_executable DESTINATION ${SKBUILD_PROJECT_NAME}/bin)
    set_target_properties(${PROJECT_NAME}_executable PROPERTIES Fortran_MODULE_DIRECTORY ${MODULE_DIR})

    target_compile_options(${PROJECT_NAME}_executable PUBLIC "$<$<CONFIG:Release>:${OPTIMFLAGS}>")
    # target_compile_options(${PROJECT_NAME}_executable PUBLIC "$<$<CONFIG:Release>:${MPFLAGS}>")
    # target_compile_options(${PROJECT_NAME}_executable PUBLIC "$<$<CONFIG:Release>:${PYTHONFLAGS}>")

    target_compile_options(${PROJECT_NAME}_executable PUBLIC "$<$<CONFIG:Dev>:${DEVFLAGS}>")

    target_compile_options(${PROJECT_NAME}_executable PUBLIC "$<$<CONFIG:Debug>:${DEBUGFLAGS}>")
    target_compile_options(${PROJECT_NAME}_executable PUBLIC "$<$<CONFIG:Debug>:${WARNFLAGS}>")
    target_compile_options(${PROJECT_NAME}_executable PUBLIC "$<$<CONFIG:Debug>:${MPFLAGS}>")
    target_compile_options(${PROJECT_NAME}_executable PUBLIC "$<$<CONFIG:Debug>:${PYTHONFLAGS}>")

    set_target_properties(${PROJECT_NAME}_executable PROPERTIES VERSION ${PROJECT_VERSION})
endif()
