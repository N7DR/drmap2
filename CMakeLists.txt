cmake_minimum_required(VERSION 3.0)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

project(drmap2)

#set (CMAKE_CXX_COMPILER "g++")
#set (CMAKE_CXX_COMPILER "/usr/local/programs/gcc-14.2.0/bin/g++-14")
set (CMAKE_CXX_COMPILER "/usr/local/programs/gcc-15.1.0/bin/g++-15")

set ( SOURCE_FILES src/command_line.cpp
                   src/diskfile.cpp
                   src/drmap.cpp
                   src/geotiff.cpp
                   src/r_figure.cpp
                   src/string_functions.cpp
    )

set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g3 -Wall -Wno-reorder -pipe -std=c++26" )

add_executable(drmap2 ${SOURCE_FILES})

target_include_directories (drmap2 PUBLIC include /usr/share/R/include /usr/lib/R/site-library/Rcpp/include /usr/lib/R/site-library/RInside/include)

#set (CMAKE_SHARED_LINKER_FLAGS "-Wl,--export-dynamic -fopenmp -Wl,-rpath,/usr/lib/R/site-library/RInside/lib")
set (CMAKE_SHARED_LINKER_FLAGS "-Wl,--export-dynamic -fopenmp")

#link_directories (/lib /usr/lib/R/lib /usr/lib/R/site-library/RInside/lib)
#target_link_directories (drmap2 PUBLIC /usr/lib/R/site-library/RInside/lib)
target_link_directories (drmap2 PUBLIC /usr/local/programs/gcc-15.1.0/lib64/ /usr/lib/R/site-library/RInside/lib)
target_link_libraries (drmap2 gdal R RInside)
#add_link_options (-Wl,--export-dynamic -fopenmp -Wl,-rpath,/usr/lib/R/site-library/RInside/lib)
#add_link_options (LINKER: --export-dynamic -fopenmp LINKER: -rpath,/usr/lib/R/site-library/RInside/lib)


set( CMAKE_VERBOSE_MAKEFILE on )

#install(TARGETS wae RUNTIME DESTINATION bin)
