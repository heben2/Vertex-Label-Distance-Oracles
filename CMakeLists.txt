cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#SET(CMAKE_CXX_FLAGS "-std=c++0x") #needed for cpp11 and cpp14

#SET(CMAKE_CXX_FLAGS "-std=c++14" CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "" FORCE) #needed for cpp11 and cpp14

SET(CMAKE_CXX_FLAGS "-std=c++14") #needed for cpp11 and cpp14

project(distOracles)
#SET(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "" FORCE)




#include_directories(${CMAKE_CURRENT_SOURCE_DIR})
#add_subdirectory(lib) #only used if sub dir contains cmakelists-file
#include_directories(lib)

SET(BASEPATH "${CMAKE_SOURCE_DIR}")
INCLUDE_DIRECTORIES("${BASEPATH}")
add_subdirectory(oracles/lib)
add_subdirectory(oracles)


add_executable(main main.cpp)
target_link_libraries(main lib oracles)







#add_library(structs STATIC lib/structs.cpp lib/structs.hpp)
#add_library(std_methods STATIC lib/std_methods.cpp lib/std_methods.hpp)
#add_library(sssp STATIC lib/sssp.cpp lib/sssp.hpp)
#add_library(parser STATIC lib/parser.cpp lib/parser.hpp)

#add_library(oracle STATIC oracle.cpp oracle.hpp)
#add_library(thorup_zwick_ora STATIC thorup_zwick_oracle.cpp thorup_zwick_oracle.hpp)
#add_library(chechik_ora STATIC chechik_oracle.cpp chechik_oracle.hpp)
#add_library(with_dijkstra_ora STATIC with_dijkstra_oracle.cpp with_dijkstra_oracle.hpp)

#target_link_libraries(thorup_zwick_ora oracle sssp)
#target_link_libraries(chechik_ora oracle sssp)
#target_link_libraries(with_dijkstra_ora oracle std_methods structs sssp)

#add_executable(main main.cpp)
#target_link_libraries(main structs sssp thorup_zwick_ora parser chechik_ora with_dijkstra_ora oracle std_methods)


########## documentation and examples:
#http://stackoverflow.com/questions/10851247/how-to-activate-c-11-in-cmake
#http://stackoverflow.com/questions/6352123/multiple-directories-under-cmake
