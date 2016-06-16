#add_library(structs STATIC structs.cpp structs.hpp)
#add_library(std_methods STATIC std_methods.cpp std_methods.hpp)
#add_library(sssp STATIC sssp.cpp sssp.hpp)
#add_library(parser STATIC parser.cpp parser.hpp)



set(lib_source_files
structs.cpp 
structs.hpp
std_methods.cpp
std_methods.hpp
sssp.hpp
sssp.cpp
tst_methods.hpp
tst_methods.cpp
parser.cpp
parser.hpp
)
add_library(lib ${lib_source_files})

target_link_libraries(lib)