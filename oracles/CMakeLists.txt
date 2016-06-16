set(oracle_source_files
spanner.cpp
spanner.hpp
oracle.cpp 
oracle.hpp
thorup_zwick_oracle.cpp
thorup_zwick_oracle.hpp
chechik_oracle.cpp
chechik_oracle.hpp
with_dijkstra_oracle.cpp
with_dijkstra_oracle.hpp
with_thorup_zwick_oracle.cpp
with_thorup_zwick_oracle.hpp
restricted_chechik_oracle.cpp
restricted_chechik_oracle.hpp
with_chechik_oracle.cpp
with_chechik_oracle.hpp
)
add_library(oracles ${oracle_source_files})

#add_subdirectory(lib)

target_link_libraries(oracles lib)





#add_library(oracle STATIC oracle.cpp oracle.hpp)
#add_library(thorup_zwick_ora STATIC thorup_zwick_oracle.cpp thorup_zwick_oracle.hpp)
#add_library(chechik_ora STATIC chechik_oracle.cpp chechik_oracle.hpp)
#add_library(with_dijkstra_ora STATIC with_dijkstra_oracle.cpp with_dijkstra_oracle.hpp)