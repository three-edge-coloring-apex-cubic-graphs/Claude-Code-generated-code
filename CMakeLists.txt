cmake_minimum_required(VERSION 3.20)
project(apex_reproducibility CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

add_library(apex_core
  src/core/graph.cpp
  src/core/rotations.cpp
  src/core/walks.cpp
  src/io/parse.cpp
  src/hom/hom.cpp
  src/hom/blocked.cpp
  src/kbar/kbar.cpp
  src/rules/combine.cpp
  src/cartwheel/cartwheel.cpp
  src/cartwheel/charge.cpp
  src/cartwheel/verify.cpp
  src/island/outer_extension.cpp
  src/island/planarity.cpp
  src/island/make_outer_extension.cpp
  src/island/island.cpp
  src/island/hom_images.cpp
  src/reduce/kempe.cpp
  src/reduce/island_graph.cpp
  src/reduce/coloring.cpp
  src/reduce/deletable.cpp
  src/reduce/semi_reducible.cpp
)
target_include_directories(apex_core PUBLIC src)
find_package(Threads REQUIRED)
target_link_libraries(apex_core PUBLIC Threads::Threads)

add_executable(verify_b1 src/bin/verify_b1.cpp)
target_link_libraries(verify_b1 apex_core)

add_executable(verify_b2 src/bin/verify_b2.cpp)
target_link_libraries(verify_b2 apex_core)

add_executable(verify_b3 src/bin/verify_b3.cpp)
target_link_libraries(verify_b3 apex_core)

add_executable(verify_b3_reduce src/bin/verify_b3_reduce.cpp)
target_link_libraries(verify_b3_reduce apex_core)

enable_testing()

add_executable(test_parse tests/test_parse.cpp)
target_link_libraries(test_parse apex_core)
add_test(NAME parse COMMAND test_parse WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
add_test(NAME lemma_b1 COMMAND verify_b1 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
