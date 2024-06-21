cmake_minimum_required(VERSION 3.10)
project(particle_reconstruction)
find_package(ROOT REQUIRED COMPONENTS RIO)

add_subdirectory(src)
add_executable(particle_reconstruction analysis.cpp)
target_link_libraries(particle_reconstruction ROOT::RIO ROOT::Hist ROOT::Gpad ROOT::Tree ROOT::TreePlayer ROOT::Core ROOT::Gui EventCollector)
target_include_directories(particle_reconstruction PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/src")