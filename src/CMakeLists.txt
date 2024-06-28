find_package(ROOT REQUIRED COMPONENTS RIO)
add_library(EventCollector EventCollector.cpp Event.cpp Particle.cpp Proton.cpp)

target_link_libraries(EventCollector ROOT::RIO ROOT::Hist ROOT::Gpad ROOT::Tree ROOT::Core ROOT::Gui)