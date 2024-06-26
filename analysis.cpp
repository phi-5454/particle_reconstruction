#include "EventCollector.h"
#include <iostream>

int main(int argc, char** argv)
{
    EventCollector evc;

    evc.initialize_events();
    evc.filter_initial_events();
    evc.init_masses_and_energy(0.13957039);
    evc.analyze("hist0.pdf");
    evc.filter();
    std::cout << "Ended filtering, started reconstructing" << std::endl;
    evc.reconstruct_particles();
    std::cout << "Ended reconstructing, started analyzing" << std::endl;
    evc.analyze_reco("hist1.pdf");
    std::cout << "Finished analyzing" << std::endl;

    return 0;
}