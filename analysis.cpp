#include "EventCollector.h"
#include <iostream>

int main(int argc, char** argv)
{
    EventCollector evc;

    evc.initialize_events();
    evc.filter_initial_events();
//    evc.init_masses_and_energy(0.13957039);
    evc.init_masses_and_energy(0.493667);
    evc.analyze("hist0.pdf");
    evc.filter();
    evc.analyze("hist0.pdf");
    evc.reconstruct_particles();
    evc.analyze_reco("hist1.pdf");

    return 0;
}