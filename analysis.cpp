#include "EventCollector.h"
#include <iostream>

int main(int argc, char** argv)
{
    EventCollector evc;

    evc.initialize_events(false);
    evc.filter_initial_events();
    evc.init_masses_and_energy(0.13957039);
//    evc.init_masses_and_energy(0.493667);
//    evc.analyze("hist11");
//    evc.analyze_new("hist21");
//    evc.filter();
//    evc.analyze("hist12");
//    evc.analyze_new("hist22");
    evc.reconstruct_particles();
//   evc.filter2();
//    evc.reconstruct_particles();
    evc.analyze_reco("hist3");

    return 0;
}