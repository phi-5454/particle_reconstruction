#include "EventCollector.h"
#include <iostream>

int main()
{
    EventCollector evc("/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/ntuples/data/TOTEM2*.root?#tree",
                       "/afs/cern.ch/user/p/ptuomola/private/particle_reconstruction_results.root");
//                     "/eos/user/y/yelberke/TOTEM_2018_ADDEDVARS_OUT/"

    evc.initialize_events(false);

    std::cout << "Filtering events" << std::endl;

    auto l1 = [](Event *e) { return e->ntracks == 4; };
    evc.filter_events(l1);
    evc.filter_events([](Event *e) {
        int j = 0;
        for (int i = 0; i < 4; ++i) {
        j += e->get_particle(0, 0, i)->q;
        }
        return j == 0;
    });

    // Primary vertex Z position
    evc.filter_events_distribution(
        [](Event *event) {
            std::vector<double> values = {event->zPV};
            return values;
        },
        "gaus", 3);

    // Particle smallest distance from the primary vertex in xy-plane
    evc.filter_events_distribution(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
            values[i] = event->get_particle(0, 0, i)->dxy;
            return values;
        },
        "gaus", 3, 200, -2, 2, "Title");

    // Particle smallest distance from the primary vertex in z-axis
    evc.filter_events_distribution(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
            values[i] = event->get_particle(0, 0, i)->dz;
            return values;
        },
        "gaus", 3, 200, -3, 3, "Title");

    std::cout << "Finished filtering" << std::endl;

    evc.init_masses_and_energy(0.13957039);

//    evc.init_masses_and_energy(0.493667);
//    evc.analyze("hist11");
//    evc.analyze_new("hist21");
//    evc.filter();
//    evc.analyze("hist12");
//    evc.analyze_new("hist22");

    std::cout << "Reconstructing particles" << std::endl;

    evc.reconstruct_particles();

    std::cout << "Finished reconstructing" << std::endl;

//    evc.filter2();
//    evc.reconstruct_particles();

    evc.analyze_reco("hist3");

    return 0;
}