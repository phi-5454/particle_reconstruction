#include "EventCollector.h"
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>

void initialize(EventCollector& evc, std::string part) {
    std::cout << "Initializing events." << std::endl;

    evc.initialize_events(false);

    if (part == "pion") evc.init_masses_and_energy(0.13957039);
    else evc.init_masses_and_energy(0.493667);

    std::cout << "Finished initializing events." << std::endl;
}

void filter(EventCollector& evc) {
    std::cout << "Filtering events." << std::endl;

    // Four-track events
    evc.filter_events([](Event *e) { return e->ntracks == 4; });

    // Net zero charge events
    evc.filter_events([](Event *e) {
        int j = 0;
        for (int i = 0; i < 4; ++i) {
        j += e->get_particle(0, 0, i)->q;
        }
        return j == 0;
    });

    // Primary vertex X position
    evc.filter_events_distribution(
        [](Event *event) {
            std::vector<double> values = {event->xPV};
            return values;
        },
        "gaus", 3);

    // Primary vertex Y position
    evc.filter_events_distribution(
        [](Event *event) {
            std::vector<double> values = {event->yPV};
            return values;
        },
        "gaus", 3);

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

    std::cout << "Finished filtering events." << std::endl;
}

void reconstruct(EventCollector& evc) {
    std::cout << "Reconstructing particles." << std::endl;
    evc.reconstruct_particles();
    std::cout << "Finished reconstructing." << std::endl;
}

void analyze_data(EventCollector& evc, std::string filename) {
    std::cout << "Analyzing data." << std::endl;
    TFile *results = TFile::Open(evc.results.c_str(), "");

    TCanvas *c1 = new TCanvas("c1", "c1");
    c1->Draw();

    c1->cd(1);
    TH1 *h11 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values = {event->zPV};
            return values;
        },
        200, -15, 15, "Primary vertex Z position", true);

    c1->SaveAs((filename + "_data.pdf").c_str());

    std::cout << "Finished analyzing data." << std::endl;
}

void analyze_reco1(EventCollector& evc, std::string filename) {
    std::cout << "Analyzing the first iteration of recreated particles." << std::endl;
    TFile *results = TFile::Open(evc.results.c_str(), "");

    TCanvas *c2 = new TCanvas("c2", "c2");
    c2->Draw();

    c2->cd(1);
    TH1 *h21 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                values[2 * i + j] = event->get_particle(1, i, j)->mass;
            }
            }
            return values;
        },
        160, 0.2, 1, "Mass of recreated particles, assumed pions",
        true);
    
    c2->SaveAs((filename + "_reco1.pdf").c_str());

    std::cout << "Finished analyzing the first iteration of recreated particles." << std::endl;
}

void analyze_reco2(EventCollector& evc, std::string filename) {
    std::cout << "Analyzing the second iteration of recreated particles." << std::endl;
    TFile *results = TFile::Open(evc.results.c_str(), "");

    TCanvas *c3 = new TCanvas("c3", "c3");
    c3->Draw();

    c3->cd(1);
    TH1 *h21 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 2; ++i) {
                values[i] = event->get_particle(2, i, 0)->mass;
            }
            return values;
        },
        200, 1.5, 2.5, "Mass of the recreated particle",
        true);
    
    c3->SaveAs((filename + "_reco2.pdf").c_str());

    std::cout << "Finished analyzing the second iteration of recreated particles." << std::endl;
}

int main()
{
    EventCollector evc(
             "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/ntuples/data/TOTEM2*.root?#tree"
//               "/eos/user/y/yelberke/TOTEM_2018_ADDEDVARS_OUT/minimal/TOTEM*.root?#tree"
               ,"/afs/cern.ch/user/p/ptuomola/private/particle_reconstruction_results.root");

    initialize(evc, "pion");
    analyze_data(evc, "histogram1");
    filter(evc);
    reconstruct(evc);
    analyze_reco1(evc, "histogram1");
    reconstruct(evc);
    analyze_reco2(evc, "histogram1");

    return 0;
}