#include "TCanvas.h"
#include "TFile.h"
#include <iostream>
#include <math.h>
#include "EventCollector.h"

void initialize(EventCollector& evc, std::string part) {
    std::cout << "Initializing events." << std::endl;

    evc.initialize_events(true);

    if (part == "pion") evc.init_masses_and_energy(0.13957039);
    else evc.init_masses_and_energy(0.493667);

    std::cout << "Finished initializing events." << std::endl;
}

Double_t CauchyDist(Double_t *x, Double_t *par) {
    return par[2] * par[0] / (2 * (pow(par[0], 2) / 4 + pow(x[0] - par[1], 2)));
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

    // Primary vertex XY position
    evc.filter_events_distribution(
        [](Event *event) {
            std::vector<double> values = {sqrt(pow(event->xPV, 2) + pow(event->yPV, 2))};
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
        "gaus", 3);

    // Particle smallest distance from the primary vertex in z-axis
    evc.filter_events_distribution(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
            values[i] = event->get_particle(0, 0, i)->dz;
            return values;
        },
        "gaus", 3);

    // No elastic protons
    evc.filter_events([](Event * event) {
        double pt = 0;
        for (int i = 0; i < 2; ++i)
            pt += event->get_proton(i)->px + event->get_proton(i)->py;
        return pt > 0.2;
    });

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

    TCanvas *c11 = new TCanvas("c11", "c11");
    c11->DivideSquare(4);
    c11->Draw();

    c11->cd(1);
    // Primary vertex XY position
    evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values = {sqrt(pow(event->xPV, 2) + pow(event->yPV, 2))};
            return values;
        },
        100, 0.1, 0.2, "Event rPV", true);

    c11->cd(2);
    // Primary vertex Z position
    TH1 *h11 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values = {event->zPV};
            return values;
        },
        100, -15, 15, "Event zPV", true);
    h11->Fit("gaus");


    c11->cd(3);
    // Particle smallest distance from the primary vertex in z-axis
    evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
                values[i] = event->get_particle(0, 0, i)->dz;
            return values;
        },
        200, -0.3, 0.3, "Particle dz", true);

    c11->cd(4);
    // Particle dz / dz error
    evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
                values[i] = event->get_particle(0, 0, i)->dz / event->get_particle(0, 0, i)->dzErr;
            return values;
        },
        200, -5, 5, "Particle dz / dzErr", true);

    c11->SaveAs((filename + "_dataA.pdf").c_str());

    TCanvas *c12 = new TCanvas("c12", "c12");
    c12->DivideSquare(4);
    c12->Draw();

    c12->cd(1);
    // Particle smallest distance from the primary vertex in xy-plane
    evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
                values[i] = event->get_particle(0, 0, i)->dxy;
            return values;
        },
        200, -0.2, 0.2, "Particle dxy", true);

    c12->cd(2);
    // Particle dxy / dxy error
    evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
                values[i] = event->get_particle(0, 0, i)->dxy / event->get_particle(0, 0, i)->dxyErr;
            return values;
        },
        200, -5, 5, "Particle dx / dxyErr", true);
    
    c12->cd(3);
    // Particle pt error / pt
    evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 4; ++i)
                values[i] = event->get_particle(0, 0, i)->ptErr / event->get_particle(0, 0, i)->pt;
            return values;
        },
        200, 0, 0.1, "Particle ptErr / pt", true);

    c12->cd(4);

/*    evc.create_2Dhistogram(
        [](Event *event) {
            double px = 0;
            for (int i = 0; i < 2; ++i)
                px += event->get_proton(i)->px;
            return std::vector<double>{px};
        }, [](Event *event) {
            double py = 0;
            for (int i = 0; i < 2; ++i)
                py += event->get_proton(i)->py;
            return std::vector<double>{py};
        }, "Total proton px vs. py", true);*/

    c12->SaveAs((filename + "_dataB.pdf").c_str());

    std::cout << "Finished analyzing data." << std::endl;
}

void analyze_reco1(EventCollector& evc, std::string filename, std::string type) {
    std::cout << "Analyzing the first iteration of recreated particles." << std::endl;
    TFile *results = TFile::Open(evc.results.c_str(), "");

    TCanvas *c21 = new TCanvas("c21", "c21");
    c21->Draw();

    c21->cd(1);

    float min = 0.8;
    float max = 1.6;
    if (type == "pion" ) {
        min = 0.2;
        max = 1;
    }

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
        160, min, max, "Mass of recreated particles, assumed " + type + "s",
        true);

    TF1* f1 = new TF1("CauchyFit", CauchyDist, -15, 15, 3);
    f1->SetParameters(0.1, 0.7, 70);
    f1->SetParNames("Sigma", "Mean", "Scale");

    h21->Fit("CauchyFit", "", "", 0.745, 0.78);
    
    c21->SaveAs((filename + "_reco1.pdf").c_str());

    std::cout << "Finished analyzing the first iteration of recreated particles." << std::endl;
}

void analyze_reco2(EventCollector& evc, std::string filename) {
    std::cout << "Analyzing the second iteration of recreated particles." << std::endl;
    TFile *results = TFile::Open(evc.results.c_str(), "");

    TCanvas *c31 = new TCanvas("c31", "c31");
    c31->Draw();

    c31->cd(1);
    TH1 *h31 = evc.create_1Dhistogram(
        [](Event *event) {
            std::vector<double> values(4);
            for (int i = 0; i < 2; ++i) {
                values[i] = event->get_particle(2, i, 0)->mass;
            }
            return values;
        },
        100, 1.5, 2.5, "Mass of the recreated particle",
        true);
    
    c31->SaveAs((filename + "_reco2.pdf").c_str());

    std::cout << "Finished analyzing the second iteration of recreated particles." << std::endl;
}

int main()
{
    const std::string part_type = "pion";
    EventCollector evc(
//             "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/ntuples/data/TOTEM2*.root?#tree"
               "/eos/user/y/yelberke/TOTEM_2018_ADDEDVARS_OUT/minimal/TOTEM*.root?#tree"
               ,"/afs/cern.ch/user/p/ptuomola/private/particle_reconstruction_results.root");

    initialize(evc, part_type);
    filter(evc);
    analyze_data(evc, "histogram1");
    reconstruct(evc);
    analyze_reco1(evc, "histogram2", part_type);
/*    reconstruct(evc);
    analyze_reco2(evc, "histogram3");
*/
    return 0;
}